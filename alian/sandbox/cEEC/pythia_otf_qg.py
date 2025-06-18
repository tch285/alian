#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import time
from itertools import permutations

import numpy as np
import ROOT
import yaml

import heppyy.util.fastjet_cppyy # noqa: F401
# import heppyy.util.pythia8_cppyy
# import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl.std import vector
from heppyy.pythia_util import configuration as pyconf

# to import more from heppyy: import the NAMESPACE declared in the .hh file
# ex: for pythiafjext/pyfjtools.hh, we import:
# from cppyy.gbl import pythiafjtools

logger = logging.getLogger(__name__)

def linbins(xmin, xmax, nbins):
    return np.linspace(xmin, xmax, nbins + 1)

def logbins(xmin, xmax, nbins):
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins + 1)

def part_w_charge(pyp):
    pj = fj.PseudoJet(pyp.px(), pyp.py(), pyp.pz(), pyp.e())
    pj.set_user_index(int(pyp.charge()))
    return pj

class PythiaOTFENC(object):
    def __init__(self, config, odir, ofn, nev, verbose):
        self.config = config
        self.verbose = verbose
        self.output_path = f"{odir}/{ofn}"
        if nev < 10:
            logger.warning(f"Number of events ({nev}) too low, setting to 10.")
            self.nev = 10
        else:
            self.nev = nev
        self.nfailed = 0
        self.nfailed_hadr = 0
        self.parse_config()

    def parse_config(self):
        with open(self.config, "r") as stream:
            config = yaml.safe_load(stream)
        self.jetR = config["jetR"] if "jetR" in config else 0.4
        self.eta_max = config["eta_max"]
        self.min_part_trk_pT = config["min_part_trk_pT"]
        self.thr = config["thr"]
        self.min_jet_pT = config["min_jet_pT"]
        self.RL_min, self.RL_max, self.RL_nbins = config["RL_binning"]
        self.pT_min, self.pT_max, self.pT_nbins = config["pT_binning"]
        self.max_match_distance = config['jet_matching_fraction'] * self.jetR
        self.parton_ids = config["parton_ids"]
        self.parton_ids_ext = self.parton_ids + [0, "0m", "2m"] # no match, two matches
        
        # self.shower = config["shower"]
        # if self.shower not in ["dire", "vincia", "simple"]:
        #     logger.warning(
        #         f"Shower model {self.shower} not recognized, defaulting to simple."
        #     )
        #     self.shower = "simple"
        # self.custom_tune = config["custom_tune"]

    def prepare(self):
        self.RL_bins = logbins(self.RL_min, self.RL_max, self.RL_nbins)
        self.pT_bins = linbins(self.pT_min, self.pT_max, self.pT_nbins)

        self.hists = {}
        self.hists['evw'] = ROOT.TH1F("evw", "evw;evw;cts", 100, -3.0, 3.0)
        self.hists['nev'] = ROOT.TH1F("hnev", "nev", 2, -0.5, 1.5)
        
        for parton_id in self.parton_ids_ext:
            self.hists[parton_id] = {}
            for ptype in ["T", "Q", "P", "M", "PM"]:
                self.hists[parton_id][ptype] = ROOT.TH2D(
                    f"E2C_{ptype}_{parton_id}_{self.thr}",
                    f"E2C_{ptype}_{parton_id}_{self.thr};jet pt;R_{{L}}",
                    self.pT_nbins,
                    self.pT_bins,
                    self.RL_nbins,
                    self.RL_bins,
                )
            self.hists[parton_id]["jet_pT"] = ROOT.TH1D(
                f"jet_pT_{parton_id}", f"jet pt id {parton_id};pt gev;cts", self.pT_nbins, self.pT_bins
            )

        self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jetR)

        self.jet_selector = fj.SelectorPtMin(self.min_jet_pT) * fj.SelectorAbsEtaMax(
            self.eta_max - self.jetR
        )

        self.part_pT_selector = fj.SelectorPtMin(self.min_part_trk_pT)
        self.thr_selector = fj.SelectorPtMin(self.thr)

        self.user_config = ['Init:showProcesses = off', 'Init:showChangedParticleData = off']
        self.user_config.append("HadronLevel:all=off") # need this to get initial partons
        # if self.shower == "vincia":
        #     self.user_config.append("PartonShowers:Model = 2")
        #     if not self.custom_tune:
        #         self.user_config.append("Vincia:Tune = -1")
        # elif self.shower == "dire":
        #     self.user_config.append("PartonShowers:Model = 3")
        #     if not self.custom_tune:
        #         self.user_config.append("Dire:Tune = -1")

    def generate(self, args):
        self.prepare()

        # pythia = Pythia8.Pythia()
        fj.ClusterSequence.print_banner()

        self.start = time.perf_counter()
        pythia = pyconf.create_and_init_pythia_from_args(args, self.user_config)

        if not pythia:
            logger.fatal("Pythia initialization failed!")
            sys.exit(1)
        else:
            logger.info("Pythia configured and initialized.")

        self.simulate(pythia)
        logger.info(f"Simulation step completed in {time.perf_counter() - self.start} sec.")

        self.finalize(pythia)
        logger.info(f"Histogram scaling/saving completed in {time.perf_counter() - self.start} sec.")

    def simulate(self, pythia):
        iev = 0
        while iev < self.nev:
            if iev % 100 == 0:
                logger.debug(f"ievt: {iev}")

            if not pythia.next():
                self.nfailed += 1
                logger.warning("Generation of event failed, skipping.")
                continue

            self.store_leading_partons(pythia.event)

            hstatus = pythia.forceHadronLevel()
            if not hstatus:
                self.nfailed_hadr += 1
                logger.warning("Hadronization of event failed, skipping.")
                continue

            parts_ch = vector[fj.PseudoJet](
                [part_w_charge(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            )

            self.evw = pythia.info.weight()
            self.analyze_event(parts_ch)

            # Some "accepted" events don't survive hadronization step -- keep track here
            self.hists['nev'].Fill(0, self.evw)
            self.hists['evw'].Fill(self.evw)

            # logger.info('NEW EVENT:')
            # logger.info(f"Sum of weights (histogram): {self.hists['nev'].GetBinContent(1)}")
            # logger.info(f"Sum of weights (Pythia): {pythia.info.weightSum()}")
            # if iev == 10:
            #     break
            iev += 1

    def analyze_event(self, parts):
        jets = fj.sorted_by_pt(
            self.jet_selector(self.jet_def(self.part_pT_selector(parts)))
        )
        for jet in jets:
            self.analyze_jet(jet)

    def store_leading_partons(self, event):
        leading_parton1 = fj.PseudoJet(event[5].px(),event[5].py(),event[5].pz(),event[5].e())
        leading_parton1.set_user_index(event[5].idAbs())
        leading_parton2 = fj.PseudoJet(event[6].px(),event[6].py(),event[6].pz(),event[6].e())
        leading_parton2.set_user_index(event[6].idAbs())
        self.parton_parents = [leading_parton1, leading_parton2]
        assert event[5].status() == event[6].status() == -23
        # event.list()
        # print([event[i].idAbs() for i in range(10)])

    def deltaR(self, p1, p2):
        return np.sqrt(p1.delta_phi_to(p2) ** 2 + (p1.eta() - p2.eta()) ** 2)

    def match_jet_to_partons(self, jet):
        candidates = [part for part in self.parton_parents if self.deltaR(jet, part) < self.max_match_distance]
        if not candidates:
            return "0m"
        elif len(candidates) == 1:
            return candidates[0].user_index()
        else:
            return "2m"

    def analyze_jet(self, jet):
        parton_ids = [0, self.match_jet_to_partons(jet)] # 0 is all jets
        jet_pT = jet.perp()

        constituents = self.thr_selector(jet.constituents())

        for parton_id in parton_ids:
            self.hists[parton_id]["jet_pT"].Fill(jet_pT, self.evw)
        for p1, p2 in permutations(constituents, 2):
            q1 = p1.user_index()
            q2 = p2.user_index()
            w = p1.perp() * p2.perp() / jet_pT / jet_pT
            rL = self.deltaR(p1, p2)
            if q1 > 0 and q2 > 0:
                ptype = "P"
            elif q1 < 0 and q2 < 0:
                ptype = "M"
            else:
                ptype = "PM"
            for parton_id in parton_ids:
                self.hists[parton_id][ptype].Fill(jet_pT, rL, w * self.evw)
                self.hists[parton_id]["T"].Fill(jet_pT, rL, w * self.evw)
                self.hists[parton_id]["Q"].Fill(jet_pT, rL, w * q1 * q2 * self.evw)

    def finalize(self, pythia):
        pythia.stat()
        self.hists['nev'].SetBinError(1, 0)

        scale_f = pythia.info.sigmaGen() / self.hists['nev'].GetBinContent(1)
        # scale_f = pythia.info.sigmaGen() / pythia.info.weightSum()

        logger.info(f"sigmaGen (Pythia) = {pythia.info.sigmaGen()} +- {pythia.info.sigmaErr()}")
        logger.info(f"Scaling factor (sigma/nev) = {scale_f} (with histogram)")
        logger.info(f"nTried: {pythia.info.nTried()}")
        logger.info(f"nSelected: {pythia.info.nSelected()}")
        logger.info(f"nAccepted: {pythia.info.nAccepted()}")
        logger.info(f"Failed (generation): {self.nfailed}")
        logger.info(f"Failed (hadronization): {self.nfailed_hadr}")
        logger.info(f"Sum of weights (histogram): {self.hists['nev'].GetBinContent(1)}")
        logger.info(f"Sum of weights (Pythia): {pythia.info.weightSum()}")
        logger.info(f"Number of successful events (histogram): {self.hists['nev'].GetEntries()}")
        # with open(f"{self.output_dir}/scales.txt", 'w') as f:
        #     f.write(f"{scale_f}")
        self._scale_hists(scale_f)
        self._save_hists()

        logger.info(
            f"N total final events: {self.hists['nev'].GetBinContent(1)} with {pythia.info.nAccepted() - self.hists['nev'].GetBinContent(1)} events rejected."
        )

    def _scale_hists(self, scale_f):
        for parton_id in self.parton_ids_ext:
            for ptype in ["T", "Q", "P", "M", "PM"]:
                self.hists[parton_id][ptype].Scale(scale_f)
            self.hists[parton_id]["jet_pT"].Scale(scale_f)

    def _save_hists(self):
        with ROOT.TFile(self.output_path, "RECREATE") as f:
            f.WriteTObject(self.hists['nev'])
            f.WriteTObject(self.hists['evw'])
            for parton_id in self.parton_ids_ext:
                for ptype in ["T", "Q", "P", "M", "PM"]:
                    f.WriteTObject(self.hists[parton_id][ptype])
                f.WriteTObject(self.hists[parton_id]["jet_pT"])

def main():
    parser = argparse.ArgumentParser(
        description="pythia8 OTF + FJ for cEEC (qg studies)", prog=os.path.basename(__file__)
    )
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument("-c", "--config", help="Config file", required=True, type=str)
    parser.add_argument(
        "-od", "--output-dir", help="Output directory", required=True, type=str
    )
    parser.add_argument(
        "-of",
        "--output-filename",
        help="Output ROOT filename",
        default="AnalysisResults.root",
        type=str,
    )
    # parser.add_argument("-n", "--nev", help="Number of events", default=10, type=int)
    parser.add_argument("-v", "--verbose", help="be verbose", action="store_true")
    args = parser.parse_args()
    
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(funcName)s - %(message)s"
    )
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    proc = PythiaOTFENC(
        args.config, args.output_dir, args.output_filename, args.nev, args.verbose
    )
    proc.generate(args)


if __name__ == "__main__":
    main()
