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
    def __init__(self, config, odir, ofn, nev):
        self.config = config
        self.output_path = f"{odir}/{ofn}"
        if nev < 10:
            logger.warning(f"Number of events ({nev}) too low, setting to 10.")
            self.nev = 10
        else:
            self.nev = nev
        self.nfailed = 0
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
        self.ew_min, self.ew_max, self.ew_nbins = config["ew_binning"]
        self.shower = config["shower"] if 'shower' in config else "simple"
        if self.shower not in ["dire", "vincia", "simple"]:
            logger.warning(
                f"Shower model {self.shower} not recognized, defaulting to simple."
            )
            self.shower = "simple"
        self.custom_tune = config["custom_tune"] if 'custom_tune' in config else False
        self.resonance_decay = config['resonance_decay'] if 'resonance_decay' in config else True
        self.strange_decay = config['strange_decay'] if 'strange_decay' in config else True
        self.shuffle_pt = config['shuffle_pt'] if 'shuffle_pt' in config else True
        self.shuffle_q = config['shuffle_q'] if 'shuffle_q' in config else True
        self.reject_tail = config['reject_tail'] if 'reject_tail' in config else False
        self.scale_by_xsec = config['scale_by_xsec'] if 'scale_by_xsec' in config else True

    def prepare(self):
        self.RL_bins = logbins(self.RL_min, self.RL_max, self.RL_nbins)
        self.pT_bins = linbins(self.pT_min, self.pT_max, self.pT_nbins)
        self.ew_bins = logbins(self.ew_min, self.ew_max, self.ew_nbins)

        self.hists = {}
        self.hists_to_scale = []
        self.hists['evw'] = ROOT.TH1F("evw", "evw;evw;cts", 100, -3.0, 3.0)
        self.hists['nev'] = ROOT.TH1F("hnev", "nev", 2, -0.5, 1.5)
        self.hists_to_scale += ['evw']

        self.hists["ew"] = ROOT.TH1D("ew", "ew;ew;cts", self.ew_nbins, self.ew_bins)
        self.hists["q"] = ROOT.TH1D("q", "q;q;cts", 3, -1.5, 1.5)
        for ptype in ["T", "Q", "P", "M", "PM"]:
            self.hists[ptype] = ROOT.TH2F(
                f"E2C_{ptype}_{self.thr}",
                f"E2C_{ptype}_{self.thr};jet pt;R_{{L}}",
                self.pT_nbins,
                self.pT_bins,
                self.RL_nbins,
                self.RL_bins,
            )
        self.hists_to_scale += ['ew', 'q', "T", "Q", "P", "M", "PM"]

        self.hists["hg_ew"] = ROOT.TH1D("hg_ew", "hg_ew;ew;cts", self.ew_nbins, self.ew_bins)
        self.hists["hg_q"] = ROOT.TH1D("hg_q", "hg_q;q;cts", 3, -1.5, 1.5)
        for ptype in ["hg_T", "hg_Q", "hg_P", "hg_M", "hg_PM"]:
            self.hists[ptype] = ROOT.TH2F(
                f"E2C_{ptype}_{self.thr}",
                f"E2C_{ptype}_{self.thr};jet pt;R_{{L}}",
                self.pT_nbins,
                self.pT_bins,
                self.RL_nbins,
                self.RL_bins,
            )
        self.hists_to_scale += ['hg_ew', 'hg_q', "hg_T", "hg_Q", "hg_P", "hg_M", "hg_PM"]
        self.hists["jet_pT"] = ROOT.TH1D(
            "jet_pT", "jet pt;pt gev;cts", self.pT_nbins, self.pT_bins
        )
        self.hists_to_scale += ['jet_pT']

        self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jetR)

        self.jet_selector = fj.SelectorPtMin(self.min_jet_pT) * fj.SelectorAbsEtaMax(
            self.eta_max - self.jetR
        )

        self.part_pT_selector = fj.SelectorPtMin(self.min_part_trk_pT)
        self.thr_selector = fj.SelectorPtMin(self.thr)

        self.user_config = ['Init:showProcesses = off', 'Init:showChangedParticleData = off']
        if self.shower == "vincia":
            logger.info("Shower model: Vincia")
            self.user_config.append("PartonShowers:Model = 2")
            if self.custom_tune:
                logger.info("Shower tune: using Vincia tune")
            else:
                logger.info("Shower tune: using default (non-Vincia) tune")
                self.user_config.append("Vincia:Tune = -1")
        elif self.shower == "dire":
            logger.info("Shower model: Dire")
            self.user_config.append("PartonShowers:Model = 3")
            if self.custom_tune:
                logger.info("Shower tune: using Dire tune")
            else:
                logger.info("Shower tune: using default (non-Dire) tune")
                self.user_config.append("Dire:Tune = -1")
        else:
            logger.info("Shower model: simple")
            logger.info("Shower tune: using Monash tune")

        if not self.resonance_decay:
            self.user_config.append("HadronLevel:Decay = off")
            logger.info("Hadronic resonance decays turned OFF.")
        else:
            # if False, do nothing since HadronLevel:Decay is on by default
            logger.info("Hadronic resonance decays turned ON.")

        if not self.strange_decay:
            pdg_codes = [
                310,# K0s
                3122,# Lambda0
                3112,# Sigma-
                3212,# Sigma0
                3222,# Sigma+
                3312,# Xi-
                3322,# Xi+
                3334,# Omega-
                1114,# Delta-
                2114,# Delta0
                2214,# Delta+
            ]
            for pdg_code in pdg_codes:
                self.user_config.append(f"{pdg_code}:mayDecay = off")
            logger.info("Strange decays turned OFF.")
        else:
            logger.info("Strange decays turned ON.")

        if self.reject_tail:
            if self.reject_tail is True: # check specifically for True value, not just truthiness
                self.reject_tail = 4.5 # default to JE criteria
                logger.info("Unphysical tail rejection turned ON.")
                logger.info("Rejection factor unspecified, set to 4.5.")
            elif isinstance(self.reject_tail, (float, int)):
                logger.info("Unphysical tail rejection turned ON.")
                logger.info(f"Rejection factor set to {self.reject_tail}.")
            else:
                logger.warning(f"Value '{self.reject_tail}' must be bool, float, or int; turning tail rejection OFF.")
                self.reject_tail = False
        else:
            self.reject_tail = False
            logger.info("Unphysical tail rejection turned OFF.")

        if self.scale_by_xsec:
            logger.info("Histogram scaling by cross-section ON.")
        else:
            logger.info("Histogram scaling by cross-section OFF.")


    def generate(self, args):
        self.prepare()

        # pythia = Pythia8.Pythia()
        fj.ClusterSequence.print_banner()

        self.start = time.perf_counter()
        pythia = pyconf.create_and_init_pythia_from_args(args, self.user_config)
        logger.info(f"Initializing RNG with seed {args.py_seed}.")
        self.rng = np.random.default_rng(args.py_seed)

        if not pythia:
            logger.fatal("Pythia initialization failed!")
            sys.exit(1)
        else:
            logger.info("Pythia configured and initialized.")

        self.simulate(pythia)
        logger.info(f"Simulation step completed in {time.perf_counter() - self.start:.2f} sec.")

        self.finalize(pythia)
        logger.info(f"Simulation completed in {time.perf_counter() - self.start:.2f} sec.")

    def simulate(self, pythia):
        iev = 0
        # ping every 10% of total events
        ping = self.nev // 10
        while iev < self.nev:
            if not pythia.next():
                self.nfailed += 1
                continue

            # parts = vector[fj.PseudoJet](
            #     [
            #         fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
            #         for p in pythia.event
            #         if p.isFinal() and p.isCharged()
            #     ]
            # )
            parts_ch = vector[fj.PseudoJet](
                [part_w_charge(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            )
            self.evw = pythia.info.weight()
            self.analyze_event(parts_ch, pythia.info.pTHat())

            # Some "accepted" events don't survive hadronization step -- keep track here
            self.hists['nev'].Fill(0, self.evw)
            self.hists['evw'].Fill(self.evw)

            # logger.info('NEW EVENT:')
            # logger.info(f"Sum of weights (histogram): {self.hists['nev'].GetBinContent(1)}")
            # logger.info(f"Sum of weights (Pythia): {pythia.info.weightSum()}")
            iev += 1
            # if iev == 10:
            #     break
            if iev % ping == 0:
                logger.log(15, f"Completed {iev} events.")

    def analyze_event(self, parts, pthat):
        jets = fj.sorted_by_pt(
            self.jet_selector(self.jet_def(self.part_pT_selector(parts)))
        )
        if self.reject_tail:
            jets = [jet for jet in jets if jet.pt() < self.reject_tail * pthat]

        for jet in jets:
            self.analyze_jet(jet)

    def analyze_jet(self, jet):
        jet_pT = jet.perp()
        self.hists["jet_pT"].Fill(jet_pT, self.evw)

        constituents = self.thr_selector(jet.constituents())
        pTs = np.array([p.perp() for p in constituents])
        qs = np.array([p.user_index() for p in constituents])
        nconst = len(constituents)
        if self.shuffle_pt:
            self.rng.shuffle(pTs)
        if self.shuffle_q:
            self.rng.shuffle(qs)
        for p1, p2 in permutations(constituents, 2):
            q1 = p1.user_index()
            q2 = p2.user_index()
            w = p1.perp() * p2.perp() / jet_pT / jet_pT
            rL = np.sqrt(p1.delta_phi_to(p2) ** 2 + (p1.eta() - p2.eta()) ** 2)
            if q1 > 0 and q2 > 0:
                ptype = "P"
                ptypeq = 1
            elif q1 < 0 and q2 < 0:
                ptype = "M"
                ptypeq = -1
            else:
                ptype = "PM"
                ptypeq = 0

            self.hists[ptype].Fill(jet_pT, rL, w * self.evw)
            self.hists["T"].Fill(jet_pT, rL, w * self.evw)
            self.hists["Q"].Fill(jet_pT, rL, w * q1 * q2 * self.evw)
            self.hists["ew"].Fill(w, self.evw)
            self.hists["q"].Fill(ptypeq, self.evw)
        for i1, i2 in permutations(range(nconst), 2):
            p1 = constituents[i1]
            p2 = constituents[i2]
            q1, q2 = qs[i1], qs[i2]
            w = pTs[i1] * pTs[i2] / jet_pT / jet_pT
            rL = np.sqrt(p1.delta_phi_to(p2) ** 2 + (p1.eta() - p2.eta()) ** 2)
            if q1 > 0 and q2 > 0:
                ptype = "hg_P"
                ptypeq = 1
            elif q1 < 0 and q2 < 0:
                ptype = "hg_M"
                ptypeq = -1
            else:
                ptype = "hg_PM"
                ptypeq = 0

            self.hists[ptype].Fill(jet_pT, rL, w * self.evw)
            self.hists["hg_T"].Fill(jet_pT, rL, w * self.evw)
            self.hists["hg_Q"].Fill(jet_pT, rL, w * q1 * q2 * self.evw)
            self.hists["hg_ew"].Fill(w, self.evw)
            self.hists["hg_q"].Fill(ptypeq, self.evw)

    def finalize(self, pythia):
        pythia.stat()
        self.hists['nev'].SetBinError(1, 0)

        # scale_f = pythia.info.sigmaGen() / self.hists['nev'].GetBinContent(1)
        # Divide by the sum of weights to normalize histograms to a per event basis
        # then multiply by the cross-section to normalize by cross-section
        scale_f = pythia.info.sigmaGen() / pythia.info.weightSum()

        logger.info(f"sigmaGen (Pythia) = {pythia.info.sigmaGen()} +- {pythia.info.sigmaErr()}")
        logger.info(f"Scaling factor (sigma/nev) = {scale_f} (with weightSum)")
        logger.info(f"nTried: {pythia.info.nTried()}")
        logger.info(f"nSelected: {pythia.info.nSelected()}")
        logger.info(f"nAccepted: {pythia.info.nAccepted()}")
        logger.info(f"Failed: {self.nfailed}")
        logger.info(f"Sum of weights (histogram): {self.hists['nev'].GetBinContent(1)}")
        logger.info(f"Sum of weights (Pythia): {pythia.info.weightSum()}")
        logger.info(f"Number of successful events (histogram): {self.hists['nev'].GetEntries()}")
        # with open(f"{self.output_dir}/scales.txt", 'w') as f:
        #     f.write(f"{scale_f}")
        logger.info(
            f"N total final events: {self.hists['nev'].GetBinContent(1)} with {pythia.info.nAccepted() - self.hists['nev'].GetBinContent(1)} events rejected."
        )

        if self.scale_by_xsec:
            self._scale_hists(scale_f)
        self._save_hists()

    def _scale_hists(self, scale_f):
        # for ptype in ["T", "Q", "P", "M", "PM"]:
        logger.info(f"Histograms to scale: {self.hists_to_scale}")
        for key in self.hists_to_scale:
            self.hists[key].Scale(scale_f)
        # self.hists["jet_pT"].Scale(scale_f)
        logger.info(f"Histograms scaled.")

    def _save_hists(self):
        with ROOT.TFile(self.output_path, "RECREATE") as f:
            for hist in self.hists.values():
                f.WriteTObject(hist)
        logger.info(f"Histograms saved to {self.output_path}.")


def main():
    parser = argparse.ArgumentParser(
        description="pythia8 OTF + FJ for cEEC", prog=os.path.basename(__file__)
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
    parser.add_argument("-t", "--track", help="Print event generation information", action="store_true")
    args = parser.parse_args()

    level = 20 # INFO level by default
    if args.track:
        level = 15
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(funcName)s - %(message)s"))
    handler.setLevel(level)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    proc = PythiaOTFENC(
        args.config, args.output_dir, args.output_filename, args.nev
    )
    proc.generate(args)


if __name__ == "__main__":
    main()
