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

        self.trk_eta_max = config.get("trk_eta_max", 0.9) # maximum eta of particles
        self.trk_pT_min = config.get("trk_pT_min", 0.15) # minimum track pT for jet finding in GeV
        self.trk_pT_max = config.get("trk_pT_max", 100) # maximum track pT for jet finding in GeV

        self.jetR = config.get('jetR', 0.4) # jet radius parameter
        self.jet_pT_min = config.get("jet_pT_min", 5.0) # minimum jet pT in GeV

        self.thr = config.get("thr", 1.0) # minimum track pT for particle pairs in GeV
        self.RL_min, self.RL_max, self.RL_nbins = config["RL_binning"]
        self.pT_min, self.pT_max, self.pT_nbins = config["pT_binning"]
        self.do_hist_scaling = config.get('do_hist_scaling', True)

        self.reject_tail = config['reject_tail'] if 'reject_tail' in config else False

    def prepare(self):
        self.RL_bins = logbins(self.RL_min, self.RL_max, self.RL_nbins)
        self.pT_bins = linbins(self.pT_min, self.pT_max, self.pT_nbins)

        self.hists = {}
        self.hists['evw'] = ROOT.TH1F("evw", "evw;evw;cts", 100, -3.0, 3.0)
        self.hists['nev'] = ROOT.TH1F("hnev", "nev", 2, -0.5, 1.5)

        for ptype in ["T", "Q", "P", "M", "PM"]:
            self.hists[ptype] = ROOT.TH2D(
                f"E2C_{ptype}_{self.thr}",
                f"E2C_{ptype}_{self.thr};jet pt;R_{{L}}",
                self.pT_nbins,
                self.pT_bins,
                self.RL_nbins,
                self.RL_bins,
            )
        self.hists["jet_pT"] = ROOT.TH1D("jet_pT", "jet pt;pt gev;cts", self.pT_nbins, self.pT_bins)
        if self.do_hist_scaling:
            self.hists_to_scale = ["T", "Q", "P", "M", "PM"] + ['jet_pT']

        self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jetR)

        self.trk_selector = fj.SelectorPtMin(self.trk_pT_min) * fj.SelectorPtMax(self.trk_pT_max) * fj.SelectorAbsEtaMax(self.trk_eta_max)
        self.jet_selector = fj.SelectorPtMin(self.jet_pT_min) * fj.SelectorAbsEtaMax(self.trk_eta_max - self.jetR)

        self.thr_selector = fj.SelectorPtMin(self.thr)

        self.user_config = ['Init:showChangedParticleData = off']

        if self.reject_tail:
            if self.reject_tail is True: # check specifically for True value, not just truthiness
                self.reject_tail = 4.5 # default to JE criteria
                logger.info("Unphysical tail rejection turned ON.")
                logger.info("Rejection factor unspecified, set to 4.5.")
            elif isinstance(self.reject_tail, (float, int)):
                logger.info("Unphysical tail rejection turned ON.")
                logger.info(f"Rejection factor set to {self.reject_tail}.")
            else:
                logger.warning(f"Value '{self.reject_tail}' must be bool, float, or int; setting tail rejection OFF.")
                self.reject_tail = False
        else:
            self.reject_tail = False
            logger.info("Unphysical tail rejection turned OFF.")

        if self.do_hist_scaling:
            logger.info("Histogram scaling by cross-section ON.")
        else:
            logger.info("Histogram scaling by cross-section OFF.")

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

            parts_ch = vector[fj.PseudoJet](
                [part_w_charge(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            )
            self.evw = pythia.info.weight()
            self.analyze_event(parts_ch, pythia.info.pTHat())

            # Some "accepted" events don't survive hadronization step -- keep track here
            self.hists['nev'].Fill(0, self.evw)
            self.hists['evw'].Fill(self.evw)

            iev += 1
            if iev % ping == 0:
                logger.log(15, f"Completed {iev} events.")

    def analyze_event(self, parts_ch, pthat):
        jets = fj.sorted_by_pt(self.jet_selector(self.jet_def(self.trk_selector(parts_ch))))
        if self.reject_tail:
            jets = [jet for jet in jets if jet.pt() < self.reject_tail * pthat]

        for jet in jets:
            self.analyze_jet(jet)

    def analyze_jet(self, jet):
        jet_pT = jet.perp()
        self.hists["jet_pT"].Fill(jet_pT, self.evw)

        constituents = self.thr_selector(jet.constituents())

        for p1, p2 in permutations(constituents, 2):
            q1 = p1.user_index()
            q2 = p2.user_index()
            w = p1.perp() * p2.perp() / jet_pT / jet_pT
            rL = np.sqrt(p1.delta_phi_to(p2) ** 2 + (p1.eta() - p2.eta()) ** 2)
            if q1 > 0 and q2 > 0:
                ptype = "P"
            elif q1 < 0 and q2 < 0:
                ptype = "M"
            else:
                ptype = "PM"

            self.hists[ptype].Fill(jet_pT, rL, w * self.evw)
            self.hists["T"].Fill(jet_pT, rL, w * self.evw)
            self.hists["Q"].Fill(jet_pT, rL, w * q1 * q2 * self.evw)

    def finalize(self, pythia):
        pythia.stat()
        self.hists['nev'].SetBinError(1, 0)

        # scale_f = pythia.info.sigmaGen() / self.hists['nev'].GetBinContent(1)
        # Divide by the sum of weights to normalize histograms to a per event basis
        # then multiply by the cross-section to normalize by cross-section
        scale_f = pythia.info.sigmaGen() / pythia.info.weightSum()

        logger.info(f"sigmaGen (Pythia) = {pythia.info.sigmaGen()} +- {pythia.info.sigmaErr()}")
        logger.info(f"Scaling factor (sigma/nev) = {scale_f}")
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

        if self.do_hist_scaling:
            self._scale_hists(scale_f)
        self._save_hists_to_file()

    def _scale_hists(self, scale_f):
        for name in self.hists_to_scale:
            self.hists[name].Scale(scale_f)
        logger.info("Histograms scaled.")

    def _save_hists_to_file(self):
        with ROOT.TFile(self.output_path, "RECREATE") as file:
            self._save_hists(self.hists, file)
    def _save_hists(self, hists, file):
        for hist in hists.values():
            if isinstance(hist, dict):
                self._save_hists(hist, file)
            else:
                file.WriteTObject(hist)


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
