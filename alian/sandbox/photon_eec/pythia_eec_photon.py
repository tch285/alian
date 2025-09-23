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

def deltaR(p1, p2):
    return np.sqrt(p1.delta_phi_to(p2) ** 2 + (p1.eta() - p2.eta()) ** 2)

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
        self.parse_config()

    def parse_config(self):
        with open(self.config, "r") as stream:
            config = yaml.safe_load(stream)

        self.trk_eta_max = config.get("trk_eta_max", 0.9) # maximum eta of particles
        self.trk_pT_min = config.get("trk_pT_min", 0.15) # minimum track pT for jet finding in GeV
        self.trk_pT_max = config.get("trk_pT_max", 100) # maximum track pT for jet finding in GeV

        self.photon_eta_max = config.get("photon_eta_max", 0.7) # maximum eta of photons (match EMCal)
        self.photon_pT_min = config.get("photon_pT_min", 15) # minimum photon pT in GeV

        self.jetR = config.get('jetR', 0.4) # jet radius parameter
        self.jet_pT_min = config.get("jet_pT_min", 5.0) # minimum jet pT in GeV

        self.isocone_R = config.get('isocone_R', self.jetR) # radius of isocone
        self.isocone_pT_max = config.get('isocone_pT_max', 2) # maximum pT in isocone for iso photons
        self.isocone_pT_density_max = self.isocone_pT_max / (np.pi * self.isocone_R * self.isocone_R) # maximum isocone pT density
        self.photon_jet_dphi_min = config.get('photon_jet_dphi_min', 0.75) * np.pi
                # minimum delta phi (in units of pi) between jets and photons to ensure back-to-back

        self.thr = config.get("thr", 1.0) # minimum track pT for particle pairs in GeV
        self.RL_min, self.RL_max, self.RL_nbins = config["RL_binning"]
        self.pT_min, self.pT_max, self.pT_nbins = config["pT_binning"]
        self.photon_pT_bin_min, self.photon_pT_bin_max, self.photon_pT_bin_nbins = config["photon_pT_binning"]
        self.max_match_distance = config.get('jet_matching_fraction', 0.6) * self.jetR
        self.do_hist_scaling = config.get('do_hist_scaling', True)

    def prepare(self):
        self.RL_bins = logbins(self.RL_min, self.RL_max, self.RL_nbins)
        self.pT_bins = linbins(self.pT_min, self.pT_max, self.pT_nbins)
        self.photon_pT_bins = linbins(self.photon_pT_bin_min, self.photon_pT_bin_max, self.photon_pT_bin_nbins)

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
        self.hists["photon_pT"] = ROOT.TH1D("photon_pT", "photon pt;pt gev;cts", self.photon_pT_bin_nbins, self.photon_pT_bins)
        self.hists_to_scale = ["T", "Q", "P", "M", "PM"] + ['jet_pT'] + ['photon_pT']

        self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jetR)

        self.trk_selector = fj.SelectorPtMin(self.trk_pT_min) * fj.SelectorPtMax(self.trk_pT_max) * fj.SelectorAbsEtaMax(self.trk_eta_max)
        self.jet_selector = fj.SelectorPtMin(self.jet_pT_min) * fj.SelectorAbsEtaMax(self.trk_eta_max - self.jetR)
        self.photon_selector = fj.SelectorPtMin(self.photon_pT_min) * fj.SelectorAbsEtaMax(self.photon_eta_max)

        self.thr_selector = fj.SelectorPtMin(self.thr)

        self.user_config = ['Init:showChangedParticleData = off']
        # self.user_config.append("HadronLevel:all=off") # need this to get initial partons

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

            parts_ch = vector[fj.PseudoJet](
                [part_w_charge(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            )
            photons = vector[fj.PseudoJet](
                [part_w_charge(p) for p in pythia.event if p.isFinal() and p.id() == 22]
            )

            self.evw = pythia.info.weight()
            self.analyze_event(parts_ch, photons)

            # initial photon is status -23, index 6, id 22

            self.hists['nev'].Fill(0, self.evw)
            self.hists['evw'].Fill(self.evw)

            iev += 1

    def analyze_event(self, parts_ch, photons):
        photons_acc = self.photon_selector(photons)
        iso_photons = self.get_iso_photons(parts_ch, photons_acc) # all isolated photons within acceptance
        if not iso_photons:
            return
        # if at least one isolated photon exists
        for photon in iso_photons:
            self.analyze_photon(photon)
        jets = fj.sorted_by_pt(self.jet_selector(self.jet_def(self.trk_selector(parts_ch))))
        jets_b2b = self.get_b2b_jets(jets, iso_photons)
        for jet in jets_b2b:
            self.analyze_jet(jet)


    def get_iso_photons(self, parts_ch, photons):
        """Return set of isolated photons from all photons and charged particles within acceptance."""
        return [ph for ph in photons if self.get_isocone_pT(ph, parts_ch) / self.get_isocone_area(ph) < self.isocone_pT_density_max]
    def get_isocone_pT(self, photon, parts_ch):
        """For a given photon and set of charged particles, get the total pT of the isolation cone."""
        return np.sum([part_ch.perp() for part_ch in parts_ch if deltaR(part_ch, photon) < self.isocone_R])
    def get_isocone_area(self, ph):
        """Get area of photon isolation cone, taking edge effects into account."""
        eta = np.abs(ph.eta())
        if eta < (self.trk_eta_max - self.isocone_R):
            return np.pi * self.isocone_R * self.isocone_R
        else:
            return (self.trk_eta_max - eta) * np.sqrt(self.isocone_R**2 - (self.trk_eta_max - eta)**2 ) + np.pi * self.isocone_R**2 * (1 - 1 / np.pi * np.arccos((self.trk_eta_max - eta) / self.isocone_R))
    def get_b2b_jets(self, jets, iso_photons):
        """Get all jets that are back-to-back with at least one photon"""
        return [jet for jet in jets if self.is_b2b(jet, iso_photons)]
    def is_b2b(self, jet, iso_photons):
        for iso_ph in iso_photons:
            if np.abs(jet.delta_phi_to(iso_ph)) > self.photon_jet_dphi_min:
                return True
        return False

    def store_leading_partons(self, event):
        leading_parton1 = fj.PseudoJet(event[5].px(),event[5].py(),event[5].pz(),event[5].e())
        leading_parton1.set_user_index(event[5].idAbs())
        leading_parton2 = fj.PseudoJet(event[6].px(),event[6].py(),event[6].pz(),event[6].e())
        leading_parton2.set_user_index(event[6].idAbs())
        self.parton_parents = [leading_parton1, leading_parton2]
        assert event[5].status() == event[6].status() == -23
        # event.list()
        # print([event[i].idAbs() for i in range(10)])
    def match_jet_to_partons(self, jet):
        candidates = [part for part in self.parton_parents if deltaR(jet, part) < self.max_match_distance]
        if not candidates:
            return "0m"
        elif len(candidates) == 1:
            return candidates[0].user_index()
        else:
            return "2m"

    def analyze_photon(self, photon):
        self.hists['photon_pT'].Fill(photon.pt())
    def analyze_jet(self, jet):
        self.hists['jet_pT'].Fill(jet.pt())
        jet_pT = jet.perp()

        constituents = self.thr_selector(jet.constituents())

        for p1, p2 in permutations(constituents, 2):
            q1 = p1.user_index()
            q2 = p2.user_index()
            w = p1.perp() * p2.perp() / jet_pT / jet_pT
            rL = deltaR(p1, p2)
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
        scale_f = pythia.info.sigmaGen() / pythia.info.weightSum()

        logger.info(f"sigmaGen (Pythia) = {pythia.info.sigmaGen()} +- {pythia.info.sigmaErr()}")
        logger.info(f"Scaling factor (sigma/nev) = {scale_f}")
        logger.info(f"nTried: {pythia.info.nTried()}")
        logger.info(f"nSelected: {pythia.info.nSelected()}")
        logger.info(f"nAccepted: {pythia.info.nAccepted()}")
        logger.info(f"Failed (generation): {self.nfailed}")
        logger.info(f"Sum of weights (histogram): {self.hists['nev'].GetBinContent(1)}")
        logger.info(f"Sum of weights (Pythia weightSum): {pythia.info.weightSum()}")
        logger.info(f"Number of successful events (histogram): {self.hists['nev'].GetEntries()}")
        # with open(f"{self.output_dir}/scales.txt", 'w') as f:
        #     f.write(f"{scale_f}")
        if self.do_hist_scaling:
            self._scale_hists(scale_f)
        self._save_hists_to_file()
        logger.info(f"Histograms saved to: {self.output_path}")

        logger.info(
            f"N total final events: {self.hists['nev'].GetBinContent(1)} with {pythia.info.nAccepted() - self.hists['nev'].GetBinContent(1)} events rejected."
        )

    def _scale_hists(self, scale_f):
        for name in self.hists_to_scale:
            self.hists[name].Scale(scale_f)

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
        description="Pythia 8 + FJ for photon-jet EEC", prog=os.path.basename(__file__)
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
    parser.add_argument("-v", "--verbose", help="Print event generation information", action="store_true")
    args = parser.parse_args()

    level = 20 # INFO level by default

    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(funcName)s - %(message)s"
    )
    handler.setFormatter(formatter)
    handler.setLevel(level)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    proc = PythiaOTFENC(
        args.config, args.output_dir, args.output_filename, args.nev, args.verbose
    )
    proc.generate(args)

if __name__ == "__main__":
    main()
