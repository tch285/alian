#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import time
# from itertools import permutations

import numpy as np
import ROOT
import yaml

import heppyy.util.fastjet_cppyy # noqa: F401
import heppyy.util.pythia8_cppyy # noqa: F401
import heppyy.util.heppyy_cppyy # noqa: F401

from cppyy.gbl import fastjet as fj
from cppyy.gbl.std import vector
from heppyy.pythia_util import configuration as pyconf
from cppyy.gbl import pythiafjtools

logger = logging.getLogger(__name__)

def linbins(xmin, xmax, nbins):
    return np.linspace(xmin, xmax, nbins + 1)

def logbins(xmin, xmax, nbins):
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins + 1)

def get_mother_status_list(pythia, p8p, mother_status_list):
    mother_idxs = p8p.motherList()
    if(len(mother_idxs) > 0):
        for imother in mother_idxs:
            get_mother_status_list(pythia, pythia.event[imother], mother_status_list)

    mother_status_list.append( p8p.status() )
    return mother_status_list

def get_src(pythia, p8p):
    mother_status_list = np.array(get_mother_status_list(pythia, p8p, []))
    types = (mother_status_list/10).astype(int)
    # 2 hard process
    # 3 UE
    # 4 ISR
    # 5 parton shower
    # 6 beam remnant
    if -3 in types:
        return 3
    elif -2 in types:
        return 2
    elif -4 in types:
        return 4
    elif -6 in types:
        return 6
    else: # catchall for missing types
        return 0

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
        # self.thr = config["thr"]
        self.min_jet_pT = config["min_jet_pT"]
        self.part_pT_min, self.part_pT_max, self.part_pT_nbins = config["part_pT_binning"]
        self.jet_pT_min, self.jet_pT_max, self.jet_pT_nbins = config["jet_pT_binning"]
        self.shower = config["shower"]
        if self.shower not in ["dire", "vincia", "simple"]:
            logger.warning(
                f"Shower model {self.shower} not recognized, defaulting to simple."
            )
            self.shower = "simple"
        self.custom_tune = config.get("custom_tune", False)
        self.resonance_decay = config.get('resonance_decay', True)
        self.strange_decay = config.get('strange_decay', True)
        self.strong_decay = config.get('strong_decay', True)
        self.Kstarphi_decay = config.get('Kstarphi_decay', True)
        self.reject_tail = config.get('reject_tail', False)
        self.scale_by_xsec = config.get('scale_by_xsec', True)

    def prepare(self):
        self.part_pT_bins = logbins(self.part_pT_min, self.part_pT_max, self.part_pT_nbins)
        self.jet_pT_bins = linbins(self.jet_pT_min, self.jet_pT_max, self.jet_pT_nbins)

        self.hists = {}
        self.hists['evw'] = ROOT.TH1F("evw", "evw;evw;cts", 100, -3.0, 3.0)
        self.hists['nev'] = ROOT.TH1F("hnev", "nev", 2, -0.5, 1.5)

        # 2 hard process
        # 3 UE
        # 4 ISR
        # 5 parton shower
        # 6 beam remnant
        for src in [2, 3, 4, 6, 0]:
            self.hists[src] = {}
            for sgn in ["P", "M"]:
                self.hists[src][sgn] = ROOT.TH2D(
                    f"pt_{src}_{sgn}",
                    f"pT_{src}_{sgn};part p_{{T}};jet p_{{T}}",
                    self.part_pT_nbins,
                    self.part_pT_bins,
                    self.jet_pT_nbins,
                    self.jet_pT_bins,
                )

        self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jetR)

        self.jet_selector = fj.SelectorPtMin(self.min_jet_pT) * fj.SelectorAbsEtaMax(
            self.eta_max - self.jetR
        )

        self.part_pT_selector = fj.SelectorPtMin(self.min_part_trk_pT)
        # self.thr_selector = fj.SelectorPtMin(self.thr)

        self.user_config = ['Init:showProcesses = off', 'Init:showChangedParticleData = on']
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

        if not self.Kstarphi_decay:
            pdg_codes = [
                313,  # K*(892)
                333,  # phi(1020)
            ]
            for pdg_code in pdg_codes:
                self.user_config.append(f"{pdg_code}:mayDecay = off")
            logger.info("Kstar and phi decays turned OFF.")
        else:
            logger.info("Kstar and phi decays turned ON.")

        if not self.strong_decay:
            pdg_codes = [
                223,  # omega
                213,  # rho+
                -213, # rho-
                113,  # rho0
                333,  # phi
                221,  # eta
                331,  # eta'
                1114, # Delta-
                2114, # Delta0
                2214, # Delta+
            ]
            for pdg_code in pdg_codes:
                self.user_config.append(f"{pdg_code}:mayDecay = off")
            logger.info("Strong decays turned OFF.")
        else:
            logger.info("Strong decays turned ON.")

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
        self.pythia = pyconf.create_and_init_pythia_from_args(args, self.user_config)

        if not self.pythia:
            logger.fatal("Pythia initialization failed!")
            sys.exit(1)
        else:
            logger.info("Pythia configured and initialized.")

        self.simulate()
        logger.info(f"Simulation step completed in {time.perf_counter() - self.start:.2f} sec.")

        self.finalize()
        logger.info(f"Simulation completed in {time.perf_counter() - self.start:.2f} sec.")

    def simulate(self):
        iev = 0
        # ping every 10% of total events
        ping = self.nev // 10
        while iev < self.nev:
            if not self.pythia.next():
                self.nfailed += 1
                continue

            # parts = vector[fj.PseudoJet](
            #     [
            #         fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
            #         for p in pythia.event
            #         if p.isFinal() and p.isCharged()
            #     ]
            # )
            # parts_ch = vector[fj.PseudoJet](
            #     [part_w_charge(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            # )
            # parts_sc = vector[fj.PseudoJet](
            #     [part_w_status(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            # )
            # parts = vector[fj.PseudoJet](
            #     [part_w_charge(p) for p in pythia.event if p.isFinal() and p.isCharged()]
            # )
            # vector[int]
            # print(type(int(pythiafjtools.Py8Part.kCharged)))
            parts = pythiafjtools.vectorize_select(self.pythia, vector[int]([pythiafjtools.kCharged, pythiafjtools.kFinal]), 0, True, 0.13957) # assume pi+- mass
            # parts = pythiafjtools.vectorize_select(self.pythia, [int(pythiafjtools.Py8Part.kCharged), int(pythiafjtools.Py8Part.kFinal)], 0, True, 0.13957) # assume pi+- mass
            self.evw = self.pythia.info.weight()
            self.pThat = self.pythia.info.pTHat()
            self.analyze_event(parts)

            # Some "accepted" events don't survive hadronization step -- keep track here
            self.hists['nev'].Fill(0, self.evw)
            self.hists['evw'].Fill(self.evw)

            iev += 1
            if iev % ping == 0:
                logger.log(15, f"Completed {iev} events.")

    def analyze_event(self, parts):
        # [self.hists["status"].Fill(part.user_index(), self.evw) for part in parts_sc]

        jets = fj.sorted_by_pt(
            self.jet_selector(self.jet_def(self.part_pT_selector(parts)))
        )
        if self.reject_tail:
            jets = [jet for jet in jets if jet.pt() < self.reject_tail * self.pThat]

        for jet in jets:
            self.analyze_cones(parts, jet)

    def analyze_cones(self, parts, jet):
        # coneR = np.sqrt(jet.area() / np.pi) if self.matched_cone else jetR
        coneR = self.jetR
        angle = np.pi / 2
        rot_jet_1 = fj.PseudoJet()
        rot_jet_1.reset_PtYPhiM(jet.pt(), jet.rap(), jet.phi() + angle, jet.m())
        rot_jet_2 = fj.PseudoJet()
        rot_jet_2.reset_PtYPhiM(jet.pt(), jet.rap(), jet.phi() - angle, jet.m())
        for part in parts:
            if np.sqrt((rot_jet_1.eta() - part.eta()) ** 2 + (rot_jet_1.delta_phi_to(part)) ** 2) <= coneR:
                self.analyze_cone_part(part, "P", jet.pt())
            # particles can only be in one cone or the other, so we can use elif
            elif np.sqrt((rot_jet_2.eta() - part.eta()) ** 2 + (rot_jet_2.delta_phi_to(part)) ** 2) <= coneR:
                self.analyze_cone_part(part, "M", jet.pt())

    def analyze_cone_part(self, part, sgn, jet_pT):
        p8part = pythiafjtools.getPythia8Particle(part)
        src = get_src(self.pythia, p8part)
        self.hists[src][sgn].Fill(part.pt(), jet_pT)


    def finalize(self):
        self.pythia.stat()
        self.hists['nev'].SetBinError(1, 0)

        # scale_f = pythia.info.sigmaGen() / self.hists['nev'].GetBinContent(1)
        # Divide by the sum of weights to normalize histograms to a per event basis
        # then multiply by the cross-section to normalize by cross-section
        scale_f = self.pythia.info.sigmaGen() / self.pythia.info.weightSum()

        logger.info(f"sigmaGen (Pythia) = {self.pythia.info.sigmaGen()} +- {self.pythia.info.sigmaErr()}")
        logger.info(f"Scaling factor (sigma/nev) = {scale_f} (with weightSum)")
        logger.info(f"nTried: {self.pythia.info.nTried()}")
        logger.info(f"nSelected: {self.pythia.info.nSelected()}")
        logger.info(f"nAccepted: {self.pythia.info.nAccepted()}")
        logger.info(f"Failed: {self.nfailed}")
        logger.info(f"Sum of weights (histogram): {self.hists['nev'].GetBinContent(1)}")
        logger.info(f"Sum of weights (Pythia): {self.pythia.info.weightSum()}")
        logger.info(f"Number of successful events (histogram): {self.hists['nev'].GetEntries()}")
        # with open(f"{self.output_dir}/scales.txt", 'w') as f:
        #     f.write(f"{scale_f}")
        logger.info(
            f"N total final events: {self.hists['nev'].GetBinContent(1)} with {self.pythia.info.nAccepted() - self.hists['nev'].GetBinContent(1)} events rejected."
        )

        if self.scale_by_xsec:
            self._scale_hists(scale_f)
        self._save_hists()

    def _scale_hists(self, scale_f):
        for src in [2, 3, 4, 6, 0]:
            for sgn in ["P", "M"]:
                self.hists[src][sgn].Scale(scale_f)
        # for ptype in ["T", "Q", "P", "M", "PM"]:
        #     self.hists[ptype].Scale(scale_f)
        # self.hists["jet_pT"].Scale(scale_f)
        logger.info("Histograms scaled.")

    def _save_hists(self):
        with ROOT.TFile(self.output_path, "RECREATE") as f:
            for src in [2, 3, 4, 6, 0]:
                for sgn in ["P", "M"]:
                    f.WriteTObject(self.hists[src][sgn])
        logger.info("Histograms saved.")


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
