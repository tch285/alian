#!/usr/bin/env python

# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)


import argparse

from alian.analysis.base import AnalysisMCBase, add_default_args

import heppyy

alian = heppyy.load_cppyy("alian")

class JetQA(AnalysisMCBase):
    def analyze_event(self):
        for jet_det in self.jets_det:
            self.hists["jet_pT_det"].Fill(jet_det.pt(), self.weight)
            self.hists["jet_eta_det"].Fill(jet_det.eta(), self.weight)
            self.hists["jet_phi_det"].Fill(jet_det.phi(), self.weight)
            self.hists["jet_eta_phi_det"].Fill(jet_det.eta(), jet_det.phi(), self.weight)
            self.hists["jet_eta_phi_pT_det"].Fill(jet_det.eta(), jet_det.phi(), jet_det.pt(), self.weight)
        for jet_gen in self.jets_gen:
            self.hists["jet_pT_gen"].Fill(jet_gen.pt(), self.weight)
            self.hists["jet_eta_gen"].Fill(jet_gen.eta(), self.weight)
            self.hists["jet_phi_gen"].Fill(jet_gen.phi(), self.weight)
            self.hists["jet_eta_phi_gen"].Fill(jet_gen.eta(), jet_gen.phi(), self.weight)
            self.hists["jet_eta_phi_pT_gen"].Fill(jet_gen.eta(), jet_gen.phi(), jet_gen.pt(), self.weight)

        pairs_py = self._get_jet_matches(self.jets_det, self.jets_gen, 0.4*0.6)
        pairs_cpp = self._get_jet_matches_cpp(self.jets_det, self.jets_gen, 0.4*0.6)
        if pairs_py != pairs_cpp:
            raise ValueError(f"Pairs don't match:\n\tPython: {pairs_py}\n\tC++: {pairs_cpp}")
        for idx_det, idx_gen in pairs_cpp:
            jet_det_matched = self.jets_det[idx_det]
            self.hists["jet_pT_det_matched"].Fill(jet_det_matched.pt(), self.weight)
            self.hists["jet_eta_det_matched"].Fill(jet_det_matched.eta(), self.weight)
            self.hists["jet_phi_det_matched"].Fill(jet_det_matched.phi(), self.weight)
            self.hists["jet_eta_phi_det_matched"].Fill(jet_det_matched.eta(), jet_det_matched.phi(), self.weight)
            self.hists["jet_eta_phi_pT_det_matched"].Fill(jet_det_matched.eta(), jet_det_matched.phi(), jet_det_matched.pt(), self.weight)

            jet_gen_matched = self.jets_gen[idx_gen]
            self.hists["jet_pT_gen_matched"].Fill(jet_gen_matched.pt(), self.weight)
            self.hists["jet_eta_gen_matched"].Fill(jet_gen_matched.eta(), self.weight)
            self.hists["jet_phi_gen_matched"].Fill(jet_gen_matched.phi(), self.weight)
            self.hists["jet_eta_phi_gen_matched"].Fill(jet_gen_matched.eta(), jet_gen_matched.phi(), self.weight)
            self.hists["jet_eta_phi_pT_gen_matched"].Fill(jet_gen_matched.eta(), jet_gen_matched.phi(), jet_gen_matched.pt(), self.weight)

            self.hists['jet_pT_resp'].Fill(jet_det_matched.pt(), jet_gen_matched.pt(), self.weight)

            res = (jet_det_matched.pt() - jet_gen_matched.pt()) / jet_gen_matched.pt() * 100
            self.hists['jet_pT_res'].Fill(jet_gen_matched.pt(), res, self.weight)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)
    args = parser.parse_args()

    ana = JetQA(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()