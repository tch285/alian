#!/usr/bin/env python3

"""Example script on using the alian framework for Run 3 analysis.

This script is used on the ROOT file (provided as an argument) to
perform the analysis of the events containing isolated photons and
jets. The analysis is configured with a config YAML file, creates
and fills histograms for the cleaned data, and stores them in the
output file.

To be used with alian/config/example/example_analysis.yaml
"""

import argparse
import itertools

import numpy as np
from alian.analysis.base import AnalysisMCBase, add_default_args, delta_R

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class AnalysisExample(AnalysisMCBase):
    _defaults = {
        'pt_min_eec': 1.0,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)

    def analyze_event(self):
        # Analyzes this event that has passed the selection criteria
        # self.event contains the selected event
        # self.tracks contains selected tracks (i.e. after selection cuts)
        # self.clusters contains selected clusters (i.e. after selection cuts)
        # self.jets contains selected jets

        [self.hists['track_pT_det'].Fill(t.pt(), self.weight) for t in self.tracks]
        [self.hists['track_pT_gen'].Fill(t.pt(), self.weight) for t in self.particles]

        for p in self.particles:
            self.hists['eff_track_pT_gen'].Fill(p.pt(), self.weight)
            if p.user_info[alian.TrackInfo]().has_match():
                self.hists['eff_track_pT_det'].Fill(p.pt(), self.weight)

                match = p.user_info[alian.TrackInfo]().match()
                match_pt = match.pt()
                res = (match_pt - p.pt()) / p.pt() * 100
                self.hists['pT_res'].Fill(p.pt(), res, self.weight)

                self.hists['ch_track_pT_gen'].Fill(p.pt(), self.weight)
                if match.user_info[alian.TrackInfo]().q() == p.user_info[alian.TrackInfo]().q():
                    self.hists['ch_track_pT_det'].Fill(p.pt(), self.weight)
        for t in self.tracks:
            self.hists['pur_track_pT_det'].Fill(t.pt(), self.weight)
            if t.user_info[alian.TrackInfo]().has_match():
                self.hists['pur_track_pT_gen'].Fill(t.pt(), self.weight)

        [self.hists['jet_pT_det'].Fill(j.pt(), self.weight) for j in self.jets_det]
        [self.hists['jet_pT_gen'].Fill(j.pt(), self.weight) for j in self.jets_gen]
        for j in self.jets_det:
            self.hists['pThat_ratio_det'].Fill(j.pt() / self.event.pThat)
            self.do_eec(j, "det")
        for j in self.jets_gen:
            self.hists['pThat_ratio_gen'].Fill(j.pt() / self.event.pThat)
            self.do_eec(j, "gen")

    def do_eec(self, jet, suffix):
        tracks = self.eec_trk_selector(jet.constituents())

        for p1, p2 in itertools.permutations(tracks, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()

            self.hists[f"eec_T_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            self.hists[f"eec_Q_{suffix}"].Fill(jet.pt(), angle, ew * self.weight * q1 * q2)
            if q1 > 0 and q2 > 0:
                self.hists[f"eec_P_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            if q1 < 0 and q2 < 0:
                self.hists[f"eec_M_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)
            else:
                self.hists[f"eec_PM_{suffix}"].Fill(jet.pt(), angle, ew * self.weight)

    def finalize(self):
        self.hists['track_pT_det'].Scale(1, "width")
        self.hists['track_pT_gen'].Scale(1, "width")

        # h_eff = self.hists['eff_track_pT_det'].Clone("h_eff")
        # h_eff.SetTitle("Efficiency;pt (GeV);Efficiency")
        # h_eff.Divide(self.hists['eff_track_pT_gen'])
        # self.hists['efficiency'] = h_eff

        # h_pur = self.hists['pur_track_pT_gen'].Clone("h_pur")
        # h_pur.SetTitle("Purity;pt (GeV);Purity")
        # h_pur.Divide(self.hists['pur_track_pT_det'])
        # self.hists['purity'] = h_pur

        # h_ch = self.hists['ch_track_pT_det'].Clone("h_ch")
        # h_ch.SetTitle("Charge reco efficiency;pt (GeV);Efficiency")
        # h_ch.Divide(self.hists['ch_track_pT_gen'])
        # self.hists['ch_eff'] = h_ch


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = AnalysisExample(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
