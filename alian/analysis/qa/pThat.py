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

import numpy as np
from alian.analysis.base import AnalysisMCBase, add_default_args

import heppyy

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")


class QaPtHat(AnalysisMCBase):
    _defaults = {
        "pThat_scale_det": 1000,
        "pThat_scale_gen": 1000
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)

    def analyze_event(self):
        self.hists['pThat'].Fill(self.pThat)
        self.hists['pThat_weighted'].Fill(self.pThat, self.weight)
        pThat_calc = 10 / np.sqrt(np.sqrt(self.weight))
        self.hists['pThat_corr'].Fill(self.pThat, pThat_calc)
        self.hists['pThat_diff'].Fill(self.pThat - pThat_calc)
        self.tracks = fj.sorted_by_pt(self.tracks)
        self.particles = fj.sorted_by_pt(self.particles)

        if self.tracks:
            self.hists['ratio_track'].Fill(self.tracks[0].pt() / self.pThat)
            for track in self.tracks:
                self.hists['track_pT'].Fill(track.pt(), self.weight)
        if self.particles:
            self.hists['ratio_particle'].Fill(self.particles[0].pt() / self.pThat)
            for particle in self.particles:
                self.hists['particle_pT'].Fill(particle.pt(), self.weight)
        if self.tracks and self.particles:
            self.hists['ratio_track_corr'].Fill(self.tracks[0].pt() / self.pThat, self.particles[0].pt() / self.pThat)

        if self.jets_det:
            self.hists['ratio_jet_det'].Fill(self.jets_det[0].pt() / self.pThat)
            for jet_det in self.jets_det:
                self.hists['jet_pT_det'].Fill(jet_det.pt(), self.weight)
        if self.jets_gen:
            self.hists['ratio_jet_gen'].Fill(self.jets_gen[0].pt() / self.pThat)
            for jet_gen in self.jets_gen:
                self.hists['jet_pT_gen'].Fill(jet_gen.pt(), self.weight)
        if self.jets_det and self.jets_gen:
            self.hists['ratio_jet_corr'].Fill(self.jets_det[0].pt() / self.pThat, self.jets_gen[0].pt() / self.pThat)

        pass_cut = True
        if self.jets_det and (self.jets_det[0].pt() / self.pThat) > self.pThat_scale_det:
            pass_cut = False
        if self.jets_gen and (self.jets_gen[0].pt() / self.pThat) > self.pThat_scale_gen:
            pass_cut = False
        if pass_cut:
            for jet_det in self.jets_det:
                self.hists['jet_pT_det_after_jet'].Fill(jet_det.pt(), self.weight)
            for jet_gen in self.jets_gen:
                self.hists['jet_pT_gen_after_jet'].Fill(jet_gen.pt(), self.weight)
        # if self.tracks and (self.tracks[0].pt() / self.pThat) < self.pThat_scale_det and \
        #    self.particles and (self.particles[0].pt() / self.pThat) < self.pThat_scale_gen:
            for track in self.tracks:
                self.hists['track_pT_det_after_jet'].Fill(track.pt(), self.weight)
            for particle in self.particles:
                self.hists['track_pT_gen_after_jet'].Fill(particle.pt(), self.weight)
        else:
            if self.jets_det or self.jets_gen:
                self.logger.warning(f"Rejecting event: det jets {len(self.jets_det)}, gen jets {len(self.jets_gen)}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = QaPtHat(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
