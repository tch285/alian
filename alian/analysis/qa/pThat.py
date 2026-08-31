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


class qaPtHat(AnalysisMCBase):
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
        if self.jets_det:
            self.hists['pThat_ratio_det'].Fill(self.jets_det[0].pt() / self.event.pThat)
            for jet_det in self.jets_det:
                self.hists['jet_pT_det'].Fill(jet_det.pt(), self.weight)
        if self.jets_gen:
            self.hists['pThat_ratio_gen'].Fill(self.jets_gen[0].pt() / self.event.pThat)
            for jet_gen in self.jets_gen:
                self.hists['jet_pT_gen'].Fill(jet_gen.pt(), self.weight)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser = add_default_args(parser)

    args = parser.parse_args()

    ana = qaPtHat(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()
