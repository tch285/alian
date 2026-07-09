#!/usr/bin/env python

# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

import argparse
# import heppyy
# fj = heppyy.load_cppyy('fastjet')
from alian.analysis.base.analysis import AnalysisBase

class PhotonPurity(AnalysisBase):
    _defaults = {
        'iso_pt_edges': [0.0, 1.5, 4.0, 25.0],
        'm02_edges': [0.0, 0.1, 0.3, 0.4, 3.0],
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)

    def analyze_event(self):
        # Analyzes this event that has passed the selection criteria
        # self.event contains the selected event
        # self.tracks contains selected tracks (i.e. after selection)
        # self.clusters contains selected clusters (i.e. after selection)

        for c in self.clusters:
            _, iso_pt, _, _ = self.selector.isolation.selects(c, self.tracks, verbose = True)
            m02 = c.m02()
            pt = c.pt()
            E = c.e()
            self.hists['photon_E'].Fill(E)
            self.hists['photon_pt'].Fill(pt)
            self.hists['purity_fine'].Fill(m02, iso_pt, pt)
            self.hists['purity_coarse'].Fill(m02, iso_pt, pt)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser.add_argument('-i', '--input-file', type=str, help="Input file or file containing a list of input files.", required = True)
    parser.add_argument("-o", "--output-file", type=str, help="File to write the analysis.", default="analysis.root")
    parser.add_argument('-c', '--config-file', type=str, help="Path to a YAML configuration file describing the cuts to be applied to the data.", default=None)
    parser.add_argument("-n", "--nev", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
    parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
    parser.add_argument('--lhc-run', type=int, help='LHC Run', default = 3)

    args = parser.parse_args()

    ana = PhotonPurity(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()