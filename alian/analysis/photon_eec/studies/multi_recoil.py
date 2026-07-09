# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python

import argparse
import numpy as np
import heppyy
from alian.analysis.base.analysis import AnalysisBase
alian = heppyy.load_cppyy("alian")

class MultiRecoilAnalysis(AnalysisBase):
    _defaults = {
        'photon_jet_angle_min': 0.875 * np.pi,
        'cluster_pt_min': 10,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)

    def analyze_event(self):
        jets = self.jet_finder.find_jets(self.tracks)

        clus_sel = [c for c in self.clusters if c.pt() > self.cluster_pt_min and self.selector.isolation.selects(c, self.tracks)]
        for j in jets:
            clus_recoil = [c for c in clus_sel if np.abs(j.delta_phi_to(c)) > self.photon_jet_angle_min]
            nrecoil = len(clus_recoil)
            self.hists['nrecoil'].Fill(nrecoil)
            self.hists['nrecoil2d'].Fill(nrecoil, j.pt())



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser.add_argument('-i', '--input-file', type=str, help="Input file or file containing a list of input files.", required = True)
    parser.add_argument("-o", "--output-file", type=str, help="File to write the analysis.", default="analysis.root")
    parser.add_argument('-c', '--config-file', type=str, help="Path to a YAML configuration file describing the cuts to be applied to the data.", default=None)
    parser.add_argument("-n", "--nev", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
    parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
    parser.add_argument('--lhc-run', type=int, help='LHC Run', default = 3)

    args = parser.parse_args()

    ana = MultiRecoilAnalysis(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()