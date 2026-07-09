# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python

import argparse
import itertools
import numpy as np
import heppyy
from alian.analysis.base import AnalysisBase
from alian.analysis.base.utils import delta_R
fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")

class AnalysisExample(AnalysisBase):
    _defaults = {
        'photon_jet_angle_min': 0.875 * np.pi,
        'cluster_pt_min': 10,
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
        # self.tracks contains selected tracks (i.e. after selection)
        # self.clusters contains selected clusters (i.e. after selection)

        # [self.hists['track_pt'].Fill(t.pt()) for t in self.tracks]
        # jets = self.jet_finder.find_jets(self.tracks)
        # [self.hists['jet_pT'].Fill(j.pt()) for j in jets]

        for c in self.clusters:
            if c.pt() < self.cluster_pt_min:
                continue
            is_iso, iso_pt, _, _ = self.selector.isolation.selects(c, self.tracks, verbose = True)
            # self.hists['tot_iso_pt'].Fill(iso_pt)
            if is_iso:
                self.hists['iso_photon_E'].Fill(c.e())
                # find back-to-back jets
                for j in self.jets:
                    dphi = np.abs(j.delta_phi_to(c))
                    if dphi > self.photon_jet_angle_min:
                        self.hists["jet_photon_dphi"].Fill(dphi)
                        # self.hists["xjy"].Fill(j.pt() / c.e(), c.e())
                        self.do_eec(j)

    def do_eec(self, jet):
        self.hists["jet_pT"].Fill(jet.pt())
        tracks = self.eec_trk_selector(jet.constituents())
        for p1, p2 in itertools.permutations(tracks, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()

            self.hists["eec_T"].Fill(jet.pt(), angle, ew)
            self.hists["eec_Q"].Fill(jet.pt(), angle, ew * q1 * q2)
            if q1 > 0 and q2 > 0:
                self.hists["eec_P"].Fill(jet.pt(), angle, ew)
            elif q1 < 0 and q2 < 0:
                self.hists["eec_M"].Fill(jet.pt(), angle, ew)
            else:
                self.hists["eec_PM"].Fill(jet.pt(), angle, ew)


    # def finalize(self):
    #     self.hists['track_pt'].Scale(1, "width")

    def dump(self):
        cfg = "\n".join([f"\t{param}: {repr(getattr(self, param))}" for param in self._defaults])
        self.logger.info(f"{type(self).__name__} configuration:\n{cfg}", stacklevel = 2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser.add_argument('-i', '--input-file', type=str, help="Input file or file containing a list of input files.", required = True)
    parser.add_argument("-o", "--output-file", type=str, help="File to write the analysis.", default="analysis.root")
    parser.add_argument('-c', '--config-file', type=str, help="Path to a YAML configuration file describing the cuts to be applied to the data.", default=None)
    parser.add_argument("-n", "--nev", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
    parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
    parser.add_argument('--lhc-run', type=int, help='LHC Run', default = 3)

    args = parser.parse_args()

    ana = AnalysisExample(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()