# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python3

import argparse
import itertools
import heppyy
fj = heppyy.load_cppyy('fastjet')
from alian.analysis.base import AnalysisBase
from alian.analysis.base.utils import delta_R
alian = heppyy.load_cppyy("alian")
import numpy as np

class EECPurityWeighted(AnalysisBase):
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

    def calc_purity(self, pt):
    # def fit(x, f, x0, b, c, d):
        # return sp.erf(a*x)+B*x*np.exp(-c*x**2)+d
        return (self.params[0]-(self.params[3]+self.params[4]*(pt-self.params[1]))*np.exp(-self.params[2]*(pt-self.params[1]))) / 100

    def analyze_event(self):
        # Analyzes this event that has passed the selection criteria
        # self.event contains the selected event
        # self.tracks contains selected tracks (i.e. after selection)
        # self.clusters contains selected clusters (i.e. after selection)

        jets = self.jet_finder.find_jets(self.tracks)

        for c in self.clusters:
            # if cluster is not above threshold or is not isolated, skip
            if c.pt() < self.cluster_pt_min or not self.selector.isolation.selects(c, self.tracks):
                continue

            # check if cluster is wide or narrow, if neither then skip
            purity = self.calc_purity(c.pt())
            if c.m02() > self.m02_narrow_min and c.m02() < self.m02_narrow_max:
                clus_type = "narrow"
                pweight = 1 / purity
            elif c.m02() > self.m02_wide_min and c.m02() < self.m02_wide_max:
                clus_type = "wide"
                pweight = (1 - purity) / purity
            else:
                continue

            # find back-to-back jets
            jets_b2b = [j for j in jets if np.abs(j.delta_phi_to(c)) > self.photon_jet_angle_min]
            if jets_b2b:
                self.hists[f"photon_pT_{clus_type}"].Fill(c.pt())
                for j in jets_b2b:
                    self.do_eec(j, clus_type, pweight)

    def do_eec(self, jet, clus_type, pweight):
        self.hists[f"jet_pT_{clus_type}"].Fill(jet.pt())
        tracks = self.eec_trk_selector(jet.constituents())
        for p1, p2 in itertools.permutations(tracks, 2):
            ew = p1.pt() * p2.pt() / jet.pt() / jet.pt()
            angle = delta_R(p1, p2)
            q1 = p1.user_info[alian.TrackInfo]().q()
            q2 = p2.user_info[alian.TrackInfo]().q()

            self.hists[f"eec_T_{clus_type}"].Fill(jet.pt(), angle, ew * pweight)
            self.hists[f"eec_Q_{clus_type}"].Fill(jet.pt(), angle, ew * q1 * q2 * pweight)
            if q1 > 0 and q2 > 0:
                self.hists[f"eec_P_{clus_type}"].Fill(jet.pt(), angle, ew * pweight)
            elif q1 < 0 and q2 < 0:
                self.hists[f"eec_M_{clus_type}"].Fill(jet.pt(), angle, ew * pweight)
            else:
                self.hists[f"eec_PM_{clus_type}"].Fill(jet.pt(), angle, ew * pweight)


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

    ana = EECPurityWeighted(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()