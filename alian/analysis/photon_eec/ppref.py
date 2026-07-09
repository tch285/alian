# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python

import itertools
import heppyy
from alian.analysis.base import AnalysisBase
from alian.analysis.base.analysis import get_default_args
from alian.analysis.base.utils import delta_R

fj = heppyy.load_cppyy('fastjet')
alian = heppyy.load_cppyy("alian")

class BasicEEC(AnalysisBase):
    _defaults = {
        'pt_min_eec': 1.0,
    }
    def init_analysis(self, analysis_cfg: dict):
        config = self._defaults | analysis_cfg
        for setting, value in config.items():
            setattr(self, setting, value)
        self.eec_trk_selector = fj.SelectorPtMin(self.pt_min_eec)

    def analyze_event(self):
        jets = self.jet_finder.find_jets(self.tracks)
        for jet in jets:
            self.do_eec(jet)

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


if __name__ == '__main__':
    parser = get_default_args()
    args = parser.parse_args()

    ana = BasicEEC(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()