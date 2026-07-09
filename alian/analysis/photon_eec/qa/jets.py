# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python

from alian.analysis.base import AnalysisBase
from alian.analysis.base.analysis import get_default_args

class JetQA(AnalysisBase):
    def analyze_event(self):
        jets = self.jet_finder.find_jets(self.tracks)
        for jet in jets:
            self.hists["pT"].Fill(jet.pt())
            self.hists["phi"].Fill(jet.phi())
            self.hists["eta"].Fill(jet.eta())
            self.hists["eta_phi"].Fill(jet.eta(), jet.phi())
            self.hists["eta_phi_pT"].Fill(jet.eta(), jet.phi(), jet.pt())

if __name__ == '__main__':
    parser = get_default_args()
    args = parser.parse_args()

    ana = JetQA(args.input_file, args.output_file, args.config_file, args.tree_struct, args.nev, args.lhc_run)
    ana.run()