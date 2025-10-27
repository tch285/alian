#!/usr/bin/env python3

"""Example script on using RDataFrame for simple QA.

An example QA script that uses RDataFrame for basic QA tasks.
This use of RDF is a little more involved than qa.py: the code
is more complicated, but it is also easier to do more with qa_ext.py
out of the box. RDF implementations can only handle basic 
selections on cluster and track quantities; for more complicated
tasks like jet finding or cluster-track matching, use example.py.
"""

import argparse
from datetime import datetime as dt
from time import perf_counter as timenow

import numpy as np
import ROOT
from alian.analysis.base import AnalysisSelector
from alian.analysis.base.logging import setup_logger


class AnalysisQA:
    def __init__(self, inf, outf, cfg):
        self.outf = outf
        self.selector = AnalysisSelector.from_file(cfg)
        self.df = ROOT.RDataFrame("eventTree", inf)
        self.histos = {}

    def note_time(self, msg):
        logger.info(f"{msg}: -------- {timenow() - self.start_time:.3f} sec. --------", stacklevel = 2)
    def note_start(self, msg):
        self.start_time = timenow()
        logger.info(f"{msg}: {dt.now().replace(microsecond=0)}", stacklevel = 2)

    def analyze(self):
        self.note_start("Starting analysis")
        # print(self.selector.event)
        self.apply_cuts()
        self.create_histos()
        self.fill_and_save()
        self.note_time("Finished analysis")

    def apply_cuts(self):
        """
        Apply cuts to events, tracks, and clusters. The event selection
        is just a Filter applied on the RDF. Cluster and track selections
        are a little more complicated since they're vector columns, so
        the implementation is hidden behind apply_to(). This method
        creates a new set of columns for tracks and clusters with the
        relevant cuts applied to them.
        For tracks, the created columns are: pt, phi and eta;
        and for clusters they are: energy, m02, time, and ncells.
        One should use these column names to fill histograms from clusters
        or tracks with the cuts applied, or use the default track_data_pt,
        track_data_phi; etc. with these
        """
        self.df = self.selector.event.apply_to(self.df)
        self.df_cluster = self.selector.cluster.apply_to(self.df)
        self.df_track = self.selector.track.apply_to(self.df)
        self.note_time("Cuts applied")

    def create_histos(self):
        """
        Define your histograms here. For Histo*D, the first argument is
        a tuple corresponding to what you would use for a normal TH*.
        The second is the column that your histogram will be Filled with.
        Note that the histograms themselves are not filled until
        fill_and_save() is called.
        """
        self.histos['track_pt'] = self.df_track.Histo1D(
            ("track_pt", 'Track #it{p}_T;track #it{p}_T (GeV);Counts', 50, 0, 30),
            "pt"
        )
        self.histos['track_eta'] = self.df_track.Histo1D(
            ("track_eta", 'Track #eta;track #eta;Counts', 50, -0.9, 0.9),
            "eta"
        )
        self.histos['track_phi'] = self.df_track.Histo1D(
            ("track_phi", 'Track #phi;track #phi;Counts', 50, 0, 2*np.pi),
            "phi"
        )
        self.histos['track_eta_phi'] = self.df_track.Histo2D(
            ("track_eta_phi", '#eta-#phi distribution;#eta;#phi', 50, -0.9, 0.9, 50, 0, 2*np.pi),
            "eta", "phi"
        )
        self.histos['cluster_E'] = self.df_cluster.Histo1D(
            ("cluster_E", 'Cluster energy distribution;#it{E} (GeV);Counts', 50, 0, 25),
            "energy"
        )
        self.note_time("Histograms defined")

    def fill_and_save(self):
        """
        Fill and save all defined histograms. Note that, for histograms
        defined in RDFs, filling is triggered when Write() is called.
        """
        with ROOT.TFile(self.outf, 'recreate'):
            for h in self.histos.values():
                h.Write()
        self.note_time("Histograms saved")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="Input file or file containing a list of input files.", required = True)
    parser.add_argument("-o", "--output", type=str, help="File to write the analysis.", required = True)
    parser.add_argument("-c", "--config", type=str, help="YAML file describing the cuts to be applied to the data.", default=None)

    args = parser.parse_args()
    logger = setup_logger(__name__)

    analysis = AnalysisQA(args.input, args.output, args.config)
    analysis.analyze()