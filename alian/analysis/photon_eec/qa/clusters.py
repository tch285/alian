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
from alian.analysis.base.selector import AnalysisSelector
from alian.analysis.base.logs import set_up_logger

def linbins(xmin, xmax, nbins):
    return np.linspace(xmin, xmax, nbins+1)
def logbins(xmin, xmax, nbins):
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)

ROOT.TH1.SetDefaultSumw2() # applies to TH2 and TH3 as well

class ClusterQA:
    def __init__(self, inf, outf, cfg):
        self.outf = outf
        self.selector = AnalysisSelector.load(cfg)
        self.df = ROOT.RDataFrame("eventTree", inf)
        self.histos = {}

    def note_time(self, msg):
        logger.info(f"{msg}: -------- {timenow() - self.start_time:.3f} sec. --------", stacklevel = 2)
    def note_start(self, msg):
        self.start_time = timenow()
        logger.info(f"{msg}: {dt.now().replace(microsecond=0)}", stacklevel = 2)

    def analyze(self):
        self.note_start("Starting analysis")
        self.apply_cuts()
        self.create_histos()
        self.fill_and_save()
        self.note_time("Finished analysis")

    def apply_cuts(self):
        self.df = self.selector.event.apply_to(self.df)
        self.df = self.df.Define("cluster_pt", "cluster_energy / cosh(cluster_eta)")
        self.note_time("Cuts applied")

    def create_histos(self):
        self.histos['E'] = self.df.Histo1D(
            ("E", 'Cluster energy distribution;#it{E} (GeV);Counts', 200, 0, 100),
            "cluster_energy"
        )
        self.histos['pT'] = self.df.Histo1D(
            ("pT", 'Cluster #it{p}_{T} distribution;#it{p}_{T} (GeV/#it{c});Counts', 200, 0, 100),
            "cluster_pt"
        )
        self.histos['m02'] = self.df.Histo1D(
            ("m02", 'Cluster shape;#it{M}_{02};Counts', 400, 0, 2),
            "cluster_m02"
        )
        self.histos['time'] = self.df.Histo1D(
            ("time", 'Cluster time;time (ns);Counts', 120, -30, 30),
            "cluster_time"
        )
        self.histos['ncells'] = self.df.Histo1D(
            ("ncells", 'Cluster cell multiplicity;n_{cells};Counts', 25, -0.5, 24.5),
            "cluster_ncells"
        )
        self.histos['exoticity'] = self.df.Histo1D(
            ("exoticity", 'Cluster exoticity;exoticity;Counts', 2, -0.5, 1.5),
            "cluster_exoticity"
        )
        self.histos['nlm'] = self.df.Histo1D(
            ("nlm", 'Number of local maxima;NLM;Counts', 10, -0.5, 9.5),
            "cluster_nlm"
        )
        self.histos['m02_pt'] = self.df.Histo2D(
            ("m02_pt", 'Cluster #it{p}_{T} vs. #it{M}_{02};#it{M}_{02};#it{p}_{T} (GeV/#it{c})', 400, 0, 2, 100, 0, 100),
            "cluster_m02", "cluster_pt"
        )
        self.histos['time_pt'] = self.df.Histo2D(
            ("time_pt", 'Cluster #it{p}_{T} vs. time;time (ns);#it{p}_{T} (GeV/#it{c})', 120, -30, 30, 100, 0, 100),
            "cluster_time", "cluster_pt"
        )
        self.histos['pt_exoticity'] = self.df.Histo2D(
            ("pt_exoticity", 'Cluster exoticity vs. #it{p}_{T};#it{p}_{T};exoticity', 100, 0, 100, 2, -0.5, 1.5),
            "cluster_pt", "cluster_exoticity"
        )
        self.histos['eta_phi'] = self.df.Histo2D(
            ("eta_phi", 'Cluster #eta-#phi distribution;#eta;#phi', 200, -0.9, 0.9, 200, 0, 2*np.pi),
            "cluster_eta", "cluster_phi"
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
        logger.info(f"Output: {self.outf}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="Input file or file containing a list of input files.", required = True)
    parser.add_argument("-o", "--output", type=str, help="File to write the analysis.", required = True)
    parser.add_argument("-c", "--config", type=str, help="YAML file describing the cuts to be applied to the data.", default=None)
    args = parser.parse_args()
    logger = set_up_logger(__name__)

    ROOT.EnableImplicitMT()
    analysis = ClusterQA(args.input, args.output, args.config)
    analysis.analyze()