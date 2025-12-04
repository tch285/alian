import heppyy
from alian.analysis.base import AnalysisSelector, Event, JetFinder, Track, Cluster, HistogramCollection
from alian.io import data_io
from alian.analysis.base.logging import setup_logger
import numpy as np
fj = heppyy.load_cppyy('fastjet')

class BaseAnalysis(heppyy.GenericObject):
    _defaults = {}

    def __init__(self, **kwargs):
        super(BaseAnalysis, self).__init__(**kwargs)
        for k, val in self.__class__._defaults.items():
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, val)

class AnalysisBase:
    def __init__(self, input_file, output_file, cfg, tree_struct, nev = -1, lhc_run = 3):
        self.input_file = input_file
        self.output_file = output_file
        self.cfg = cfg
        self.nev = nev
        self.tree_struct = tree_struct
        self.lhc_run = lhc_run

    def run(self):
        self.init()
        for e in self.data_source.next_event():
            event = Event(e)
            if not self.selector.selects(event, verbose = False):
                continue
            self.analyze(event)
        self._finalize()

    def init(self):
        # initialize data input
        self.data_source = data_io.DataInput(self.input_file,
                                             lhc_run = self.lhc_run,
                                             yaml_file = self.tree_struct,
                                             n_events = self.nev)
        # initialize selectors for events, tracks, clusters
        self.selector = AnalysisSelector.from_file(self.cfg)
        self.selector.dump()
        # initialize jet finder
        self.jet_finder = JetFinder.from_file(self.cfg)
        self.logger = setup_logger(__name__)
        # initialize analysis parameters
        self.init_analysis(self.cfg['analysis'])
        self.init_output(self.cfg['bins'], self.cfg['histograms'])

    def init_analysis(self, analysis_cfg):
        """Initialize analysis configuration parameters."""
        raise NotImplementedError("Must implement an init_analysis method!")

    def analyze(self, event):
        tracks = [t for t in event.tracks if self.selector.track.selects(t, verbose = False)]
        clusters = [c for c in event.clusters if self.selector.cluster.selects(c, tracks)]
        self.jet_finder.find_jets(tracks)
        self.analyze_event(self, tracks, clusters, self.jet_finder.jets)

    def analyze_event(self, tracks, clusters, jets):
        raise NotImplementedError("Must implement an analyze_event method!")

    def _finalize(self):
        self.finalize()
        self.save_histograms()

    def finalize(self):
        # finalize analysis before saving histograms
        pass