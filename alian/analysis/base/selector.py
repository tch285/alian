import yaml
import numpy as np
import alian.analysis.base.selection as sel

class Selector:
    _defaults = {}
    _nosel = {}
    _selfid = ''
    def __init__(self, use_defaults = False, **selections):
        self.set_selection(use_defaults = use_defaults, **selections)

    def set_selection(self, use_defaults, **selections):
        """Set selections. Unspecified selections are set to extremal (no cut) values."""
        self.reset_selection(use_defaults)
        self.update_selection(**selections)

    def reset_selection(self, use_defaults = False):
        """Reset selection to default or no-selection values."""
        if use_defaults:
            for selection, value in self._defaults.items():
                setattr(self, selection, value)
        else:
            for selection, value in self._nosel.items():
                setattr(self, selection, value)

    def update_selection(self, **selections):
        valid_sels = self._validate_selection(**selections)
        self._update_selection(**valid_sels)

    def _update_selection(self, **selections):
        raise NotImplementedError

    def apply_to(self, df):
        """Apply cuts to the event RDataFrame."""
        raise NotImplementedError

    def _validate_selection(self, **selections):
        """Validate the passed selection, and raise warnings if invalid selections were given."""
        bad_sels = [sel for sel in selections.keys() if sel not in self._defaults]
        if bad_sels:
            print(f"{self._selfid}: invalid cuts given, will be ignored: {', '.join(bad_sels)}")
        return {sel: val for sel, val in selections.items() if sel in self._defaults}

    def selects(self, obj):
        raise NotImplementedError

    def dump(self):
        print(f"{self._selfid}:\n" + "\n".join([f"\t{key}: {getattr(self, key)}" for key in self._defaults]))

class EventSelector(Selector):
    _defaults = {
        "event_selection": sel.EvSel.sel8,
        "triggersel": sel.TrgSel.GammaHighPtEMCAL
    }
    _nosel = {
        "event_selection": ~sel.EvSel(0),
        "triggersel": ~sel.TrgSel(0)
    }
    _selfid = "Event selector"

    def _update_selection(self, **selections):
        if "event_selection" in selections:
            self.event_selection = sel.EvSel(selections["event_selection"])
        if "triggersel" in selections:
            self.triggersel = sel.TrgSel(selections["triggersel"])

    def apply_to(self, df):
        """Apply event selections to the event RDataFrame."""
        return df.Filter(f"(event_selection & {self.event_selection}) && (triggersel & {self.triggersel})")

    def selects(self, ev, verbose = False):
        evsel_pass = bool(self.event_selection & ev.event_selection)
        trgsel_pass = bool(self.triggersel & ev.triggersel)
        if verbose:
            return evsel_pass and trgsel_pass, evsel_pass, trgsel_pass
        else:
            return evsel_pass and trgsel_pass

class ClusterSelector(Selector):
    _defaults = {
        # cluster cuts
        "E_min": 0.7, #  min cluster energy (GeV)
        "m02_min": 0.1, # min m02
        "m02_max": 0.3, # max m02
        "time_min": -30, # min relative timing (ns)
        "time_max": 35, # max relative timing (ns)
        "ncells_min": 2, # min number of cells
        # track-cluster matching cuts
        "delta_phi": 0.05, # min angle in phi to track match
        "delta_eta": 0.05, # min distance in pseudorapidity to track match
        "e_pt_max": 1.75, # max Eclus/track pT ratio to track match
        # isolation cuts
        "iso_cone_R": 0.4, # isocone radius
        "iso_pt_max": 1.5, # max isocone pT (GeV)
        # FIXME: these cuts are filled with dummy values in JE derived data,
        # so these should be optimized once they are properly implemented
        # in O2Physics.
        #   distancebadchannel currently filled with 1024
        #   nlm currently filled with zeros
        "distancebadchannel": 0, # min distance to a bad EMC channel
        "nlm": 1, # max number of cluster local maxima
    }
    _nosel = {
        # cluster cuts
        "E_min": -1, #  min cluster energy (GeV)
        "m02_min": -1, # min m02
        "m02_max": 10000, # max m02
        "time_min": -5000, # min relative timing (ns)
        "time_max": 5000, # max relative timing (ns)
        "ncells_min": -1, # min number of cells
        # FIXME: these cuts are filled with dummy values in JE derived data,
        # so these should be optimized once they are properly implemented
        # in O2Physics.
        #   distancebadchannel currently filled with 1024
        #   nlm currently filled with zeros
        "distancebadchannel": -1, # min distance to a bad EMC channel
        "nlm": 100, # max number of cluster local maxima
        # track-cluster matching cuts
        "delta_phi": -1, # min angle in phi to track match
        "delta_eta": -1, # min distance in pseudorapidity to track match
        "e_pt_max": -1, # max Eclus/track pT ratio to track match
        # isolation cuts
        "iso_cone_R": 0.4, # isocone radius
        "iso_pt_max": 5000, # max isocone pT (GeV)
    }
    _selfid = 'Cluster selector'

    def _update_selection(self, **selections):
        for name, val in selections.items():
            setattr(self, name, val)

    def apply_to(self, df):
        mask = (f"(cluster_data_energy >= {self.E_min}) "
             f"&& (cluster_data_m02 >= {self.m02_min}) "
             f"&& (cluster_data_m02 <= {self.m02_max}) "
             f"&& (cluster_data_time >= {self.time_min}) "
             f"&& (cluster_data_time <= {self.time_max}) "
             f"&& (cluster_data_ncells >= {self.ncells_min}) "
             f"&& (cluster_data_distancebadchannel >= {self.distancebadchannel}) "
             f"&& (cluster_data_nlm <= {self.nlm}) "
             )
        ret = df.Define("cluster_mask", mask)
        for colname, newname in {"cluster_data_energy": "energy",
                                 "cluster_data_m02":    "m02",
                                 "cluster_data_time":   "time",
                                 "cluster_data_ncells": "ncells",
                                 }.items():
            ret = ret.Define(newname, f"{colname}[cluster_mask]")
        return ret

    def is_geo_matched(self, cluster, track):
        # TODO: phi* matching instead of phi
        # TODO: KDTree based matching
        # TODO: use matchedTrackIndex from AO2D
        """Checks if a cluster and track match in eta and phi."""
        return np.abs(cluster.delta_eta(track)) <= self.delta_eta and \
            np.abs(cluster.delta_phi(track)) <= self.delta_phi
    def is_e_pt_matched(self, cluster, track):
        """Checks if a cluster and track match in energy."""
        return cluster.energy / track.p <= self.e_pt_max
    def is_not_in_sector_edge(self, cluster, track):
        # TODO: implement TPC sector edge cut from LF
        return True
    def is_matched(self, cluster, track):
        """Checks if a cluster and track match."""
        return self.is_geo_matched(cluster, track) and self.is_e_pt_matched(cluster, track) \
            and self.is_not_in_sector_edge(cluster, track)
    def pass_shape(self, cluster):
        """Checks if the cluster passes the shape cut."""
        return cluster.m02 > self.m02_min and cluster.m02 < self.m02_max
    def pass_time(self, cluster):
        """Checks if the cluster passes the time range cut."""
        return cluster.time > self.time_min and cluster.time < self.time_max
    def pass_ncells(self, cluster):
        """Checks if the cluster has at least the minimum number of cells."""
        return cluster.ncells >= self.ncells_min
    def pass_E(self, cluster):
        """Checks if the cluster has enough energy."""
        return cluster.energy > self.E_min
    def pass_exotic(self, cluster):
        """Checks if the cluster passes the exoticity cut."""
        return not cluster.isexotic
    def pass_dbc(self, cluster):
        """Checks if the cluster is far enough away from a bad EMC channel."""
        return cluster.distancebadchannel >= self.distancebadchannel
    def pass_nlm(self, cluster):
        """Checks if the cluster has few enough local maxima."""
        return cluster.nlm <= self.nlm
    def is_cluster(self, cluster):
        return self.pass_shape(cluster) and self.pass_time(cluster) \
           and self.pass_ncells(cluster) and self.pass_E(cluster) \
           and self.pass_exotic(cluster) and self.pass_dbc(cluster) \
           and self.pass_nlm(cluster)
    def is_iso(self, cluster, tracks, verbose = False):
        """Determines if cluster is isolated."""
        iso_cone_pt = np.sum([t.pt for t in tracks if cluster.delta_R(t) < self.iso_cone_R])
        if not verbose:
            return iso_cone_pt < self.iso_pt_max
        else:
            return iso_cone_pt < self.iso_pt_max, iso_cone_pt
    def selects(self, cluster, tracks, do_iso = False):
        cluster_check = self.is_cluster(cluster)
        # cluster that does not match any track will return True for track_check
        track_check = all([not self.is_matched(cluster, t) for t in tracks])
        if not do_iso:
            return cluster_check and track_check
        else:
            return cluster_check and track_check and self.is_iso(cluster, tracks)

class TrackSelector(Selector):
    _defaults = {
        'pt_min': .150, # GeV
        'tracksel': sel.TrackSel.globalTrack
    }
    _nosel = {
        'pt_min': -1, # GeV
        'tracksel': ~sel.TrackSel(0)
    }
    _selfid = 'Track selector'

    def _update_selection(self, **selections):
        if "pt_min" in selections:
            self.pt_min = selections["pt_min"]
        if "tracksel" in selections:
            self.tracksel = sel.TrackSel(selections["tracksel"])

    def apply_to(self, df):
        mask = (f"(track_data_pt >= {self.pt_min}) "
             f"&& (track_data_tracksel & {self.tracksel}) "
             )
        ret = df.Define('track_mask', mask)
        for colname, newname in {"track_data_eta": "eta",
                                 "track_data_phi": "phi",
                                 "track_data_pt":  "pt"
                                 }.items():
            ret = ret.Define(newname, f"{colname}[track_mask]")
        return ret

    def selects(self, track, verbose = True):
        tracksel_pass = bool(self.tracksel & track.tracksel)
        pt_min_pass = track.pt > self.pt_min
        if verbose:
            return tracksel_pass and pt_min_pass, tracksel_pass, pt_min_pass
        else:
            return tracksel_pass and pt_min_pass

class AnalysisSelector:
    def __init__(self, event = None, track = None, cluster = None):
        self.event = EventSelector()
        self.track = TrackSelector()
        self.cluster = ClusterSelector()

    @classmethod
    def from_file(cls, file: str):
        selector = cls()
        if file is not None:
            with open(file) as stream:
                cfg = yaml.safe_load(stream)
        selections = cfg['selections']
        use_defaults = cfg.get('defaults', False)
        if "event" in selections:
            selector.event.set_selection(use_defaults = use_defaults, **selections["event"])
        if "track" in selections:
            selector.track.set_selection(use_defaults = use_defaults, **selections["track"])
        if "cluster" in selections:
            selector.cluster.set_selection(use_defaults = use_defaults, **selections["cluster"])
        return selector

    @property
    def active(self):
        sels = ['event', 'track', 'cluster']
        return [sel for sel in sels if getattr(self, sel) is not None]

    def dump(self):
        [getattr(self, sel).dump() for sel in self.active]