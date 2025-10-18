import yaml

import alian.analysis.base.selections as sel

class Selector:
    _defaults = {}
    _selfid = ''
    def __init__(self, **selections):
        self.set_selection(**selections)

    def set_selection(self, **selections):
        """Set selections. Unspecified selections are set to default values."""
        self.reset_selection()
        self.update_selection(**selections)

    def reset_selection(self):
        """Reset selection to defaults."""
        for selection, value in self._defaults.items():
            setattr(self, selection, value)

    def update_selection(self, **selections):
        valid_sels = self._validate(**selections)
        self._update_selection(**valid_sels)
    def _update_selection(self, **selections):
        raise NotImplementedError

    def _validate(self, **selections):
        bad_sels = [sel for sel in selections.keys() if sel not in self._defaults]
        if bad_sels:
            print(f"{self._selfid}: invalid cuts given, will be ignored: {', '.join(bad_sels)}")
        return {sel: val for sel, val in selections.items() if sel in self._defaults}

    def dump(self):
        print(f"{self._selfid}:\n" + "\n".join([f"\t{key}: {getattr(self, key)}" for key in self._defaults]))


class EventSelector(Selector):
    _defaults = {
        "event_selection": sel.EvSel.sel8,
        "triggersel": sel.TrgSel.GammaHighPtEMCAL
    }
    _selfid = "Event selector"

    def _update_selection(self, **selections):
        if "event_selection" in selections:
            self.event_selection = sel.EvSel(selections["event_selection"])
        if "triggersel" in selections:
            self.triggersel = sel.TrgSel(selections["triggersel"])

class ClusterSelector(Selector):
    _defaults = {
        "E_min": 0.7, # GeV
        "m02_min": 0.1,
        "m02_max": 0.3,
        "time_min": -30, # ns
        "time_max": 35, # ns
        "ncells_min": 2,
        "delta_phi": 0.05, # radians
        "delta_eta": 0.05, # pseudorapidity
        "energy_pt_ratio_threshold": 1.75,
        "iso_cone_radius": 0.4,
        "iso_pt_threshold": 1.5, # total GeV
    }
    _selfid = 'Cluster selector'

    def _update_selection(self, **selections):
        for name, val in selections.items():
            setattr(self, name, val)

class TrackSelector(Selector):
    _defaults = {
        'min_pt': .150, # GeV
        'tracksel': sel.TrackSel.globalTrack
    }
    _selfid = 'Track selector'

    def _update_selection(self, **selections):
        if "min_pt" in selections:
            self.min_pt = selections["min_pt"]
        if "tracksel" in selections:
            self.tracksel = sel.TrackSel(selections["tracksel"])

class AnalysisSelector:
    def __init__(self, event = None, track = None, cluster = None):
        self.event = event
        self.track = track
        self.cluster = cluster

    @classmethod
    def from_file(cls, file: str):
        selector = cls()
        with open(file) as stream:
            cuts = yaml.safe_load(stream)
        if "event_cuts" in cuts:
            selector.event = EventSelector(**cuts["event_cuts"])
        if "track_cuts" in cuts:
            selector.track = TrackSelector(**cuts["track_cuts"])
        if "cluster_cuts" in cuts:
            selector.cluster = ClusterSelector(**cuts["cluster_cuts"])
        return selector

    @property
    def active(self):
        sels = ['event', 'track', 'cluster']
        return [sel for sel in sels if getattr(self, sel) is not None]

    def dump(self):
        [getattr(self, sel).dump() for sel in self.active]