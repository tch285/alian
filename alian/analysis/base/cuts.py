import yaml

import alian.analysis.base.selections as sel

class EventCuts:
    def __init__(self,
                 ev_sel: sel.EvSel = sel.EvSel.sel8,
                 trg_sel: sel.TrgSel = sel.TrgSel.GammaHighPtEMCAL
                 ):
        self.ev_sel = ev_sel
        self.trg_sel = trg_sel

    def set_cuts(self, cuts: dict):
        if "ev_sel" in cuts:
            self.ev_sel = sel.EvSel(0)
            for esel in cuts["ev_sel"].split(","):
                self.ev_sel |= sel.EvSel[esel]
        if "trg_sel" in cuts:
            self.trg_sel = sel.TrgSel(0)
            for tsel in cuts["trg_sel"].split(","):
                self.trg_sel |= sel.TrgSel[tsel]

    def dump(self):
        print( "Event cuts:\n"
              f"\tEvent selection: {self.ev_sel}\n"
              f"\tTrigger selection: {self.trg_sel}")

class ClusterCuts:
    def __init__(self, **cuts):
        _defaults = {
            "E_min": 0.7, # GeV
            "m02_min": 0.1,
            "m02_max": 0.3,
            "time_min": -30, # ns
            "time_max": 35, # ns
            "ncells_min": 2,
            "delta_phi": 0.05,
            "delta_eta": 0.05,
            "energy_pt_ratio_threshold": 1.75,
            "iso_cone_radius": 0.4,
            "iso_pt_threshold": 1.5, # GeV
        }
        self.cut_names = _defaults.keys()
        valid_cuts = self._validate_cut_names(cuts)
        _defaults.update(valid_cuts)

        for name, val in _defaults.items():
            setattr(self, name, val)

    def _validate_cut_names(self, cuts: dict):
        invalid_cut_names = [name for name in cuts.keys() if name not in self.cut_names]
        if invalid_cut_names:
            print(f"Invalid cluster cuts given, will be ignored: {', '.join(invalid_cut_names)}")
        return {name: val for name, val in cuts.items() if name in self.cut_names}

    def update_cuts(self, cuts: dict):
        valid_cuts = self._validate_cut_names(cuts)

        for name, val in valid_cuts.items():
            setattr(self, name, val)

    def dump(self):
        print(   "Cluster cuts:\n"
              f"\tE_min: {self.E_min} GeV\n"
              f"\tm02_min: {self.m02_min}\n"
              f"\tm02_max: {self.m02_max}\n"
              f"\ttime_min: {self.time_min} ns\n"
              f"\ttime_max: {self.time_max} ns\n"
              f"\tncells_min: {self.ncells_min}\n"
              f"\tdelta_phi: {self.delta_phi}\n"
              f"\tdelta_eta: {self.delta_eta}\n"
              f"\tenergy_pt_ratio_threshold: {self.energy_pt_ratio_threshold}\n"
              f"\tiso_cone_radius: {self.isolation_cone_radius}\n"
              f"\tiso_pt_threshold: {self.isolation_pt_threshold} GeV" )

class TrackCuts:
    def __init__(self, min_pt: float = 0.150, tracksel: sel.TrackSel = sel.TrackSel.globalTrack):
        self.min_pt = min_pt # GeV
        self.tracksel = tracksel

    def update_cuts(self, cuts: dict):
        # simpler validation
        if "min_pt" in cuts:
            self.min_pt = cuts["min_pt"]
        if "tracksel" in cuts:
            self.tracksel = cuts['tracksel']

    def update_tracksel(self, newsel):
        # to add a selection rather than override it
        self.tracksel |= newsel

    def dump(self):
        print(    "Track cuts:\n"
              f"\tTrack selection: {self.tracksel}\n"
              f"\tMin pT: {self.min_pt}")


class AnalysisCuts:
    def __init__(self, file):
        self.extract_from_file(file)

    def extract_from_file(self, file: str):
        try:
            with open(file) as stream:
                yaml_struct = yaml.safe_load(stream)
                if "event_cuts" in yaml_struct:
                    self.event.set_cuts(yaml_struct["event_cuts"])
                else:
                    self.event = EventCuts()
                if "track_cuts" in yaml_struct:
                    self.track.set_cuts(yaml_struct["track_cuts"])
                else:
                    self.track = TrackCuts()
                if "cluster_cuts" in yaml_struct:
                    self.cluster.set_cuts(yaml_struct["cluster_cuts"])
                else:
                    self.cluster = ClusterCuts()
        except yaml.YAMLError:
            print(f"Error processing file {file}")
            exit(1)
        except FileNotFoundError:
            print(f"Could not find cuts YAML file at {file}")
            exit(1)

    def dump(self):
        self.event.dump()
        self.track.dump()
        self.cluster.dump()