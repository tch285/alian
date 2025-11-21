from alian.analysis.base.track import get_tracks
from alian.analysis.base.cluster import get_clusters
from alian.analysis.base.selection import EvSel, TrgSel

class Event:
    def __init__(self, event_struct):
        self.run_number = event_struct.data['run_number']
        self.event_selection = EvSel(event_struct.data['event_selection'])
        self.triggersel = TrgSel(event_struct.data['triggersel'])
        self.centrality = event_struct.data['centrality']
        self.multiplicity = event_struct.data['multiplicity']
        self.tracks = get_tracks(event_struct)
        self.clusters = get_clusters(event_struct)