from alian.analysis.base.track import Track
from alian.analysis.base.cluster import Cluster

class Event:
    def __init__(self, event_struct):
        self.run_number = event_struct.data['run_number']
        self.event_selection = event_struct.data['event_selection']
        self.triggersel = event_struct.data['triggersel']
        self.centrality = event_struct.data['centrality']
        self.multiplicity = event_struct.data['multiplicity']
        self.tracks = Track.get_event_tracks(event_struct)
        self.clusters = Cluster.get_event_clusters(event_struct)