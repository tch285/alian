from collections import namedtuple as nt
import numpy as np

_fields = ["energy",
           "eta",
           "phi",
           "m02",
           "m20",
           "ncells",
           "time",
           "isexotic",
           "distancebadchannel",
           "nlm",
           "clusterdef",
           "matchedTrackN",
           "matchedTrackEta",
           "matchedTrackPhi",
           "matchedTrackP"
           ]

class Cluster(nt('Cluster', _fields)):
    __slots__ = ()
    def delta_phi(self, other):
        dphi = other.phi - self.phi
        if dphi > np.pi:
            return dphi - 2*np.pi
        if dphi < -np.pi:
            return dphi + 2*np.pi
        return dphi
    def delta_eta(self, other):
        return other.eta - self.eta
    def delta_R(self, other):
        return np.sqrt(self.delta_phi(other) ** 2 + self.delta_eta(other) ** 2)
    def has_geo_ep_matched_track(self, deta, dphi, ep):
        for track_eta, track_phi, track_p in zip(self.matchedTrackEta, self.matchedTrackPhi, self.matchedTrackP):
            if np.abs(track_eta - self.eta) < deta and np.abs(track_phi - self.phi) < dphi and self.energy / track_p < ep:
                return True
        return False
    def has_geo_matched_track(self, deta, dphi):
        for track_eta, track_phi in zip(self.matchedTrackEta, self.matchedTrackPhi):
            if np.abs(track_eta - self.eta) < deta and np.abs(track_phi - self.phi) < dphi:
                return True
        return False
    def has_ep_matched_track(self, ep):
        for track_p in self.matchedTrackP:
            if self.energy / track_p < ep:
                return True
        return False

def get_clusters(ev):
    arrays = [ev.data[f'cluster_data_{field}'] for field in _fields]
    return [Cluster(*values) for values in zip(*arrays)]