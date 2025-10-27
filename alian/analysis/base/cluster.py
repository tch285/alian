from collections import namedtuple as nt
import numpy as np

_fields = ["energy", "eta", "phi", "m02", "m20", "ncells", "time", "isexotic", "distancebadchannel", "nlm", "clusterdef", "matchedTrackIndex"]

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


def get_clusters(ev):
    arrays = [ev.data[f'cluster_data_{field}'] for field in _fields]
    return [Cluster(*values) for values in zip(*arrays)]