import numpy as np
from collections import namedtuple as nt

from alian.analysis.base.selection import TrackSel

_fields = ['pt', 'eta', 'phi', 'label', 'tracksel']

class Track(nt('Track', _fields)):
    __slots__ = ()
    B = 0.5 # signed B field, + for positive polarity and vice versa
    R = 4.5 # reference radius in meters, TPC = 1.1, EMCal = 4.5
    @property
    def p(self):
        """Return track total momentum."""
        return self.pt * np.cosh(self.eta)
    @property
    def phistar(self):
        """Return angle propagated to given reference radius."""
        return (self.phi - self.q * np.arcsin(0.15 * Track.B * Track.R / self.pt)) % (2 * np.pi)
    @property
    def q(self):
        """Return charge."""
        return 1 if self.tracksel & TrackSel.trackSign else -1
    def delta_phi(self, other):
        """Return difference in phi (other - self) in range -pi to pi."""
        dphi = other.phi - self.phi
        if dphi > np.pi:
            return dphi - 2*np.pi
        if dphi < -np.pi:
            return dphi + 2*np.pi
    def delta_eta(self, other):
        """Return difference in eta (other - self)."""
        return other.eta - self.eta
    def delta_R(self, other):
        """Return angular distance in eta-phi."""
        return np.sqrt(self.delta_phi(other) ** 2 + self.delta_eta(other) ** 2)

def get_tracks(ev):
    return [Track(pt, eta, phi, label, TrackSel(tracksel))
              for pt, eta, phi, label, tracksel
              in zip(ev.data['track_data_pt'],
                     ev.data['track_data_eta'],
                     ev.data['track_data_phi'],
                     ev.data['track_data_label'],
                     ev.data['track_data_tracksel'])]