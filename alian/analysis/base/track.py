import numpy as np
from collections import namedtuple as nt

from alian.analysis.base.selections import TrackSel

_fields = ['pt', 'eta', 'phi', 'label', 'tracksel']

class Track(nt('Track', _fields)):
    __slots__ = ()
    @property
    def p(self):
        return self.pt * np.cosh(self.eta)
    def delta_phi(self, track):
        dphi = track.phi - self.phi
        if dphi > np.pi:
            return dphi - 2*np.pi
        if dphi < -np.pi:
            return dphi + 2*np.pi
    def delta_eta(self, track):
        return track.eta - self.eta
    def deltaR(self, track):
        return np.sqrt(self.delta_phi(track) ** 2 + self.delta_eta(track) ** 2)

def get_tracks(ev):
    return [Track(pt, eta, phi, label, TrackSel(tracksel))
              for pt, eta, phi, label, tracksel
              in zip(ev.data['track_data_pt'],
                     ev.data['track_data_eta'],
                     ev.data['track_data_phi'],
                     ev.data['track_data_label'],
                     ev.data['track_data_tracksel'])]