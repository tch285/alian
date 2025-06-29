import numpy as np

class Track:
    def __str__(self):
        return f"Track(pT={self.pt}, eta={self.eta}, phi={self.phi})"

    def __init__(self, pt = 0, eta = 0, phi = 0, label = 0, tracksel = 0):
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.label = label
        self.tracksel = tracksel
        self.p = self.pt * np.cosh(self.eta)
    
    def delta_phi(self, track):
        dphi = track.phi - self.phi
        if dphi > np.pi:
            return dphi - 2*np.pi
        if dphi < -np.pi:
            return dphi + 2*np.pi
        # return (track.phi - self.phi + np.pi) % (2 * np.pi) - np.pi

    def delta_eta(self, track):
        return track.eta - self.eta

    def deltaR(self, track):
        return np.sqrt(self.delta_phi(track) ** 2 + self.delta_eta(track) ** 2)

    @classmethod
    def get_event_tracks(cls, ev_struct):
        return [cls(pt, eta, phi, label, tracksel) for pt, eta, phi, label, tracksel
                in zip(ev_struct.data['track_data_pt'],
                       ev_struct.data['track_data_eta'],
                       ev_struct.data['track_data_phi'],
                       ev_struct.data['track_data_label'],
                       ev_struct.data['track_data_tracksel'])]

class TrackSelector:
    