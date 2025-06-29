class Cluster:
    def __str__(self):
        return (f"Cluster(energy={self.energy}, "
                f"eta={self.eta}, "
                f"phi={self.phi}, "
                f"m02={self.m02}, "
                f"m20={self.m20}, "
                f"ncells={self.ncells})")

    def __init__(self, energy = 0, eta = 0, phi = 0, m02 = 0, m20 = 0, ncells = 0, time = 0, isexotic = 0, distancebadchannel = 0, nlm = 0, clusterdef = 0, matchedTrackIndex = 0):
        self.energy = energy
        self.eta = eta
        self.phi = phi
        self.m02 = m02
        self.m20 = m20
        self.ncells = ncells
        self.time = time
        self.isexotic = isexotic
        self.distancebadchannel = distancebadchannel
        self.nlm = nlm
        self.clusterdef = clusterdef
        self.matchedTrackIndex = matchedTrackIndex

    def get_event_clusters(cls, ev_struct):
        return [cls(energy, eta, phi, m02, m20, ncells, time, isexotic, distancebadchannel, nlm, clusterdef, matchedTrackIndex)
                for energy, eta, phi, m02, m20, ncells, time, isexotic, distancebadchannel, nlm, clusterdef, matchedTrackIndex
                in zip(ev_struct.data['cluster_data_energy'],
                       ev_struct.data['cluster_data_eta'],
                       ev_struct.data['cluster_data_phi'],
                       ev_struct.data['cluster_data_m02'],
                       ev_struct.data['cluster_data_m20'],
                       ev_struct.data['cluster_data_ncells'],
                       ev_struct.data['cluster_data_time'],
                       ev_struct.data['cluster_data_isexotic'],
                       ev_struct.data['cluster_data_distancebadchannel'],
                       ev_struct.data['cluster_data_nlm'],
                       ev_struct.data['cluster_data_clusterdef'],
                       ev_struct.data['cluster_data_matchedTrackIndex'])]