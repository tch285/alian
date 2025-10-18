from collections import namedtuple as nt

_fields = ["energy", "eta", "phi", "m02", "m20", "ncells", "time", "isexotic", "distancebadchannel", "nlm", "clusterdef", "matchedTrackIndex"]

Cluster = nt('Cluster', _fields)

def get_clusters(ev):
    arrays = [ev.data[f'cluster_data_{field}'] for field in _fields]
    return [Cluster(*values) for values in zip(*arrays)]
    # return [Cluster(energy, eta, phi, m02, m20, ncells, time, isexotic, distancebadchannel, nlm, clusterdef, matchedTrackIndex)
    #             for energy, eta, phi, m02, m20, ncells, time, isexotic, distancebadchannel, nlm, clusterdef, matchedTrackIndex
    #             in zip(ev.data['cluster_data_energy'],
    #                    ev.data['cluster_data_eta'],
    #                    ev.data['cluster_data_phi'],
    #                    ev.data['cluster_data_m02'],
    #                    ev.data['cluster_data_m20'],
    #                    ev.data['cluster_data_ncells'],
    #                    ev.data['cluster_data_time'],
    #                    ev.data['cluster_data_isexotic'],
    #                    ev.data['cluster_data_distancebadchannel'],
    #                    ev.data['cluster_data_nlm'],
    #                    ev.data['cluster_data_clusterdef'],
    #                    ev.data['cluster_data_matchedTrackIndex'])]