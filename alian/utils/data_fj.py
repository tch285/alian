import heppyy
import math

fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')

def psjv_from_tracks_run3_slow(event, index_offset=0, m=0.13957):
    rv = std.vector[fj.PseudoJet]()
    for i in range(len(event['track_data_pt'])):
        px = event['track_data_pt'][i] * math.cos(event['track_data_phi'][i])
        py = event['track_data_pt'][i] * math.sin(event['track_data_phi'][i])
        pz = event['track_data_pt'][i] * math.sinh(event['track_data_eta'][i])
        E = math.sqrt(px * px + py * py + pz * pz + m * m)
        psj = fj.PseudoJet(px, py, pz, E)
        # note that this is wrong because here we treat eta==y
        # psj.reset_PtYPhiM(event['track_data_pt'][i], event['track_data_eta'][i], event['track_data_phi'][i], m)
        psj.set_user_index(i + index_offset)
        rv.push_back(psj)
    return rv

alian = heppyy.load_cppyy("alian")
# def psjv_from_tracks_run3_cpp(event, index_offset=0, m=0.13957):
def psjv_from_tracks_run3(event, index_offset=0, m=0.13957):
    rv = alian.numpy_ptetaphi_to_pseudojets(event['track_data_pt'], event['track_data_eta'], event['track_data_phi'], m, index_offset)
    return rv

def psjv_from_tracks_run2(event, index_offset=0, m=0.13957):
    rv = std.vector[fj.PseudoJet]()
    for i in range(len(event['ParticlePt'])):
        px = event['ParticlePt'][i] * math.cos(event['ParticlePhi'][i])
        py = event['ParticlePt'][i] * math.sin(event['ParticlePhi'][i])
        pz = event['ParticlePt'][i] * math.sinh(event['ParticleEta'][i])
        E = math.sqrt(px * px + py * py + pz * pz + m * m)
        psj = fj.PseudoJet(px, py, pz, E)
        # note that this is wrong because here we treat eta==y
        # psj.reset_PtYPhiM(event['ParticlePt'][i], event['ParticleEta'][i], event['ParticlePhi'][i], m)
        psj.set_user_index(i + index_offset)
        rv.push_back(psj)
    return rv

def data_tracks_to_pseudojets(event, index_offset=0, m=0.13957, lhc_run=3):
		if lhc_run == 3:
				return psjv_from_tracks_run3(event, index_offset, m)
		elif lhc_run == 2:
				return psjv_from_tracks_run2(event, index_offset, m)
		else:
				raise ValueError(f'Invalid run number: {lhc_run}')
		return None
  
