import argparse
import yaml
import uproot

import math
from yasp import GenericObject
from tqdm import tqdm
from root_output import SingleRootFile
from analysis import BaseAnalysis

import ROOT

import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')

def psjv_run3_slow(event, m=0.13957):
    rv = std.vector[fj.PseudoJet]()
    for i in range(len(event['track_data_pt'])):
        px = event['track_data_pt'][i] * math.cos(event['track_data_phi'][i])
        py = event['track_data_pt'][i] * math.sin(event['track_data_phi'][i])
        pz = event['track_data_pt'][i] * math.sinh(event['track_data_eta'][i])
        E = math.sqrt(px * px + py * py + pz * pz + m * m)
        psj = fj.PseudoJet(px, py, pz, E)
        # note that this is wrong because here we treat eta==y
        # psj.reset_PtYPhiM(event['track_data_pt'][i], event['track_data_eta'][i], event['track_data_phi'][i], m)
        psj.set_user_index(i)
        rv.push_back(psj)
    return rv

alian = heppyy.load_cppyy("alian")
# def psjv_run3_cpp(event, m=0.13957):
def psjv_run3(event, m=0.13957):
    rv = alian.numpy_ptetaphi_to_pseudojets(event['track_data_pt'], event['track_data_eta'], event['track_data_phi'], m)
    return rv

# note default is Run3 data
class JetAnalysisRoot(BaseAnalysis):
    _defaults = { 'jet_R': 0.4, 
                  'jet_algorithm': fj.antikt_algorithm, 
                  'jet_eta_max': 0.5,
                  'bg_y_max': 0.9,
                  'bg_grid_spacing': 0.1,
                  'data_input_name': 'aliceRun3',
                  'centrality_min': -1,
                  'centrality_max': 1001,
                  'n_accepted_jets': 0,
                  'n_accepted_events': 0
                  }
    
    def user_init(self):
        self.root_output = SingleRootFile()
        self.tn_mult = ROOT.TNtuple("jet_ev", "jet_ev", "mult:track_count:centr:jet_count:bgrho:bgsigma")
        self.tn_jets = ROOT.TNtuple("jet_v", "jet_v","emult:track_count:centr:pt:eta:phi:m:e:jmult:nlead:leadpt:area:rho")
        self.root_output.add(self.tn_mult)
        self.root_output.add(self.tn_jets)
        self.results.append(self.root_output)

        fj.ClusterSequence().print_banner()
        self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
        self.area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.jet_eta_max + self.jet_R, 1, 0.01))
        self.jet_selector = fj.SelectorAbsEtaMax(self.jet_eta_max)
        self.bg_estimator = fj.GridMedianBackgroundEstimator(self.bg_y_max, self.bg_grid_spacing)
                
    def set_centrality(self, centmin, centmax):
        self.centrality_min = centmin
        self.centrality_max = centmax

    def build_psjv(self, event_struct):
        if event_struct.psjv is None:
          event_struct.psjv = GenericObject()
        if event_struct.psjv.run3 is None:
          event_struct.psjv.run3 = psjv_run3(event_struct.data[self.data_input_name])
    
    def accept_event(self, event_struct):
        if event_struct is None:
            return False
        if not isinstance(event_struct, GenericObject):
            return False
        if self.data_input_name not in event_struct.data:
            return False
        event = event_struct.data[self.data_input_name]
        if event is None:
            return False
        if 'multiplicity' not in event:
            return False
        if 'centrality' not in event:
            return False
        if 'track_data_pt' not in event:
            return False
        if event['centrality'] < self.centrality_min or event['centrality'] > self.centrality_max:
            return False
        return True
      
    def process_event(self, event_struct):
        if not self.accept_event(event_struct):
            return
        event = event_struct.data[self.data_input_name]
        self.multiplicity = event['multiplicity']
        self.centrality = event['centrality']
        self.n_accepted_events += 1
        self.track_count = len(event['track_data_pt'])
        self.build_psjv(event_struct)
        self.psjv = event_struct.psjv.run3
        self.jet_analysis(event_struct)

    def jet_analysis(self, event_struct):
        # estimate event background rho with grid estimator
        self.bg_estimator.set_particles(self.psjv)
        rho = self.bg_estimator.rho()
        sigma = self.bg_estimator.sigma()
        # print('rho:', rho, 'centrality:', self.centrality, 'track_count:', self.track_count, 'self.psjv.size():', self.psjv.size())
        self.ca = fj.ClusterSequenceArea(self.psjv, self.jet_def, self.area_def)
        self.jets = fj.sorted_by_pt(self.jet_selector(self.ca.inclusive_jets()))
        jet_count = 0
        for ij, j in enumerate(self.jets):
            jconstits = [jconst for jconst in j.constituents() if jconst.perp() > 0.1]
            if len(jconstits) == 0:
                continue
            leadpt = fj.sorted_by_pt(jconstits)[0].perp()
            self.tn_jets.Fill(self.multiplicity, self.track_count, self.centrality, j.perp(), j.eta(), j.phi(), j.m(), j.e(), len(jconstits), ij, leadpt, j.area(), rho)
            jet_count += 1
            self.n_accepted_jets += 1
        self.tn_mult.Fill(self.multiplicity, self.track_count, self.centrality, jet_count, rho, sigma)
        

def psjv_run2(event, m=0.13957):
    rv = std.vector[fj.PseudoJet]()
    for i in range(len(event['ParticlePt'])):
        px = event['ParticlePt'][i] * math.cos(event['ParticlePhi'][i])
        py = event['ParticlePt'][i] * math.sin(event['ParticlePhi'][i])
        pz = event['ParticlePt'][i] * math.sinh(event['ParticleEta'][i])
        E = math.sqrt(px * px + py * py + pz * pz + m * m)
        psj = fj.PseudoJet(px, py, pz, E)
        # note that this is wrong because here we treat eta==y
        # psj.reset_PtYPhiM(event['ParticlePt'][i], event['ParticleEta'][i], event['ParticlePhi'][i], m)
        psj.set_user_index(i)
        rv.push_back(psj)
    return rv


class JetAnalysisRootRun2(JetAnalysisRoot):
    _defaults = JetAnalysisRoot._defaults.copy()
    _defaults['data_input_name'] = 'aliceRun2'
    
    def user_init(self):
        # call the base class user_init
        super(JetAnalysisRootRun2, self).user_init()
                
    def build_psjv(self, event_struct):
        if event_struct.psjv is None:
          event_struct.psjv = GenericObject()
        if event_struct.psjv.run2 is None:
          event_struct.psjv.run2 = psjv_run2(event_struct.data[self.data_input_name])
      
    def accept_event(self, event_struct):
        if event_struct is None:
            return False
        if not isinstance(event_struct, GenericObject):
            return False
        if self.data_input_name not in event_struct.data:
            return False
        event = event_struct.data[self.data_input_name]
        if event is None:
            return False
        if 'V0Amult' not in event:
            return False
        if 'centrality' not in event:
            return False
        if 'ParticlePt' not in event:
            return False
        if list(set(event['centrality']))[0] < self.centrality_min or list(set(event['centrality']))[0] > self.centrality_max:
            return False
        return True

    def process_event(self, event_struct):
        if not self.accept_event(event_struct):
            return
        event = event_struct.data[self.data_input_name]
        self.multiplicity = list(set(event['V0Amult']))[0]
        self.centrality = list(set(event['centrality']))[0]
        self.track_count = len(event['ParticlePt'])
        self.n_accepted_events += 1
        self.build_psjv(event_struct)
        self.psjv = event_struct.psjv.run2
        self.jet_analysis(event_struct)
