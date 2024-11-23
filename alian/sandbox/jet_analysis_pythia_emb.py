import math
from yasp import GenericObject
from tqdm import tqdm
from root_output import SingleRootFile
from jet_analysis import JetAnalysisRoot
from pythia_io import psj_from_particle_with_index
import numba

import ROOT

import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')


class JetAnalysisPythiaEmbRoot(JetAnalysisRoot):
    _defaults = { 'jet_R': 0.4, 
                  'jet_algorithm': fj.antikt_algorithm, 
                  'jet_eta_max': 0.5,
                  'jet_pt_min': 10,
                  'bg_y_max': 0.9,
                  'bg_grid_spacing': 0.1,
                  'data_input_name': 'aliceRun3',
                  'centrality_min': -1,
                  'centrality_max': 1001,
                  'n_accepted_jets': 0,
                  'n_accepted_events': 0
                  }
    
    def init(self):
        self.root_output = SingleRootFile()
        self.tn_mult = ROOT.TNtuple("emb_jet_ev", "emb_jet_ev", "mult:track_count:centr:jet_count:bgrho:bgsigma")
        # self.tn_jets = ROOT.TNtuple("emb_jet_v", "emb_jet_v","emult:track_count:centr:pt:eta:phi:m:e:jmult:nlead:leadpt:area:rho:pt_sig:eta_sig:phi_sig:m_sig:e_sig:jmult_sig:nlead_sig:leadpt_sig:area_sig")
        self.tn_jets = ROOT.TNtuple("emb_jet_v", "emb_jet_v","centr:pt:jmult:nlead:leadpt:area:rho:pt_sig:jmult_sig:nlead_sig:leadpt_sig:area_sig:deltaR:deltapt:deltapt_rho")
        self.root_output.add(self.tn_mult)
        self.root_output.add(self.tn_jets)
        self.results.append(self.root_output)

        fj.ClusterSequence().print_banner()
        self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
        self.area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.jet_eta_max + self.jet_R, 1, 0.01))
        self.jet_selector = fj.SelectorAbsEtaMax(self.jet_eta_max)
        self.jet_selector_pt = fj.SelectorPtMin(self.jet_pt_min)
        self.bg_estimator = fj.GridMedianBackgroundEstimator(self.bg_y_max, self.bg_grid_spacing)
        
        self._build_psjv = self.build_psjv_run3
        self._accept_event = self.accept_event_run3
        if self.lhc_run:
            if self.lhc_run == 2:
                self._build_psjv = self.build_psjv_run2
                self._accept_event = self.accept_event_run2
        if self.pythia_offset_index is None:
          self.pythia_offset_index = 100000 # unlikely a particle from bg event with this index


    def get_pythia_jets_to_embed(self, event_struct, ntries=10000):
        while ntries > 0:
          ntries -= 1
          # pythia_check = event_struct.data['PythiaIO']
          pythiaIO = event_struct.source['PythiaIO']
          # _pythia = pythiaIO.pythia
          _pythia = next(pythiaIO.next_event())
          self.pythia_parts = std.vector[fj.PseudoJet]([psj_from_particle_with_index(p, i + self.pythia_offset_index) 
                                                          for i, p in enumerate(_pythia.event) if p.isFinal() and p.isVisible() and p.isCharged()])

          self.pythia_ca = fj.ClusterSequenceArea(self.pythia_parts, self.jet_def, self.area_def)
          self.pythia_jets = fj.sorted_by_pt(self.jet_selector_pt(self.jet_selector(self.pythia_ca.inclusive_jets())))
            
          # print('pythia_parts.size():', self.pythia_parts.size())
          # print('pythia_jets.size():', self.pythia_jets.size()) 
          if len(self.pythia_jets) > 0:
            return True
        RuntimeError('Failed to find jets to embed after {} tries'.format(ntries))
        return False

    def embed_jet(self, j):
        # print('Embedding jets')
        _rv = std.vector[fj.PseudoJet]()
        for jconst in j.constituents():
            _rv.push_back(jconst)
        for bg in self.psjv:
            _rv.push_back(bg)
        return _rv            

    def pt_match_jets(self, jsig, jemb):
        ptmatch = 0
        for jsigc in jsig.constituents():
            for jembc in jemb.constituents():
                if jsigc.user_index() == jembc.user_index():
                    ptmatch += jsigc.pt()
        return ptmatch / jsig.pt()

    def analysis(self, event_struct):
        if self.get_pythia_jets_to_embed(event_struct) is False:
          RuntimeError('Failed to find jets to embed')
        
        for ijsig, jsig in enumerate(self.pythia_jets): 
          leadptsig = fj.sorted_by_pt(jsig.constituents())[0].perp()
          _parts = self.embed_jet(jsig)
          # print(f'parts in bg: {self.psjv.size()} parts in jet {jsig.constituents().size()}, total: {_parts.size()}')
          # estimate event background rho with grid estimator
          self.bg_estimator.set_particles(_parts)
          rho = self.bg_estimator.rho()
          sigma = self.bg_estimator.sigma()
          # print('rho:', rho, 'centrality:', self.centrality, 'track_count:', self.track_count, 'self.psjv.size():', self.psjv.size())
          self.ca = fj.ClusterSequenceArea(_parts, self.jet_def, self.area_def)
          self.jets = fj.sorted_by_pt(self.jet_selector(self.ca.inclusive_jets()))
          jet_count = 0
          for ij, j in enumerate(self.jets):
              if self.pt_match_jets(jsig, j) < 0.5:
                continue
              leadpt = fj.sorted_by_pt(j.constituents())[0].perp()
              # self.tn_jets = ROOT.TNtuple("emb_jet_v", "emb_jet_v","emult:track_count:centr:pt:eta:phi:m:e:jmult:nlead:leadpt:area:rho:pt_sig:eta_sig:phi_sig:m_sig:e_sig:jmult_sig:nlead_sig:leadpt_sig:area_sig")
              # self.tn_jets.Fill(self.multiplicity, self.track_count, self.centrality, j.perp(), j.eta(), j.phi(), j.m(), j.e(), j.constituents().size(), ij, leadpt, j.area(), rho,
              #                  jsig.perp(), jsig.eta(), jsig.phi(), jsig.m(), jsig.e(), jsig.constituents().size(), ijsig, leadptsig, jsig.area())
              self.tn_jets.Fill(self.centrality, j.perp(), j.constituents().size(), ij, leadpt, j.area(), rho,
                                                jsig.perp(), jsig.constituents().size(), ijsig, leadptsig, jsig.area(), jsig.delta_R(j), j.perp() - jsig.perp(), j.perp() - jsig.perp() - rho * j.area())
              jet_count += 1
              self.n_accepted_jets += 1
              break
          self.tn_mult.Fill(self.multiplicity, self.track_count, self.centrality, -1, rho, sigma)
