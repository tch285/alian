import alian
import heppyy
import numpy as np
import yaml
from copy import deepcopy
fj = heppyy.load_cppyy('fastjet')

class JetFinder:
    _defaults = {
        'R': 0.4,
        'alg': fj.antikt_algorithm,
        'eta_max': None, # effectively 0.9 - R
        'pT_min': 10,
    }
    def __init__(self, alg = fj.antikt_algorithm, R = 0.4, eta_max = None, pT_min = 10):
        self.R = R
        if eta_max is None:
            self.eta_max = 0.9 - self.R
        self.jet_def = fj.JetDefinition(self.alg, self.R)
        self.jet_selector = fj.SelectorAbsEtaMax(self.eta_max) & fj.SelectorPtMin(pT_min)
        # self.area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.eta_max + self.jet_R))
        # self.bg_estimator = fj.GridMedianBackgroundEstimator(self.bg_y_max, self.bg_grid_spacing)

    @classmethod
    def from_file(cls, file: str):
        options = deepcopy(cls._defaults)
        with open(file) as stream:
            cfg = yaml.safe_load(stream)['jet_finder']
        if 'alg' in cfg:
            if cfg['alg'] in ['antikt', 'anti-kt', 'antikT', 'anti-kT']:
                cfg['alg'] = fj.antikt_algorithm
            elif cfg['alg'] in ['kt', 'kT']:
                cfg['alg'] = fj.kt_algorithm
            elif cfg['alg'] in ['ca', 'CA', 'cambridge' 'cambridge-aachen', 'cambridge/aachen']:
                cfg['alg'] = fj.cambridge_algorithm
            else:
                print(f"Did not recognize jet algorithm '{cfg['alg']}', defaulting to anti-kT.")
                cfg['alg'] = fj.antikt_algorithm
        options.update(cfg)

        return cls(**options)

    def find_jets(self, tracks):
        self.get_psj_tracks(tracks)
        self.get_jets(tracks)

    def get_psj_tracks(self, tracks, m = 0.13957, index_offset = 0):
        pt = np.array([track.pt for track in tracks])
        eta = np.array([track.eta for track in tracks])
        phi = np.array([track.phi for track in tracks])
        self.tracks = alian.numpy_ptetaphi_to_pseudojets(pt, eta, phi, m, index_offset)

    def get_jets(self):
        ca = fj.ClusterSequence(self.tracks, self.jet_def)
        self.jets = fj.sorted_by_pt(self.jet_selector(ca.inclusive_jets()))
