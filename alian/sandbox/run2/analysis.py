#!/usr/bin/env python

import argparse
import yaml
import uproot
import pandas as pd

from yasp import GenericObject
from tqdm import tqdm

import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')

from alian.sandbox.root_output import SingleRootFile

import ROOT

class BaseAnalysis(GenericObject):
    _defaults = {}
    
    def __init__(self, **kwargs):
        super(BaseAnalysis, self).__init__(**kwargs)
        self.results = []
        for k, val in self.__class__._defaults.items():
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, val)
                # print(f'[i] Setting {k} to default value {val}', self.__class__, self.__getattr__(k))
        self.user_init()
        print(self)

    def user_init(self):
        # Implement the user-defined initialization logic here
        pass
      
    def process_event(self, event):
        # Implement the event processing logic here
        pass

    def write_to_file(self):
        if self.output_file is not None:
            with open(self.output_file, 'w') as file:
                # Implement the logic to write results to file
                file.write(str(self.results))

def charged_particles_psjv(event, m=0.13957):
    rv = std.vector[fj.PseudoJet]()
    for i in range(len(event['ParticlePt'])):
      psj = fj.PseudoJet()
      psj.reset_PtYPhiM(event['ParticlePt'][i], event['ParticleEta'][i], event['ParticlePhi'][i], m)
      psj.set_user_index(i)
      rv.push_back(psj)
    return rv

class JetAnalysisRoot(BaseAnalysis):
    _defaults = {'jet_R': 0.4, 
                  'jet_algorithm': fj.antikt_algorithm, 
                  'jet_eta_max': 0.5
                  }
    
    def user_init(self):
        print('defaults are:', self.__class__._defaults)
        self.root_output = SingleRootFile()
        self.tn_mult = ROOT.TNtuple("jet_ev", "jet_ev", "mult:track_count:centr:jet_count")
        self.tn_jets = ROOT.TNtuple("jet_v", "jet_v","emult:track_count:centr:pt:eta:phi:m:e:jmult:nlead:leadpt:area:rho")
        self.root_output.add(self.tn_mult)
        self.root_output.add(self.tn_jets)
        self.results.append(self.root_output)

        fj.ClusterSequence().print_banner()
        self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
        # self.ca = fj.ClusterSequenceArea(charged_particles_psjv(event), self.jet_def, fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.jet_eta_max)))
        self.area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.jet_eta_max + self.jet_R, 1, 0.01))
        self.jet_selector = fj.SelectorAbsEtaMax(self.jet_eta_max)
        self.bg_estimator = fj.GridMedianBackgroundEstimator(1.0, 0.1)
        
    def process_event(self, event):
        #multiplicity = event['multiplicity']
        multiplicity = list(set(event['V0Amult']))[0]
        centrality = list(set(event['centrality']))[0]
        # print("Event:", multiplicity, centrality)
        
        # track_count = len(event['track_data_pt'])
        track_count = len(event['ParticlePt'])
        # if track_count < 10:
        #     return
        # estimate event background rho with grid estimator
        self.bg_estimator.set_particles(charged_particles_psjv(event))
        rho = self.bg_estimator.rho()
        if rho < 20.:
            return
        psjv = charged_particles_psjv(event)
        # jets = fj.sorted_by_pt(self.jet_selector(self.jet_def(psjv)))        
        ca = fj.ClusterSequenceArea(psjv, self.jet_def, self.area_def)
        jets = fj.sorted_by_pt(self.jet_selector(ca.inclusive_jets()))
        jet_count = 0
        for ij, j in enumerate(jets):
            jconstits = [jconst for jconst in j.constituents() if jconst.perp() > 0.1]
            if len(jconstits) == 0:
                continue
            leadpt = fj.sorted_by_pt(jconstits)[0].perp()
            self.tn_jets.Fill(multiplicity, track_count, centrality, j.perp(), j.eta(), j.phi(), j.m(), j.e(), len(jconstits), ij, leadpt, j.area(), rho)
            jet_count += 1
        self.tn_mult.Fill(multiplicity, track_count, centrality, jet_count)
        
        
class AnalysisEngine:
    def __init__(self, root_file_path, tree_particle_name, tree_event_char_name, branches, **kwargs):
        self.root_file_path = root_file_path
        self.tree_particle_name = tree_particle_name
        self.tree_event_char_name = tree_event_char_name
        self.branches = branches
        self.analyses = []

    def add_analysis(self, analysis):
        self.analyses.append(analysis)

    def run(self):
        # Open the ROOT file
        file = uproot.open(self.root_file_path)
        # Access the trees
        tree_particle = file[self.tree_particle_name]
        tree_event_char = file[self.tree_event_char_name]
        
        # Convert the trees to pandas DataFrames
        df_particle = tree_particle.arrays(library="pd")
        df_event_char = tree_event_char.arrays(library="pd")
        
        # Merge the DataFrames on 'ev_id', 'ev_id_ext', and 'run_number'
        merged_df = pd.merge(df_particle, df_event_char, on=['ev_id', 'ev_id_ext', 'run_number'])
        
        # Filter out events where 'is_ev_rej' is not zero and 'centrality' is not less than 10
        filtered_df = merged_df[(merged_df['is_ev_rej'] == 0) & (merged_df['centrality'] < 10)]
        
        # Group by 'ev_id', 'ev_id_ext', and 'run_number'
        grouped_df = filtered_df.groupby(['ev_id', 'ev_id_ext', 'run_number'])
        
        # Efficiently iterate over the grouped DataFrame with a progress bar
        total_entries = len(grouped_df)
        with tqdm(total=total_entries, desc="Processing events") as pbar:
            for name, group in grouped_df:
                event = {branch: group[branch].tolist() for branch in self.branches if branch in group}
                self.process_event(event)
                pbar.update(1)

    def process_event(self, event):
        for analysis in self.analyses:
            analysis.process_event(event)

    def write_results(self):
        for analysis in self.analyses:
            print(f'Writing results for {analysis.output_file}')
            analysis.write_to_file()

def load_default_trees_and_branches(yaml_file_path):
    with open(yaml_file_path, 'r') as file:
        config = yaml.safe_load(file)
    
    tree_particle_name = 'PWGHF_TreeCreator/tree_Particle'
    tree_event_char_name = 'PWGHF_TreeCreator/tree_event_char'
    
    branches_particle = list(config['PWGHF_TreeCreator']['tree_Particle'].keys())
    branches_event_char = list(config['PWGHF_TreeCreator']['tree_event_char'].keys())
    
    branches = branches_particle + branches_event_char
    
    return tree_particle_name, tree_event_char_name, branches

if __name__ == '__main__':
    # Example usage
    import argparse

    parser = argparse.ArgumentParser(description='Run the analysis engine.')
    parser.add_argument('root_file_path', help='Path to the ROOT file')
    parser.add_argument('yaml_file_path', help='Path to the YAML file')
    parser.add_argument('-o', '--output', help='Output file name', default='analysis_results.root')
    args = parser.parse_args()

    tree_particle_name, tree_event_char_name, branches = load_default_trees_and_branches(args.yaml_file_path)

    engine = AnalysisEngine(args.root_file_path, tree_particle_name, tree_event_char_name, branches)
    # Add your analyses to the engine here
    # engine.add_analysis(your_analysis)
    root_output = SingleRootFile(args.output)
    jet_analysis = JetAnalysisRoot()
    engine.add_analysis(jet_analysis)
    engine.run()
    engine.write_results()
    root_output.close()