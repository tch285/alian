#!/usr/bin/env python

import argparse
import yaml
import uproot

from yasp import GenericObject
from tqdm import tqdm
from root_output import SingleRootFile

import ROOT

import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')


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

    def summary(self):
        # Implement the logic to summarize the results
        print(self)
        pass
    
    def write_to_file(self):
        if self.output_file is not None:
            with open(self.output_file, 'w') as file:
                # Implement the logic to write results to file
                file.write(str(self.results))
        self.summary()
            
class PrintEventAnalysis(BaseAnalysis):
    def process_event(self, event):
        # print("Event:", event)
        self.results.append(event)

class MultiplicityAnalysis(BaseAnalysis):              
    def process_event(self, event):
        multiplicity = event['multiplicity']
        # print("Multiplicity:", multiplicity)
        self.results.append(multiplicity)

class MultiplicityAnalysisRoot(BaseAnalysis):
    def user_init(self):
        self.root_output = SingleRootFile()
        self.tn_mult = ROOT.TNtuple("multiplicity", "multiplicity", "multiplicity:track_count")
        self.root_output.add(self.tn_mult)
        self.results.append(self.root_output)

    def process_event(self, event):
        multiplicity = event['multiplicity']
        track_count = len(event['track_data_pt'])
        self.tn_mult.Fill(multiplicity, track_count)
        # print("Multiplicity:", multiplicity)

def charged_particles_psjv(event, m=0.13957):
    rv = std.vector[fj.PseudoJet]()
    for i in range(len(event['track_data_pt'])):
      psj = fj.PseudoJet()
      psj.reset_PtYPhiM(event['track_data_pt'][i], event['track_data_eta'][i], event['track_data_phi'][i], m)
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
        bg_ymax = 1.0
        bg_grid_spacing = 0.1
        self.bg_estimator = fj.GridMedianBackgroundEstimator(bg_ymax, bg_grid_spacing)
        
        self.data_input_name = 'aliceRun3'
        self.n_accepted_jets = 0
        self.n_accepted_events = 0
        self.set_centrality(-1, 1001)
        
    def set_centrality(self, centmin, centmax):
        self.centrality_min = centmin
        self.centrality_max = centmax
        
    def process_event(self, event_struct):
        event = event_struct
        if isinstance(event_struct, GenericObject):
            event = event_struct.data[self.data_input_name]
        multiplicity = event['multiplicity']
        centrality = event['centrality']
        if centrality < self.centrality_min or centrality > self.centrality_max:
            return
        self.n_accepted_events += 1
        track_count = len(event['track_data_pt'])
        # estimate event background rho with grid estimator
        part_psjv = charged_particles_psjv(event)
        self.bg_estimator.set_particles(part_psjv)
        rho = self.bg_estimator.rho()
        # print('rho:', rho, 'centrality:', centrality, 'track_count:', track_count, 'part_psjv.size():', part_psjv.size())
        psjv = charged_particles_psjv(event)
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
            self.n_accepted_jets += 1
        self.tn_mult.Fill(multiplicity, track_count, centrality, jet_count)
        
def parse_yaml(file_path):
    with open(file_path, 'r') as file:
        tree_structure = yaml.safe_load(file)
    return tree_structure

class AnalysisEngine(GenericObject):
    def __init__(self, root_file_path, tree_name, branches, **kwargs):
        super(AnalysisEngine, self).__init__(**kwargs)
        self.root_file_path = root_file_path
        self.tree_name = tree_name
        self.branches = branches
        self.analyses = []

    def add_analysis(self, analysis):
        self.analyses.append(analysis)

    # Efficiently iterate over the tree
    def run(self, step_size=1000):
        # Open the ROOT file
        file = uproot.open(self.root_file_path)
        # Access the tree
        tree = file[self.tree_name]
        
        # Efficiently iterate over the tree with a progress bar
        total_entries = tree.num_entries
        if self.user_entries:
          total_entries = self.user_entries
        # Efficiently iterate over the tree with a progress bar
        with tqdm(total=total_entries, desc="Processing events") as pbar:
            # for data in tree.iterate(self.branches, library="np", step_size="100MB"):
            for data in tree.iterate(self.branches, library="np", step_size=step_size):
                # Iterate over the events in the chunk
                for i in range(len(next(iter(data.values())))):
                    event = {branch: data[branch][i] for branch in self.branches}
                    self.process_event(event)
                    pbar.update(1)
                if pbar.n >= total_entries:
                    break
                                
    def process_event(self, event):
        for analysis in self.analyses:
            analysis.process_event(event)

    def write_results(self):
        for analysis in self.analyses:
            print(f'Writing results for {analysis.output_file}')
            analysis.write_to_file()

if __name__ == '__main__':
  # Example usage
  # read from command line the yaml file path and the root file path
  # then run the analysis engine
  parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
  parser.add_argument("yaml_file_path", type=str, help="Path to the YAML file.")
  parser.add_argument("root_file_path", type=str, help="Path to the ROOT file.")
  parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.")
  parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
  args = parser.parse_args()
  
  yaml_file_path = args.yaml_file_path
  root_file_path = args.root_file_path

  yaml_file_structure = parse_yaml(yaml_file_path)

  print(yaml_file_structure)

  # Loop over the top-level items in yaml_file_structure - the tree names...
  for tname, tinfo in yaml_file_structure.items():
      print(tname)
      
      # Ensure tname is a string and tinfo contains the branches
      tree_name = tname
      branches = tinfo['branches']
      
      engine = AnalysisEngine(root_file_path, tree_name, branches, user_entries=args.entries)

      # Add analyses
      root_output = SingleRootFile(args.output)
      print_event_analysis = PrintEventAnalysis(output_file='print_event_results.txt')
      multiplicity_analysis = MultiplicityAnalysis(output_file = 'multiplicity_results.txt')
      multiplicity_analysis_root = MultiplicityAnalysisRoot()
      jana = JetAnalysisRoot()

      # engine.add_analysis(print_event_analysis)
      engine.add_analysis(multiplicity_analysis)
      engine.add_analysis(multiplicity_analysis_root)
      engine.add_analysis(jana)

      # Run the engine
      engine.run()

      # Write results to files
      engine.write_results()
      
      root_output.close()