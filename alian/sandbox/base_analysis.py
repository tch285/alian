#!/usr/bin/env python

import argparse
import yaml
import uproot

from yasp import GenericObject
from tqdm import tqdm
from root_output import SingleRootFile

import ROOT

import heppyy

class BaseAnalysis(GenericObject):
    _defaults = {}
    
    def __init__(self, **kwargs):
        super(BaseAnalysis, self).__init__(**kwargs)
        self.results = []
        for k, val in self.__class__._defaults.items():
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, val)
                # print(f'[i] Setting {k} to default value {val}', self.__class__, self.__getattr__(k))
        self.init()
        print(self)

    def init(self):
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
    def run(self):
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
            for data in tree.iterate(self.branches, library="np", step_size=1):
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

