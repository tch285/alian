#!/usr/bin/env python

import uproot
import yaml
import pandas as pd
from yasp import GenericObject
from tqdm import tqdm

class Run3FileInput(GenericObject):
    def __init__(self, file_list, yaml_file=None, tree_name=None, branches=None, **kwargs):
        super(Run3FileInput, self).__init__(**kwargs)
        self.file_list = file_list
        if isinstance(self.file_list, str):
            if self.file_list.endswith(".txt"):
                with open(self.file_list, "r") as file:
                    self.file_list = file.readlines()
                self.file_list = [x.strip() for x in self.file_list]
            else:
                self.file_list = [self.file_list]
        self.tree_name = tree_name
        self.branches = branches
        if yaml_file is not None:
            with open(yaml_file, 'r') as file:
                tree_structure = yaml.safe_load(file)
                for tname, tinfo in tree_structure.items():
                    self.tree_name = tname
                    self.branches = tinfo['branches']
        if self.name is None:
            self.name = "FileIO"
        self.event = None
        self.event_count = 0
        self.n_events = kwargs.get('n_events', -1)

    # Efficiently iterate over the tree as a generator
    def next_event(self, step_size=1000):
        self.event = None
        pbar_total = None
        if self.n_events > 0:
            pbar_total = tqdm(total=self.n_events, desc="Total events")
        for root_file_path in tqdm(self.file_list, desc="Files"):
            # Open the ROOT file
            file = uproot.open(root_file_path)
            # Access the tree
            tree = file[self.tree_name]            
            # Efficiently iterate over the tree with a progress bar
            total_entries = tree.num_entries
            # Efficiently iterate over the tree with a progress bar
            with tqdm(total=total_entries, desc=f'...{root_file_path[-25:]}') as pbar:
                for data in tree.iterate(self.branches, library="np", step_size=step_size):
                    # Iterate over the events in the chunk
                    for i in range(len(next(iter(data.values())))):
                        self.event = {branch: data[branch][i] for branch in self.branches}
                        pbar.update(1)
                        self.event_count += 1
                        if pbar_total is not None:
                            pbar_total.update(1)
                            if pbar_total.n >= self.n_events:
                                break
                        yield self.event
                    if pbar.n >= total_entries:
                        break
                    if pbar_total is not None:
                        if pbar_total.n >= self.n_events:
                            break
            if pbar_total is not None:
                if pbar_total.n >= self.n_events:
                    pbar_total.close()
                    break
class Run2FileInput(GenericObject):
    def __init__(self, file_list, yaml_file, **kwargs):
        super(Run2FileInput, self).__init__(**kwargs)
        self.file_list = file_list
        if isinstance(self.file_list, str):
            if self.file_list.endswith(".txt"):
                with open(self.file_list, "r") as file:
                    self.file_list = file.readlines()
                self.file_list = [x.strip() for x in self.file_list]
            else:
                self.file_list = [self.file_list]
        self.setup_default_trees_and_branches(yaml_file)
        if self.name is None:
            self.name = "FileIORun2"
        self.event = None
        self.event_count = 0

    def setup_default_trees_and_branches(self, yaml_file_path):
        with open(yaml_file_path, 'r') as file:
            config = yaml.safe_load(file)
        
        self.tree_particle_name = 'PWGHF_TreeCreator/tree_Particle'
        self.tree_event_char_name = 'PWGHF_TreeCreator/tree_event_char'
        
        self.branches_particle = list(config['PWGHF_TreeCreator']['tree_Particle'].keys())
        self.branches_event_char = list(config['PWGHF_TreeCreator']['tree_event_char'].keys())
        
        self.branches = self.branches_particle + self.branches_event_char

    # Efficiently iterate over the tree as a generator
    def next_event(self):
        self.event = None
        pbar_total = None
        if self.n_events > 0:
            pbar_total = tqdm(total=self.n_events, desc="Total events")
        for root_file_path in tqdm(self.file_list, desc="Files"):
            # Open the ROOT file
            file = uproot.open(root_file_path)
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
            with tqdm(total=total_entries, desc=f'...{root_file_path[-25:]}') as pbar:
                for name, group in grouped_df:
                    self.event = {branch: group[branch].tolist() for branch in self.branches if branch in group}
                    yield self.event
                    pbar.update(1)
                    self.event_count += 1
                    if pbar.n >= total_entries:
                        break
                    if pbar_total is not None:
                        pbar_total.update(1)
                        if pbar_total.n >= self.n_events:
                            break
            if pbar_total is not None:
                if pbar_total.n >= self.n_events:
                    pbar_total.close()
                    break

            
class AnalysisEngine(GenericObject):
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(AnalysisEngine, cls).__new__(cls)
        return cls._instance

    def __init__(self, **kwargs):
        # if not hasattr(self, 'initialized'):  # Ensure __init__ is only called once
        if self.initialized:
            return
        super(AnalysisEngine, self).__init__(**kwargs)
        if self.sources is None:
            self.sources = {}
        if self.analyses is None:
            self.analyses = []
        self.initialized = True

    def add_source(self, source, is_driver=False, rename=None):
        if rename is not None:
            source.name = rename
        if is_driver:
            source.is_driver = True
        self.sources[source.name] = source
        count_drivers = sum([1 for source in self.sources.values() if source.is_driver])
        if count_drivers > 1:
            raise ValueError("AnalysisEngine: Only one driver can be specified")
        driver_source_list = [source for source in self.sources.values() if source.is_driver]
        if len(driver_source_list) > 0:
            self.driver_source = driver_source_list[0]
            
    def add_analysis(self, analysis):
        self.analyses.append(analysis)
        
    def run(self):
        if self.driver_source is None:
            raise ValueError("AnalysisEngine: No driver source specified")
        for e in self.driver_source.next_event():
            event = GenericObject()
            event.data = {}
            event.data[self.driver_source.name] = e
            for source in self.sources.values():
                if source.is_driver:
                    continue
                event.data[source.name] = source.next_event()
            for analysis in self.analyses:
                analysis.process_event(event)

    def write_results(self):
        for analysis in self.analyses:
            print(f'AnalysisEngine: Writing results for {analysis.output_file}')
            analysis.write_to_file()
