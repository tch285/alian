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
            event.source = {}
            event.source[self.driver_source.name] = self.driver_source
            for source in self.sources.values():
                if source.is_driver:
                    continue
                # print(f'AnalysisEngine: Getting next event from {source.name}')
                if source.auto_next:
                  event.data[source.name] = next(source.next_event())
                else:
                  event.data[source.name] = source.next_event()
                event.source[source.name] = source
            for analysis in self.analyses:
                analysis.process_event(event)

    def write_results(self):
        for analysis in self.analyses:
            print(f'AnalysisEngine: Writing results for {analysis.output_file}')
            analysis.write_to_file()
