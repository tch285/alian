#!/usr/bin/env python

import sys 

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

from data_io import Run3FileInput, AnalysisEngine
# from analysis import JetAnalysisRoot
from jet_analysis import JetAnalysisRoot, JetAnalysisRootRun2

if __name__ == '__main__':
    # Example usage
    # read from command line the yaml file path and the root file path
    # then run the analysis engine
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser.add_argument("yaml_file_path", type=str, help="Path to the YAML file.")
    parser.add_argument("root_file_path", type=str, help="Path to the ROOT file.")
    parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1)
    parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
    args = parser.parse_args()

    yaml_file_path = args.yaml_file_path
    root_file_path = args.root_file_path

    jana2 = JetAnalysisRootRun2()
    sys.exit(0)
    
    # Add analyses
    root_output = SingleRootFile(args.output)
    jana = JetAnalysisRoot()
    jana.set_centrality(0, 10)
    
    run3data = Run3FileInput(root_file_path, yaml_file=yaml_file_path, n_events=args.entries)

    engine = AnalysisEngine()
    engine.add_source(run3data, is_driver=True, rename=jana.data_input_name)
    engine.add_analysis(jana)

    print(engine)
    # Run the engine
    engine.run()

    # Write results to files
    engine.write_results()
    root_output.close()
