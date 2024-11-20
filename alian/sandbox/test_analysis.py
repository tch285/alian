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

from data_io import Run3FileInput, Run2FileInput, AnalysisEngine
# from analysis import JetAnalysisRoot
from jet_analysis import JetAnalysisRoot

if __name__ == '__main__':
    # Example usage
    # read from command line the yaml file path and the root file path
    # then run the analysis engine
    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser.add_argument("yaml_file_path", type=str, help="Path to the YAML file.")
    parser.add_argument("input_file", type=str, help="Input file or file containing a list of input files.")
    parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1)
    parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
    parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)
    args = parser.parse_args()

    yaml_file_path = args.yaml_file_path
    input_file = args.input_file

    # Add analyses
    root_output = SingleRootFile(args.output)

    jana = None
    data_source = None
    jana = JetAnalysisRoot(lhc_run=args.lhc_run)
    if args.lhc_run == 3:
        data_source = Run3FileInput(input_file, yaml_file=yaml_file_path, n_events=args.entries)
    elif args.lhc_run == 2:
        data_source = Run2FileInput(input_file, yaml_file=yaml_file_path, n_events=args.entries)

    if not jana or not data_source:
        print("Invalid run number. Please specify either 2 or 3.")
        sys.exit(1)

    jana.set_centrality(0, 10)
    engine = AnalysisEngine()
    engine.add_source(data_source, is_driver=True, rename=jana.data_input_name)
    engine.add_analysis(jana)

    print(engine)
    # Run the engine
    engine.run()

    # Write results to files
    engine.write_results()
    root_output.close()
