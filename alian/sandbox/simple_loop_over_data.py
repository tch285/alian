#!/usr/bin/env python

import argparse
import heppyy
from alian.io import data_io
from alian.utils import data_fj


def main():
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
	parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
	parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)

	args = parser.parse_args()

	# initialize the data input
	# data_source = None
	# if args.tree_struct == None:
	# 	args.tree_struct = data_io.get_default_tree_structure(lhc_run=args.lhc_run)
	# if args.lhc_run == 3:
	# 	data_source = data_io.Run3FileInput(args.input_file, yaml_file=args.tree_struct, n_events=args.entries)
	# elif args.lhc_run == 2:
	# 	data_source = data_io.Run2FileInput(args.input_file, yaml_file=args.tree_struct, n_events=args.entries)
  
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	print(data_source)

	# event loop using the data source directly
	for e in data_source.next_event():
		# get all the tracks into a vector of pseudojets - prep for jet finding
		psjv = data_fj.data_tracks_to_pseudojets(e, lhc_run=args.lhc_run)
		# do what you will with these...
		# print(e.counter, e.multiplicity, e.centrality, len(psjv))
 
if __name__ == '__main__':
		main()