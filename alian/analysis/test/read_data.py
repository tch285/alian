#!/usr/bin/env python

import argparse
import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
from alian.io.root_io import SingleRootFile
from alian.io import data_io
from alian.utils import data_fj
from alian.steer.glob import globals
from alian.utils.treewriter import RTreeWriter
from alian.analysis.base import BaseAnalysis
from alian.analysis.base import CEventSubtractor

import ROOT

import numpy as np
import array


def main():
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
	parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)
	parser.add_argument('--cent-min', type=int, help='Minimum centrality', default=-1)
	parser.add_argument('--cent-max', type=int, help='Maximum centrality', default=101)
	parser.add_argument('--no-tqdm', action='store_true', help='Disable tqdm progress bars', default=False)
	parser.add_argument('--part-eta-max', default=0.9, type=float)

	args = parser.parse_args()
	print(args)

	globals.tqdm_silent = args.no_tqdm
 
	if args.tree_struct is None:
		args.tree_struct = data_io.get_default_tree_structure(args.lhc_run)
	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# get the fj banner out of the way
	fj.ClusterSequence().print_banner()
	parts_selector = fj.SelectorAbsEtaMax(args.part_eta_max)

	# event loop using the data source directly
	for i,e in enumerate(data_source.next_event()):
		if e.centrality < args.cent_min or e.centrality > args.cent_max:
			continue
		psjv = data_fj.data_tracks_to_pseudojets(e, lhc_run=args.lhc_run)
		e.psjv = parts_selector(psjv)
		print(f'Event {i}, centrality = {e.centrality}, multiplicity = {e.multiplicity}, n_tracks = {len(e.psjv)}, z_vtx = {e.z_vtx}')
		for p in e.psjv:
			print(p, f'p.perp() = {p.perp()}, p.eta() = {p.eta()}, p.phi() = {p.phi()}')
		if i > args.entries and args.entries > 0:
			break

  # finalize the analyses - write the output, close the file, etc.
	# an.finalize()
	# if cs:
	# 	an_cs.analyze(e)
	root_file.close()

if __name__ == '__main__':
		main()
