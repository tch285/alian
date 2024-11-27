#!/usr/bin/env python

import argparse
import heppyy
from alian.io import data_io
from alian.utils import data_fj
import pandas as pd
import numpy as np
import ROOT
import math

def main():
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
	parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
	parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)

	args = parser.parse_args()

	# initialize an output root file
	tfile = ROOT.TFile(args.output, "recreate")
	tfile.cd()  # making sure we are in the right 'directory' to create some histograms
	h_track_pt 	= ROOT.TH1F("track_pt", "track pt; p_{T} (GeV/c); counts", 200, 0, 200)
	h_track_pt_cut = ROOT.TH1F("track_pt_cut", "track pt with cut; p_{T} (GeV/c); counts", 200, 0, 200)
	h_track_phi = ROOT.TH1F("track_phi", "track phi", 360, 0, 2.*math.pi)
	h_track_eta = ROOT.TH1F("track_eta", "track eta", 100, -1.5, 1.5)
	h_track_pt_cent = ROOT.TH2F("track_pt_cent", "centrality vs track pt;p_T (GeV/c);centrality", 200, 0, 200, 100, 0, 100)

	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# show some internal information of the data source
	print(data_source)

	# event loop using the data source directly
	for e in data_source.next_event():
		# get all the tracks into a vector of pseudojets - prep for jet finding
		# psjv = data_fj.data_tracks_to_pseudojets(e, lhc_run=args.lhc_run)
		# do what you will with these...
		# print(e.counter, e.multiplicity, e.centrality, len(psjv))

		# fill the histograms in a standard loop - not very efficient but simple
		for i in range(len(e.data['track_data_pt'])):
			h_track_pt.Fill(e.data['track_data_pt'][i])
			h_track_phi.Fill(e.data['track_data_phi'][i])
			h_track_eta.Fill(e.data['track_data_eta'][i])
			h_track_pt_cent.Fill(e.data['track_data_pt'][i], e.centrality)
			if e.centrality <= 10:
				h_track_pt_cut.Fill(e.data['track_data_pt'][i])

	# write the histograms to the output file
	# h_track_pt.Write()
	# ...
	# or more efficient way to write all histograms in the file
	tfile.Write()
	# close the file
	tfile.Close()


if __name__ == '__main__':
		main()
