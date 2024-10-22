#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import yasp

import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from heppyy.util.mputils import logbins
from heppyy.pythia_util import configuration as pyconf

import ROOT
import math
import array

from yasp import GenericObject
from alian.sandbox.root_output import SingleRootFile


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='pythia_jet_axis_output.root', type=str)
	parser.add_argument('--jet-pt-min', help='jet pt min', default=20.0, type=float)
	args = parser.parse_args()

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	Rs = [0.2, 0.4]
	jet_defs = {}
	jet_selectors = {}
	for R in Rs:
		jet_defs[R] = fj.JetDefinition(fj.antikt_algorithm, R)
		jet_selectors[R] = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(1 - R * 1.05)
  
	jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
	reclusterer_wta =  fj.contrib.Recluster(jet_def_wta)

	sd01 = fj.contrib.SoftDrop(0, 0.1, 1.0)
	sd02 = fj.contrib.SoftDrop(0, 0.2, 1.0)
 
	fout = SingleRootFile(args.output)
	fout.root_file.cd()
	tn = ROOT.TNtuple("tn", "tn", "pt:eta:phi:mass:wtastd:wtasd01:wtasd02:R")
	fout.add(tn)
 
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return
	if args.nev < 10:
		args.nev = 10
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		for R in Rs:
			jet_def = jet_defs[R]
			jet_selector = jet_selectors[R]
			jets = jet_selector(jet_def(parts))
			for j in jets:
				jet_wta = reclusterer_wta.result(j)
				jet_sd01 = sd01.result(j)
				jet_sd02 = sd02.result(j)
				wtastd = jet_wta.delta_R(j)
				wtasd01 = jet_sd01.delta_R(j)
				wtasd02 = jet_sd02.delta_R(j)
				tn.Fill(j.perp(), j.eta(), j.phi(), j.m(), wtastd, wtasd01, wtasd02, R)

	pythia.stat()

	print(type(pythia))

	fout.close()

if __name__ == '__main__':
	main()