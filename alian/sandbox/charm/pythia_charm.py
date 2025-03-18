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

def make_psj(p):
	psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
	psj.set_user_index(p.index())
	return psj

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='pythia_charm.root', type=str)
	parser.add_argument('--jet-pt-min', help='jet pt min', default=10.0, type=float)
	args = parser.parse_args()

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	Rs = [0.4]
	jet_defs = {}
	jet_selectors = {}
	for R in Rs:
		jet_defs[R] = fj.JetDefinition(fj.antikt_algorithm, R)
		jet_selectors[R] = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(1 - R * 1.05)

	# output
	fout = SingleRootFile(args.output)
	fout.root_file.cd()
	tn_events = ROOT.TNtuple("tn_events", "tn_events", "sigma:weight:code")
	tn = ROOT.TNtuple("tn", "tn", "sigma:weight:code:pt:eta:phi:mass:R:leadpid")
	tn_partons = ROOT.TNtuple("tn_partons", "tn_partons", "sigma:weight:code:pid:pt:eta:phi:mass")
	fout.add(tn)
 
	# configure pythia
	mycfg = ['HadronLevel:all = off', f'PhaseSpace:pThatMin = {args.jet_pt_min}']
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	# event loop
	if args.nev < 100:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		abs_fs_part_codes = [abs(p.id()) for p in pythia.event if p.isFinal()]
		if 4 not in abs_fs_part_codes:
			continue
		# parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		pythia_info = Pythia8.getInfo(pythia)
		sigma = pythia_info.sigmaGen()
		weight = pythia.info.weight()
		leading_process_code = pythia_info.code()
		tn_events.Fill(sigma, weight, leading_process_code)
		_ = [tn_partons.Fill(sigma, weight, leading_process_code, p.id(), p.pT(), p.eta(), p.phi(), p.m()) for p in pythia.event if p.isFinal() and p.isParton()]

		# jet finding
		parts = vector[fj.PseudoJet]([make_psj(p) for p in pythia.event if p.isFinal() and p.isParton()])
		for R in Rs:
			jet_def = jet_defs[R]
			jet_selector = jet_selectors[R]
			jets = jet_selector(jet_def(parts))
			for j in jets:
				leading = fj.sorted_by_pt(j.constituents())[0]
				leading_pythia_id = leading.user_index()
				leading_pythia = pythia.event[leading_pythia_id]
				tn.Fill(sigma, weight, leading_process_code, j.perp(), j.eta(), j.phi(), j.m(), R, leading_pythia.id())

	pythia.stat()

	print(type(pythia))

	fout.close()

if __name__ == '__main__':
	main()