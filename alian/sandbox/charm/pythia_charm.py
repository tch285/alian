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
	parser.add_argument('--charm-pt-min', help='charm pt min', default=5.0, type=float)
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
	tn_events = ROOT.TNtuple("tn_events", "tn_events", "sigma:weight:code:out1pid:out2pid")
	tn_jets = ROOT.TNtuple("tn_jets", "tn_jets", "sigma:weight:code:pt:eta:phi:mass:R:leadpid:has_charm:leadpt")
	tn_partons = ROOT.TNtuple("tn_partons", "tn_partons", "sigma:weight:code:pid:pt:eta:phi:mass")
 
	# configure pythia
	mycfg = ['HadronLevel:all = off', 
          f'PhaseSpace:pThatMin = {args.jet_pt_min}',
          'PhaseSpace:bias2Selection = off']
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	# event loop
	if args.nev < 100:
		args.nev = 100
	accepted = 0
	pbar = tqdm.tqdm(total=args.nev)
	while accepted < args.nev:
		if not pythia.next():
			continue
		abs_fs_part_codes = [abs(p.id()) for p in pythia.event if p.isFinal()]
		if 4 not in abs_fs_part_codes:
			continue
		# charm above a pt
		charm_above_pt = [p for p in pythia.event if p.isFinal() and abs(p.id()) == 4 and p.pT() > args.jet_pt_min and abs(p.eta()) < 1]
		if len(charm_above_pt) == 0:
			continue
		# parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		pythia_info = Pythia8.getInfo(pythia)
		sigma = pythia_info.sigmaGen()
		weight = pythia.info.weight()
		leading_process_code = pythia_info.code()
		tn_events.Fill(sigma, weight, leading_process_code, pythia.event[4].id(), pythia.event[5].id())
		_ = [tn_partons.Fill(sigma, weight, leading_process_code, p.id(), p.pT(), p.eta(), p.phi(), p.m()) for p in pythia.event if p.isFinal() and p.isParton()]

		# jet finding
		parts = vector[fj.PseudoJet]([make_psj(p) for p in pythia.event if p.isFinal() and p.isParton()])
		for R in Rs:
			jet_def = jet_defs[R]
			jet_selector = jet_selectors[R]
			jets = jet_selector(jet_def(parts))
			if R==0.4 and len(jets) > 0:
				accepted += 1
				pbar.update(1)
			for j in jets:
				leading = fj.sorted_by_pt(j.constituents())[0]
				leading_pythia_id = leading.user_index()
				leading_pythia = pythia.event[leading_pythia_id]
				has_charm = any([abs(pythia.event[p.user_index()].id()) == 4 for p in j.constituents()])
				tn_jets.Fill(sigma, weight, leading_process_code, j.perp(), j.eta(), j.phi(), j.m(), R, leading_pythia.id(), int(has_charm), leading.perp())

	pythia.stat()

	print(type(pythia))

	fout.close()

if __name__ == '__main__':
	main()
 
 