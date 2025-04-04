#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import yasp
import tqdm
import cppyy
import argparse
import itertools
import ROOT

import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from heppyy.pythia_util import configuration as pyconf

from yasp import GenericObject
from alian.sandbox.root_output import SingleRootFile

def count_lhe_events(lhe_filename):
		"""Count the number of events in an LHE file."""
		count = 0
		with open(lhe_filename, 'r') as f:
				for line in f:
						if '<event>' in line:
								count += 1
		return count
	
def make_psj(p):
	psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
	psj.set_user_index(p.index())
	return psj

def find_Z0_daughter_pairs(parts, event):
	"""find Z0 from muons"""
	muons = []
	for p in parts:
		pdg_code = event[p.user_index()].id()
		if abs(pdg_code) == 13:
			pythia_part = event[p.user_index()]
			mother1 = pythia_part.mother1()
			mother2 = pythia_part.mother2()
			if mother1 > 0 and event[mother1].id() == 23:
				muons.append(p)
			elif mother2 > 0 and event[mother2].id() == 23:
				muons.append(p)
	mu_pairs = itertools.combinations(muons, 2)
	Z0daughters = []
	Z0pairs = []
	Z0s = []
	for pair in mu_pairs:
		# check if they are from the same Z0
		p1_mother1 = event[pair[0].user_index()].mother1()
		p2_mother1 = event[pair[1].user_index()].mother1()
		p1_mother2 = event[pair[0].user_index()].mother2()
		p2_mother2 = event[pair[1].user_index()].mother2()
		if p1_mother1 == p2_mother1 and p1_mother2 == p2_mother2:
			Z0pairs.append(pair)
			Z0daughters.append(pair[0].user_index())
			Z0daughters.append(pair[1].user_index())
			Z0s.append(make_psj(event[p1_mother1]))
			break
	return Z0daughters, Z0pairs, Z0s

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet taking lhe file [from powheg for example]', prog=os.path.basename(__file__))
	parser.add_argument('input', help="lhe file [from powheg would be likely pwgevents.lhe]", type=str)
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('-o','--output', help='root output filename', default='pythia_zjet_from_lhe.root', type=str)

	args = parser.parse_args()

	nevents_lhe = count_lhe_events(args.input)
	if nevents_lhe < 0:
		print(f"[e] error counting events in {args.input}")
		return
	if args.nev > nevents_lhe:
		args.nev = nevents_lhe

	pythia = Pythia8.Pythia()

	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(5) * fj.SelectorAbsEtaMax(20)
	print(jet_def)

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	lund_gen = fj.contrib.LundGenerator(jet_def_lund)
	print('making lund diagram for all jets...')
	print(f' {lund_gen.description()}')

	# output
	fout = SingleRootFile(args.output)
	fout.root_file.cd()
	tn_events = ROOT.TNtuple("tn_events", "tn_events", "iev:sigma:weight:code:out1pid:out2pid")
	tn_Z0 = ROOT.TNtuple("tn_Z0", "tn_Z0", "iev:pt:eta:phi:mass:pid:status")
	tn_Zjet = ROOT.TNtuple("tn_Zjet", "tn_Zjet", "iev:pt:eta:phi:mass:Zpt:Zeta:Zphi:Zmass")

	mycfg = [	'HardQCD:all=off', f'Beams:LHEF = {args.input}', 'Beams:frameType = 4', 
						'Next:numberCount = 0', 'Next:numberShowEvent = 0', 'Next:numberShowInfo = 0', 'Next:numberShowProcess = 0', 'Stat:showProcessLevel = on']
	# pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	pythia = Pythia8.Pythia()
	for s in mycfg:
		pythia.readString(s)
	if not pythia.init():
		print("[e] pythia initialization failed.")
		return
	if not pythia:
		print("[e] pythia initialization failed.")
		return
	if args.nev < 10:
		args.nev = 10
	# event loop
	if args.nev < 100:
		args.nev = nevents_lhe
	accepted = 0
	pbar = tqdm.tqdm(total=args.nev)
	total_events = 0
	iev = 0
	while accepted < args.nev:
		total_events += 1
		iev += 1
		if not pythia.next():
			if total_events < nevents_lhe:
				print(f"[w] event {total_events} failed but continuing")
				continue
			print(f"[i] processed {total_events} events, accepted {accepted}")
			print("[i] end of LHE file reached")
			break
		# event processing...
		pythia_info = Pythia8.getInfo(pythia)
		sigma = pythia_info.sigmaGen()
		# weight = pythia.info.weight()
		weight = 1
		leading_process_code = pythia_info.code()
		tn_events.Fill(iev, sigma, weight, leading_process_code, pythia.event[4].id(), pythia.event[5].id())
		parts = vector[fj.PseudoJet]([make_psj(p) for p in pythia.event if p.isFinal()])

		all_parts = vector[fj.PseudoJet]([make_psj(p) for p in pythia.event])
		Z0s = [p for p in all_parts if abs(pythia.event[p.user_index()].id()) == 23]
		# print(f'Z0: {Z0}')
		if len(Z0s) > 0:
			accepted += 1
			pbar.update(1)
		else:
			Z0s = None
			continue
		for z in Z0s:
			tn_Z0.Fill(iev, z.perp(), z.eta(), z.phi(), z.m(), pythia.event[z.user_index()].id(), pythia.event[z.user_index()].status())

		# find Z0 from muons
		Z0daughters, Z0pairs, Z0s = find_Z0_daughter_pairs(parts, pythia.event)
		if len(Z0pairs) > 0:
			jf_parts = vector[fj.PseudoJet]([p for p in parts if p.user_index() not in Z0daughters])
			jets = jet_selector(jet_def(jf_parts))
			for j in jets:
				for Z in Z0s:
					tn_Zjet.Fill(iev, j.perp(), j.eta(), j.phi(), j.m(), Z.perp(), Z.eta(), Z.phi(), Z.m())

	pythia.stat()
	fout.close()

if __name__ == '__main__':
	main()
