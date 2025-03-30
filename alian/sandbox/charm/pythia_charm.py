#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import json
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


def idx_match_to_out_parton(jet, out_parton1, out_parton2, jet_R0=0.4):
    # print(jet, out_parton1, out_parton2)
    # print(jet.delta_R(out_parton1), jet.delta_R(out_parton2))
    # if jet.delta_R(out_parton1) < jet_R0:
    # 		return out_parton1.user_index()
    # if jet.delta_R(out_parton2) < jet_R0:
    # 		return out_parton2.user_index()
	if jet.delta_R(out_parton1) < jet.delta_R(out_parton2):
		if jet.delta_R(out_parton1) < jet_R0:
			return out_parton1.user_index(), jet.delta_R(out_parton1)
		else:
			if jet.delta_R(out_parton2) < jet_R0:
				return out_parton2.user_index(), jet.delta_R(out_parton2)
	return -1, -1


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
		jet_selectors[R] = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_min + 10.) * fj.SelectorAbsEtaMax(1 - R * 1.05)

	# output
	fout = SingleRootFile(args.output)
	fout.root_file.cd()
	tn_jets = ROOT.TNtuple("tn_jets", "tn_jets", "nev:xsev:ev_weight:code:pt:eta:phi:mass:R:leadpid:has_charm:leadpt:outpid:dRmatch")
	tn_partons = ROOT.TNtuple("tn_partons", "tn_partons", "nev:xsev:ev_weight:code:pid:pt:eta:phi:mass")
	tn_norm   = ROOT.TNtuple('tn_norm', 'tn_norm', 'nev:xsec:xsec_err:sum_of_weights')
	tn_hard 	= ROOT.TNtuple('tn_hard', 'tn_hard', 'nev:xsec:ev_weight:x1:x2:QFac:id1:id2:id3:pt3:eta3:id4:pt4:eta4')
	tn_events = ROOT.TNtuple("tn_events", "tn_events", "nev:xsev:ev_weight:code:out1pid:out2pid:nparts:x1:x2:QFac")
	# tn_events = ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts:x1:x2:QFac:ncoll')

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
	iev = 0
	sum_weights = 0.0
	pbar = tqdm.tqdm(total=args.nev)
	while accepted < args.nev:
		if not pythia.next():
			continue
		iev += 1
		pythia_info = Pythia8.getInfo(pythia)
		sigma = pythia_info.sigmaGen()
		ev_weight = pythia_info.weight()
		leading_process_code = pythia_info.code()
		sum_weights += ev_weight
		abs_fs_part_codes = [abs(p.id()) for p in pythia.event if p.isFinal()]
		if 4 not in abs_fs_part_codes:
			continue
		# charm above a pt
		charm_above_pt = [p for p in pythia.event if p.isFinal() and abs(p.id()) == 4 and p.pT() > args.charm_pt_min and abs(p.eta()) < 1]
		if len(charm_above_pt) == 0:
			continue
		# parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		pythia_info = Pythia8.getInfo(pythia)
		sigma = pythia_info.sigmaGen()
		ev_weight = pythia_info.weight()
		leading_process_code = pythia_info.code()
		# tn_events.Fill(iev, sigmaGen, ev_weight, pythia.event.size(), _info.x1(), _info.x2(), _info.QFac(), ncoll)
		# tn_events = ROOT.TNtuple("tn_events", "tn_events", "nev:xsev:ev_weight:code:out1pid:out2pid:nparts:x1:x2:QFac")
		tn_events.Fill(iev, sigma, ev_weight, leading_process_code, pythia.event[4].id(), pythia.event[5].id(), len(pythia.event), pythia_info.x1(), pythia_info.x2(), pythia_info.QFac())
		_ = [tn_partons.Fill(iev, sigma, ev_weight, leading_process_code, p.id(), p.pT(), p.eta(), p.phi(), p.m()) for p in pythia.event if p.isFinal() and p.isParton()]

		# jet finding
		parts = vector[fj.PseudoJet]([make_psj(p) for p in pythia.event if p.isFinal() and p.isParton()])

		in_parton1 = pythia.event[3]
		in_parton2 = pythia.event[4]
		out_parton1 = pythia.event[5]
		out_parton2 = pythia.event[6]
		fj_out_partons = vector[fj.PseudoJet]()
		for p in [out_parton1, out_parton2]:
			fjp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
			fjp.set_user_index(p.index())
			fj_out_partons.push_back(fjp)
		# tn_hard = ROOT.TNtuple('tn_hard', 'tn_hard', 'nev:xsec:ev_weight:x1:x2:QFac:id1:id2:id3:id4:pt3:eta3:pt4:eta4')
		tn_hard.Fill(iev, sigma, ev_weight, pythia_info.x1(), pythia_info.x2(), pythia_info.QFac(), in_parton1.id(), in_parton2.id(), out_parton1.id(), out_parton1.pT(), out_parton1.eta(), out_parton2.id(), out_parton2.pT(), out_parton2.eta())

		for R in Rs:
			jet_def = jet_defs[R]
			jet_selector = jet_selectors[R]
			jets = jet_selector(jet_def(parts))
			# if R==0.4 and len(jets) > 0:
			# 	accepted += 1
			# 	pbar.update(1)
			for j in jets:
				leading = fj.sorted_by_pt(j.constituents())[0]
				leading_pythia_id = leading.user_index()
				leading_pythia = pythia.event[leading_pythia_id]
				has_charm = any([abs(pythia.event[p.user_index()].id()) == 4 for p in j.constituents()])
				_idx_match_to_out_parton, dRmatch = idx_match_to_out_parton(j, fj_out_partons[0], fj_out_partons[1], R)
				out_parton_id = 0
				if _idx_match_to_out_parton != -1:
					out_parton = pythia.event[_idx_match_to_out_parton]
					out_parton_id = out_parton.id()
				if out_parton_id != 0:
					tn_jets.Fill(iev, sigma, ev_weight, leading_process_code, j.perp(), j.eta(), j.phi(), j.m(), R, leading_pythia.id(), int(has_charm), leading.perp(), out_parton_id, dRmatch)
					accepted += 1
					pbar.update(1)

	pythia.stat()

	# Get the total cross section and weight sum
	sigma_gen = pythia_info.sigmaGen()
	sigma_gen_err = pythia_info.sigmaErr()
	weight_sum = pythia_info.weightSum()  # Same as sum_weights
	nAccepted = pythia_info.nAccepted()

	# tn_norm =   ROOT.TNtuple('tn_norm', 'tn_norm', 'nev:xsec:xsec_err:sum_of_weights')
	tn_norm.Fill(nAccepted, sigma_gen, sigma_gen_err, weight_sum)

	# Save to JSON summary
	json_file = args.output.replace('.root', '.json')
	with open(json_file, "w") as f:
			json.dump({
					"n_accepted": nAccepted,
					"sigma_gen": sigma_gen,
					"sum_weights": weight_sum
			}, f, indent=2)

	# Output for verification
	print(f"sigma_gen = {sigma_gen:.3f} +- {sigma_gen_err:.3f} mb")
	print(f"Sum of weights = {weight_sum:.3f} [check: {sum_weights:.3f}]")
	print(f"Number of accepted events = {nAccepted}")

	fout.close()

if __name__ == '__main__':
	main()
