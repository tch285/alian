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

def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr

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

def deep_mother(idx, pythia):
	# find the mother of the particle with index idx
	mother = pythia.event[idx].mother1()
	if mother == -1:
		print('[?] mother is -1 for', idx)
		return 0
	if mother in [5, 6]:
		return pythia.event[mother].id()
	mother_pid = pythia.event[mother].id()
	if abs(mother_pid) == 4:
		return deep_mother(mother, pythia)
	else:
		return mother_pid

def follow_two_mothers(mother1, mother2, pythia, indent=0, debug=False):
	sindent = ' ' * indent
	m1x = deep_parton_mother(mother1, pythia, indent + 1, debug=debug)
	m2x = deep_parton_mother(mother2, pythia, indent + 1, debug=debug)
	return [m1x, m2x]
 
def deep_parton_mother(idx, pythia, indent=0, debug=False):
	# find the mother of the particle with index idx
	sindent = ' ' * indent
	mother1 = pythia.event[idx].mother1()
	mother2 = pythia.event[idx].mother2()
	sidx = f'[{idx} - {pythia.event[idx].id()}]'
	if mother2 > 0:
		if mother1 != mother2 and mother2 != 0:
			if debug: print(f'{sindent} [?] {sidx} different parton mothers', mother1, mother2)
			return follow_two_mothers(mother1, mother2, pythia, indent, debug=debug)
	if type(mother1) == list:
		if debug: print(f'{sindent} [!] {sidx} mother is not int', mother1)
		return mother1
	while mother1 > 0:
		if debug: print(f'{sindent} - {sidx} mother1: {mother1}, mother2: {mother2}')
		if mother1 in [5, 6]:
			return mother1
		mother1 = deep_parton_mother(mother1, pythia, indent + 1)
		if type(mother1) == list:
			len_mother1 = len(mother1)
			if len_mother1 == 1:
				mother1 = mother1[0]
			else:
				# get the lowest number non zero but zero if only zero left
				mother1 = [m for m in mother1 if m != 0]
				if len(mother1) == 0:
					mother1 = 0
				else:
					mother1 = min(mother1)
	if mother1 < 4:
		mother1 = 0
	return mother1

def deep_parton_mother_pid(idx, pythia, debug=False):
	mother_pid = 0
	mothers = deep_parton_mother(idx, pythia, debug=debug)
	if type(mothers) == list:
		if len(mothers) == 1:
			mother = mothers[0]
			if mother == 0:
				if debug: print(f'[!] {idx} mother is None for')
				return mother_pid
		if len(mothers) == 2:
			if debug: print(f'[!] {idx} complex mother is {mothers}')
			mother_pid = -1
			return mother_pid
		if len(mothers) == 0:
			if debug: print(f'[!] {idx} mother is None for')
			return mother_pid
	else:
		mother = mothers
	mother_pid = pythia.event[mother].id()
	if debug: print(f'[=>] {idx} {pythia.event[idx].pT()} <- mother is {mother} with pid {mother_pid}')
	return mother_pid


def deep_mother_hadron(idx, pythia):
	# find the mother of the particle with index idx
	mother = pythia.event[idx].mother1()
	if mother == -1:
		print('[?] mother is -1 for', idx)
		return 0
	mother2 = pythia.event[idx].mother2()
	if mother2 == -1:
		mother2_pid = 0
		return 0
	if mother in [5, 6]:
		return pythia.event[mother].id()
	mother_pid = pythia.event[mother].id()
	if pythia.event[mother].isHadron():
		return deep_mother(mother, pythia)
	if abs(mother_pid) == 4:
		return deep_mother(mother, pythia)
	else:
 		return mother_pid

def deep_hadron_parton_mother(idx, pythia):
  # find the parton mother of the hadron with index idx
	mother1 = pythia.event[idx].mother1()
	mother2 = pythia.event[idx].mother2()
	
def trace_production_tree(idx, pythia, level=0):
    """
    Recursively print the production tree of the particle at index `idx`
    down to the parton level.
    """
    if idx <= 0:
      # print("Invalid index:", idx)
      return
    particle = pythia.event[idx]
    indent = "  " * level
    print(f"{indent}Idx: {idx}, pid: {particle.id()}, pT: {particle.pT():.2f}")
    
    # Get mothers (mother1 and mother2). If both are -1, we've reached the initial particle.
    mother1 = particle.mother1()
    mother2 = particle.mother2()
    if mother1 == -1 and mother2 == -1:
        print(f"{indent}No mother (primary particle)")
        return
    
    # Trace the first mother
    if mother1 > 0 or mother1 != 90:
        print(f"{indent}Mother1:")
        trace_production_tree(mother1, pythia, level + 1)
    
    # If a second distinct mother exists, trace it as well
    if mother2 > 0 and mother2 != mother1 and mother2 != 90:
        print(f"{indent}Mother2:")
        trace_production_tree(mother2, pythia, level + 1)
    
def print_ascii_production_tree(idx, pythia, prefix=""):
    """
    Recursively prints the production tree of the particle at index `idx`
    as an ASCII tree.
    """
    particle = pythia.event[idx]
    # Print current node information
    print(f"{prefix}{idx}: pid {particle.id()}, pT {particle.pT():.2f}")
    
    # Retrieve mothers and filter out invalid indices (-1)
    mothers = []
    mother1 = particle.mother1()
    mother2 = particle.mother2()
    if mother1 != -1:
        mothers.append(mother1)
    if mother2 != -1 and mother2 != mother1:
        mothers.append(mother2)
    
    # Determine the drawing prefix for each child
    for i, m in enumerate(mothers):
        if i == len(mothers) - 1:
            new_prefix = prefix + "    "
            branch = "└── "
        else:
            new_prefix = prefix + "│   "
            branch = "├── "
        print_ascii_production_tree(m, pythia, prefix + branch)

# Example usage:
# After generating an event and choosing a particle index (e.g., charm hadron),
# call this function:
#   print_ascii_production_tree(chosen_idx, pythia)
# This will print the ASCII tree to the terminal.


def top_down_production_tree(idx, pythia, prefix="", selected_idx=None, highlight=False):
    """
    Recursively prints the production tree of the particle at index `idx`
    as an ASCII tree, highlighting in red the branch starting from the selected index.
    
    Parameters:
      idx         : current particle index
      pythia      : Pythia event object
      prefix      : string prefix for drawing the tree structure
      selected_idx: if provided, the branch starting at this particle index is highlighted
      highlight   : internal flag to propagate the highlight to descendants
    """
    particle = pythia.event[idx]
    node_str = f"{idx}: pid {particle.id()}, pT {particle.pT():.2f}"
    
    # Check if this node should be printed in red
    if highlight or (selected_idx is not None and idx == selected_idx):
        node_str = f"\033[91m{node_str}\033[0m"  # ANSI red color
        new_highlight = True  # propagate highlight to children
    else:
        new_highlight = False

    print(f"{prefix}{node_str}")

    # Retrieve mothers and filter out invalid indices (-1)
    mothers = []
    mother1 = particle.mother1()
    mother2 = particle.mother2()
    if mother1 > 0:
        mothers.append(mother1)
    if mother2 > 0 and mother2 != mother1:
        mothers.append(mother2)

    # Recursively print each branch with proper tree branch symbols
    for i, m in enumerate(mothers):
        branch = "├── " if i < len(mothers) - 1 else "└── "
        top_down_production_tree(m, pythia, prefix + branch, selected_idx, new_highlight)
        
# Example usage:
# After generating an event and choosing a particle index to highlight (e.g., charm hadron index 42),
# call:
#   print_ascii_production_tree(chosen_idx, pythia, prefix="", selected_idx=42)
# The branch starting at particle 42 (and all its descendants) will be printed in red.

def print_ascii_daughter_tree(idx, pythia, prefix="", selected_idx=None, highlight=False):
    """
    Recursively prints the daughter tree of the particle at index `idx`
    as an ASCII tree. If `selected_idx` is provided, the branch starting at
    that particle is printed in red.

    Parameters:
      idx         : current particle index
      pythia      : Pythia event object
      prefix      : prefix string for formatting the tree
      selected_idx: if provided, the branch from this index onward is highlighted in red
      highlight   : propagate the highlight to descendants
    """
    particle = pythia.event[idx]
    node_str = f"{idx}: pid {particle.id()}, pT {particle.pT():.2f}"

    if highlight or (selected_idx is not None and idx == selected_idx):
        node_str = f"\033[91m{node_str}\033[0m"  # ANSI escape code for red
        new_highlight = True
    else:
        new_highlight = False

    print(f"{prefix}{node_str}")

    # Get daughters from the Pythia event record.
    # In Pythia, daughter1() and daughter2() return the first and last daughters.
    d1 = particle.daughter1()
    d2 = particle.daughter2()
    daughters = []
    if d1 > 0:
        daughters.append(d1)
    # d2 might be equal to d1 if there's only one daughter.
    if d2 > 0 and d2 != d1:
        daughters.append(d2)

    # Recursively print each daughter with proper branch formatting.
    for i, d in enumerate(daughters):
        branch = "├── " if i < len(daughters) - 1 else "└── "
        print_ascii_daughter_tree(d, pythia, prefix + branch, selected_idx, new_highlight)


# Example usage in main():
#
# After initializing pythia and generating an event, print the daughter trees
# starting from the incoming partons (indices 1 and 2). Optionally, you can
# highlight a particular branch by providing selected_idx.
#
#   print_ascii_daughter_tree(1, pythia, prefix="", selected_idx=42)
#   print_ascii_daughter_tree(2, pythia, prefix="", selected_idx=42)
#
# This will print the ASCII tree of daughters (i.e. the downstream event structure)
# of the earlier indices in the terminal.
        
def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='pythia_charm_track.root', type=str)
	parser.add_argument('--jet-pt-min', help='jet pt min', default=10.0, type=float)
	parser.add_argument('--charm-pt-min', help='charm pt min', default=1.0, type=float)
	parser.add_argument('--debug', help='debug mode', default=False, action='store_true')
	args = parser.parse_args()

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	R = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, R)
	jet_selector= fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_min + 10.) * fj.SelectorAbsEtaMax(1 - R * 1.05)

	parton_jet_R=0.4
	parton_jet_selector= fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_min + 10.) * fj.SelectorAbsEtaMax(1 - R * 1.05)
	parton_jet_def = fj.JetDefinition(fj.antikt_algorithm, parton_jet_R)

	# output
	fout = SingleRootFile(args.output)
	fout.root_file.cd()
	tn_jets = ROOT.TNtuple("tn_jets", "tn_jets", "nev:xsev:ev_weight:code:pt:eta:phi:mass:R:leadpid:has_charm:leadpt:outpid:dRmatch")
	tn_partons = ROOT.TNtuple("tn_partons", "tn_partons", "nev:xsev:ev_weight:code:pid:pt:eta:phi:mass")
	tn_norm   = ROOT.TNtuple('tn_norm', 'tn_norm', 'nev:xsec:xsec_err:sum_of_weights')
	tn_hard 	= ROOT.TNtuple('tn_hard', 'tn_hard', 'nev:xsec:ev_weight:x1:x2:QFac:id1:id2:id3:pt3:eta3:id4:pt4:eta4')
	tn_events = ROOT.TNtuple("tn_events", "tn_events", "nev:xsev:ev_weight:code:out1pid:out2pid:nparts:x1:x2:QFac")
	# tn_events = ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts:x1:x2:QFac:ncoll')

	h_z_charm = ROOT.TH1F('h_z_charm', 'h_z_charm', 10, 0, 1)
	h_z_gsplit = ROOT.TH1F('h_z_gsplit', 'h_z_gsplit', 10, 0, 1)
	h_z_other = ROOT.TH1F('h_z_other', 'h_z_other', 10, 0, 1)
	h_z_total = ROOT.TH1F('h_z_total', 'h_z_total', 10, 0, 1)

	h_zD0_charm = ROOT.TH1F('h_zD0_charm', 'h_zD0_charm', 10, 0, 1)
	h_zD0_gsplit = ROOT.TH1F('h_zD0_gsplit', 'h_zD0_gsplit', 10, 0, 1)
	h_zD0_other = ROOT.TH1F('h_zD0_other', 'h_zD0_other', 10, 0, 1)
	h_zD0_total = ROOT.TH1F('h_zD0_total', 'h_zD0_total', 10, 0, 1)
 
	h_z_charm_log = ROOT.TH1F('h_z_charm_log', 'h_z_charm_log', 10, logbins(0.01, 1.0, 10))
	h_z_gsplit_log = ROOT.TH1F('h_z_gsplit_log', 'h_z_gsplit_log', 10, logbins(0.01, 1.0, 10))
	h_z_other_log = ROOT.TH1F('h_z_other_log', 'h_z_other_log', 10, logbins(0.01, 1.0, 10))
	h_z_total_log = ROOT.TH1F('h_z_total_log', 'h_z_total_log', 10, logbins(0.01, 1.0, 10))

	# configure pythia
	# pThatmin = args.jet_pt_min - 5
	pThatmin = args.charm_pt_min
	pThatmax = 10000.
	mycfg = ['HadronLevel:all = off',
          f'PhaseSpace:pThatMin = {pThatmin}',
          f'PhaseSpace:pThatMax = {pThatmax}',
          'PhaseSpace:bias2Selection = off']
	# disable D0 decay
	mycfg += ['421:onMode = off']
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	# event loop
	if args.nev < 10:
		args.nev = 10
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

		in_parton1 = pythia.event[3]
		in_parton2 = pythia.event[4]
		out_parton1 = pythia.event[5]
		out_parton2 = pythia.event[6]
		fj_out_partons = vector[fj.PseudoJet]()
		for p in [out_parton1, out_parton2]:
			fjp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
			fjp.set_user_index(p.index())
			fj_out_partons.push_back(fjp)

		# jet finding
		parts = vector[fj.PseudoJet]([make_psj(p) for p in pythia.event if p.isFinal() and p.isParton()])
		jets = parton_jet_selector(parton_jet_def(parts))
		jet_labels = [0 for _ in range(len(jets))]
		for i, j in enumerate(jets):
			has_charm = any([abs(pythia.event[p.user_index()].id()) == 4 for p in j.constituents()])
			if has_charm == False:
				continue
			# we have a jet with charm
			leading = fj.sorted_by_pt(j.constituents())[0]
			leading_pythia_part = pythia.event[leading.user_index()]
			# for each charm parton check if it has a mother that is not charm
			charm_in_jet = [p for p in j.constituents() if abs(pythia.event[p.user_index()].id()) == 4]
			if len(charm_in_jet) > 1:
				if args.debug: print('[!] more than one charm in jet')
			for c in charm_in_jet:
				# trace_production_tree(c.user_index(), pythia, level=0)
				# print_ascii_production_tree(c.user_index(), pythia, prefix="-")
				# top_down_production_tree(c.user_index(), pythia)
				# print_ascii_daughter_tree(1, pythia, prefix="-", selected_idx=c.user_index())
				# deep_mother_pid = deep_mother(c.user_index(), pythia)
				deep_mother_pid = deep_parton_mother_pid(c.user_index(), pythia, debug=args.debug)
				print('[=>] deep mother pid', deep_mother_pid)
				# if abs(deep_mother_pid) not in [4, 90]:
				# 	pythia.event.list(False, False, 2)
				# 	return
				h_z_total.Fill(c.pt() / j.pt())
				h_z_total_log.Fill(c.pt() / j.pt())
				jet_labels[i] = 0
				if abs(deep_mother_pid) == 4:
					# print('deep mother is charm')
					h_z_charm.Fill(c.pt() / j.pt())
					h_z_charm_log.Fill(c.pt() / j.pt())
					jet_labels[i] = 4
				if abs(deep_mother_pid) == 21:
					# print('deep mother is gluon')
					h_z_gsplit.Fill(c.pt() / j.pt())
					h_z_gsplit_log.Fill(c.pt() / j.pt())
					jet_labels[i] = 21
				if abs(deep_mother_pid) in [90, 1, 2, 3, 5, 6]:
					# print('deep mother is 0')
					h_z_other.Fill(c.pt() / j.pt())
					h_z_other_log.Fill(c.pt() / j.pt())
					jet_labels[i] = abs(deep_mother_pid)

			# force hadronization and analyze the hadron event
			pythia.forceHadronLevel()
			charged_parts = [p for p in pythia.event if p.isFinal() and (p.isCharged() or abs(p.id()) == 421)]
			# check for D0 mesons
			d0_parts = [p for p in charged_parts if p.isFinal() if abs(p.id()) == 421]
			if len(d0_parts) == 0:
				# print('[!] no D0 mesons found')
				continue
			jets_h = jet_selector(jet_def(vector[fj.PseudoJet]([make_psj(p) for p in charged_parts])))
			for j_h in jets_h:
				d0_in_jet = [p for p in j_h.constituents() if abs(pythia.event[p.user_index()].id()) == 421]
				if len(d0_in_jet) == 0:
					# print('[!] no D0 mesons in jet')
					continue
				for d0 in d0_in_jet:
					h_zD0_total.Fill(d0.pt() / j_h.pt())
					for i, jp in enumerate(jets):
						# if d0.delta_R(jp) < 0.4:
						if j_h.delta_R(jp) < 0.4:
							print('[!] D0 in jet', i, d0.delta_R(jp), jet_labels[i])
							if jet_labels[i] == 4:
								h_zD0_charm.Fill(d0.pt() / j_h.pt())
							if jet_labels[i] == 21:
								h_zD0_gsplit.Fill(d0.pt() / j_h.pt())
							if jet_labels[i] in [90, 1, 2, 3, 5, 6]:
								h_zD0_other.Fill(d0.pt() / j_h.pt())
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
