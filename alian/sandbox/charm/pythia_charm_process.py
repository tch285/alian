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

def deep_charm_mother_pid(idx, pythia):
	# find the mother of the particle with index idx
	mother_idx = pythia.event[idx].mother1()
	mother_pid = pythia.event[mother_idx].id()
	mother2_idx = pythia.event[idx].mother2()
	mother2_pid = pythia.event[mother2_idx].id()
	# follow charm if one of the mothers is charm
	if abs(mother2_idx) > 0:
		if abs(mother2_pid) == 4:
			mother_pid = mother2_pid
			mother_idx = mother2_idx
	# if mother is outgpong parton
	# if mother_idx in [5, 6]:
	if mother_idx in [1, 2, 3, 4, 5, 6, 90]:
		return mother_pid
	else:
		mother_pid = deep_charm_mother_pid(mother_idx, pythia)
	return mother_pid

########### OUTPUT

class Output(GenericObject):
	def __init__(self, foutname, pythia=None):
		self.pythia = pythia
		self.fout = SingleRootFile(foutname)
		self.fout.root_file.cd()
		self.tn_norm = ROOT.TNtuple('tn_norm', 'tn_norm', 'nev:xsec:xsec_err:sum_of_weights')
		self.tn_parton = ROOT.TNtuple('tn_parton', 'tn_parton', 'mother_pid:code:weight:pt:eta:phi:jpt:jeta:jphi')
		self.tn_D0 = ROOT.TNtuple('tn_D0', 'tn_D0', 'mother_pid:code:weight:pt:eta:phi:jpt:jeta:jphi')

		self.no_pid = 777

		self.f_hists_z = {}
		self.f_hists_pt = {}

		self.D0_hists_z = {}
		self.D0_hists_pt = {}
  
		self.all_hists = []

	def get_hists(self, motherpid):
		if self.pythia:
			pythia_info = Pythia8.getInfo(self.pythia)
			self.process = pythia_info.code()
		abs_motherpid = abs(motherpid)
		if self.process not in self.f_hists_z:
			self.f_hists_z[self.process] = {}
			self.f_hists_pt[self.process] = {}
			self.D0_hists_z[self.process] = {}
			self.D0_hists_pt[self.process] = {}
		if abs_motherpid not in self.f_hists_z[self.process]:
			self.fout.root_file.cd()
			self.f_hists_z[self.process][abs_motherpid] = ROOT.TH1F(f'h_z_{self.process}_{abs_motherpid}', f'h_z_{self.process}_{abs_motherpid}', 10, 0, 1)
			self.f_hists_pt[self.process][abs_motherpid] = ROOT.TH1F(f'h_pt_{self.process}_{abs_motherpid}', f'h_pt_{self.process}_{abs_motherpid}', 100, 0, 100)
			self.D0_hists_z[self.process][abs_motherpid] = ROOT.TH1F(f'h_D0_z_{self.process}_{abs_motherpid}', f'h_D0_z_{self.process}_{abs_motherpid}', 10, 0, 1)
			self.D0_hists_pt[self.process][abs_motherpid] = ROOT.TH1F(f'h_D0_pt_{self.process}_{abs_motherpid}', f'h_D0_pt_{self.process}_{abs_motherpid}', 100, 0, 100)
		return [self.f_hists_z[self.process][abs_motherpid], self.f_hists_pt[self.process][abs_motherpid], self.D0_hists_z[self.process][abs_motherpid], self.D0_hists_pt[self.process][abs_motherpid]]

	def fill_parton_hists(self, motherpid, part, jet, weight=1.0):
		hists = self.get_hists(motherpid)
		h_z = hists[0]
		h_pt = hists[1]
		h_z.Fill(part.pt() / jet.pt(), weight)
		h_pt.Fill(part.pt(), weight)
		if motherpid != self.no_pid:
			self.fill_parton_hists(self.no_pid, part, jet)
		self.tn_parton.Fill(motherpid, self.process, weight, part.pt(), part.eta(), part.phi(), jet.pt(), jet.eta(), jet.phi())

	def fill_D0_hists(self, motherpid, part, jet, weight=1.0):
		hists = self.get_hists(motherpid)
		h_z = hists[2]
		h_pt = hists[3]
		h_z.Fill(part.pt() / jet.pt(), weight)
		h_pt.Fill(part.pt(), weight)
		if motherpid != self.no_pid:
			self.fill_D0_hists(self.no_pid, part, jet)
		self.tn_D0.Fill(motherpid, self.process, weight, part.pt(), part.eta(), part.phi(), jet.pt(), jet.eta(), jet.phi())

	def set_norm(self, pythia):
		pythia_info = Pythia8.getInfo(pythia)
		xsec = pythia_info.sigmaGen()
		xsec_err = pythia_info.sigmaErr()
		sum_weights = pythia_info.weightSum()  # Same as sum_weights
		nev = pythia_info.nAccepted()
		self.tn_norm.Fill(nev, xsec, xsec_err, sum_weights)
		for process in self.f_hists_z:
			for abs_motherpid in self.f_hists_z[process]:
				hs = self.get_hists(abs_motherpid)
				for h in hs:
					h.Scale(xsec / sum_weights)
					self.all_hists.append(h)

	def close(self):
		print('closing output', self.fout.root_file.GetName())
		self.fout.close()

	def __delete__(self):
		self.close()

########### MY PYTHIA INFO

class MyPythiaInfo(GenericObject):
	def __init__(self, pythia):
		self.pythia = pythia	
		self.sum_weights = 0.0

	def update(self):
		pythia_info = Pythia8.getInfo(self.pythia)
		self.sigma = pythia_info.sigmaGen()
		self.ev_weight = pythia_info.weight()
		self.process = pythia_info.code()
		self.sum_weights += self.ev_weight
		self.abs_fs_part_codes = [abs(p.id()) for p in self.pythia.event if p.isFinal()]
		self.in_parton1 = self.pythia.event[3]
		self.in_parton2 = self.pythia.event[4]
		self.out_parton1 = self.pythia.event[5]
		self.out_parton2 = self.pythia.event[6]

	def has_charm(self):
		self.has_charm_flag = False
		if 4 in self.abs_fs_part_codes:
			self.has_charm_flag = True
		return self.has_charm_flag

	def charm_above_pt(self, pt_min):
		# charm above a pt
		self.charm_above_pt_list = [p for p in self.pythia.event if p.isFinal() and abs(p.id()) == 4 and p.pT() > pt_min and abs(p.eta()) < 1]
		if len(self.charm_above_pt_list) == 0:
			return False
		return True

	def get_out_partons(self):
		self.fj_out_partons = vector[fj.PseudoJet]()
		for p in [self.out_parton1, self.out_parton2]:
			fjp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
			fjp.set_user_index(p.index())
			self.fj_out_partons.push_back(fjp)
		return self.fj_out_partons

	def get_fs_partons_psj(self):
		self.fs_partons = vector[fj.PseudoJet]([make_psj(p) for p in self.pythia.event if p.isFinal() and p.isParton()])
		return self.fs_partons

	def write_info_json(self, foutname):
   		# Save to JSON summary
		# Get the total cross section and weight sum
		pythia_info = Pythia8.getInfo(self.pythia)
		sigma_gen = pythia_info.sigmaGen()
		sigma_gen_err = pythia_info.sigmaErr()
		weight_sum = pythia_info.weightSum()  # Same as sum_weights
		nAccepted = pythia_info.nAccepted()
		json_file = foutname.replace('.root', '.json')
		with open(json_file, "w") as f:
				json.dump({
						"n_accepted": nAccepted,
						"sigma_gen": sigma_gen,
						"sigma_gen_err": sigma_gen_err,
						"sum_weights": weight_sum
				}, f, indent=2)

########### MAIN PROGRAM

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='pythia_charm_process.root', type=str)
	parser.add_argument('--jet-pt-min', help='jet pt min', default=10.0, type=float)
	parser.add_argument('--charm-pt-min', help='charm pt min', default=0.1, type=float)
	parser.add_argument('--debug', help='debug mode', default=False, action='store_true')
	parser.add_argument('--D0required', help='require D0 in jet', default=False, action='store_true')
	args = parser.parse_args()

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	R = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, R)
	# jet_selector= fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_min + 10.) * fj.SelectorAbsEtaMax(1 - R * 1.05)
	jet_selector= fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(1 - R * 1.05)
	parton_jet_R = 0.4
	# parton_jet_selector= fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_min + 10.) * fj.SelectorAbsEtaMax(1 - R * 1.05)
	parton_jet_selector= fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(1 - R * 1.05)
	parton_jet_def = fj.JetDefinition(fj.antikt_algorithm, parton_jet_R)

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

	pinfo = MyPythiaInfo(pythia)


	# output
	file_output = Output(args.output, pythia)
	print(f'output file: {args.output}, {file_output.fout.root_file.GetName()}')

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

		pinfo.update()
		if pinfo.has_charm() is False:
			continue
		if pinfo.charm_above_pt(args.charm_pt_min) is False:
			continue

		out_partons = pinfo.get_out_partons()
		psj_fs_partons = pinfo.get_fs_partons_psj()		

		jets = parton_jet_selector(parton_jet_def(psj_fs_partons))
		jet_labels = [0 for _ in range(len(jets))]

		for i, j in enumerate(jets):
			charm_in_jet = [p for p in j.constituents() if abs(pythia.event[p.user_index()].id()) == 4]
			has_charm = any(charm_in_jet)
			if has_charm == False:
				continue
			if len(charm_in_jet) > 1:
				if args.debug: print('[!] more than one charm in jet')

			# now we have a jet with charm
			leading_part = fj.sorted_by_pt(j.constituents())[0]
			leading_part_pid = pythia.event[leading_part.user_index()].id()
			# for each charm parton check if it has a mother that is not charm
			for c in charm_in_jet:
				deep_mother_pid = deep_charm_mother_pid(c.user_index(), pythia)
				# print('[=>] deep mother pid' , deep_mother_pid) #, 'fout', file_output)
				file_output.fill_parton_hists(deep_mother_pid, c, j, weight=pinfo.ev_weight)
				jet_labels[i] = deep_mother_pid
			if not args.D0required:
				accepted += 1
				pbar.update(1)
				continue

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
					jet_label = file_output.no_pid
					for i, jp in enumerate(jets):
						if j_h.delta_R(jp) < 0.4:
							jet_label = jet_labels[i]
							print('[!] D0 in jet', i, d0.delta_R(jp), jet_labels[i])
							break
					file_output.fill_D0_hists(jet_label, d0, j_h, weight=pinfo.ev_weight)
				if args.D0required:
					accepted += 1
					pbar.update(1)

	pythia.stat()
	pinfo.write_info_json(args.output)
	# Output for verification
	file_output.set_norm(pythia)
	file_output.close()

if __name__ == '__main__':
	main()
