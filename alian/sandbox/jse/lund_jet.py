#!/usr/bin/env python
# -*- coding: utf-8 -*-

from yasp import GenericObject
import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from heppyy.util.mputils import logbins
from heppyy.pythia_util import configuration as pyconf

import argparse
import os 
import math
import tqdm
import pandas as pd
import numpy as np

def generate_boltzmann_particles_mult(n, mean_pt, mass=0.139, eta_range=(-2, 2)):
		"""
		Generate n particles with transverse momenta (pT) drawn from a Boltzmann distribution
		with a given mean pT. Optionally specify the particle mass (default: pion mass in GeV).
		Returns a list of dicts with px, py, pz, E for each particle.
		"""
		# The Boltzmann distribution for pT: f(pT) ~ pT * exp(-pT/T)
		# The mean pT = 2*T, so T = mean_pt / 2
		# T = mean_pt / 2.0
		T = mean_pt
		# Sample pT
		pT = np.random.exponential(scale=T, size=n)
		# Sample phi uniformly
		phi = np.random.uniform(0, 2*np.pi, size=n)
		px = pT * np.cos(phi)
		py = pT * np.sin(phi)
		# Sample eta uniformly in some range (e.g., -1 to 1)
		eta = np.random.uniform(eta_range[0], eta_range[1], size=n)
		pz = pT * np.sinh(eta)
		E = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
		# particles = [{'px': px[i], 'py': py[i], 'pz': pz[i], 'E': E[i] , 'pT': pT[i], 'phi': phi[i], 'eta': eta[i]} for i in range(n)]
		particles = {'px': px, 'py': py, 'pz': pz, 'E': E , 'pT': pT, 'phi': phi, 'eta': eta}
		return particles

def generate_boltzmann_particles_dndeta(dndeta, mean_pt, mass=0.139, eta_range=(-2, 2)):
	n = int(dndeta * (eta_range[1] - eta_range[0]))
	return generate_boltzmann_particles_mult(n, mean_pt, mass, eta_range)

def generate_boltzmann_pseudojets(dndeta, mean_pt, mass=0.139, eta_range=(-2, 2)):
	ps = generate_boltzmann_particles_dndeta(dndeta, mean_pt, mass, eta_range)
	psjv = vector[fj.PseudoJet]([fj.PseudoJet(ps['px'][i], ps['py'][i], ps['pz'][i], ps['E'][i]) for i in range(len(ps['px']))])
	return psjv

# make a singleton class for JetAlgoHelper
class JetAlgoHelper(object):
	_instance = None

	@staticmethod
	def get_instance():
		if JetAlgoHelper._instance is None:
			JetAlgoHelper()
		return JetAlgoHelper._instance

	def __init__(self):
		if JetAlgoHelper._instance is not None:
			raise Exception("This class is a singleton!")
		else:
			JetAlgoHelper._instance = self
			self.jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
			self.jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
			self.reclusterer_wta 	= fj.contrib.Recluster(self.jet_def_wta)
			self.sd01 = fj.contrib.SoftDrop(0, 0.1, 1.0)
			self.sd02 = fj.contrib.SoftDrop(0, 0.2, 1.0)
			self.lund_gen = fj.contrib.LundGenerator()
			print('[i] creating LundGenerator:', self.lund_gen)

	@classmethod
	def angularity(self, jet, a, k, jetR):
		ang = 0.0
		for p in jet.constituents():
			dr = jet.delta_R(p) / jetR
			pt = p.perp() / jet.perp()
			ang += ((dr)**a) * ((pt)**k)
		return ang

	def mass(self, jet):
		m2 = jet.e()**2 - jet.px()**2 - jet.py()**2 - jet.pz()**2
		if m2 > 0:
			return math.sqrt(m2)
		return 0.0

	def lund_delta_kt(self, jet):
		return [[l.Delta(), l.kt()] for l in self.lund_gen.result(jet)]

	def lund_log(self, jet):
		return [[math.log(1./l.Delta()), math.log(l.kt())] for l in self.lund_gen.result(jet)]

	def lunds_dict_list(self, jet):
		lunds = []
		for i, l in enumerate(self.lund_gen.result(jet)):
			lunds.append({'i': i, 'pt': l.pair().perp(), 
                 		'pt1': l.harder().perp(), 'pt2': l.softer().perp(), 'eta': l.pair().eta(), 
                   	'kt': l.kt(), 'delta': l.Delta(), 'kappa': l.kappa(), 'psi': l.psi(), 'z': l.z(), 'm': l.m()})
		return lunds

class LundJet(GenericObject):
	def __init__(self, jet, jetR, label=None, **kwargs):
		super().__init__(**kwargs)
		self.jet = jet
		self.pt = jet.pt()
		self.eta = jet.eta()
		self.y = jet.rap()
		self.phi = jet.phi()
		self.e = jet.e()
		self.m = jet.m()
		self.nconst = jet.constituents().size()
		self.jetR = jetR
		self.label = label

		self._jalgo = JetAlgoHelper.get_instance()

		self.jet_wta 	= self._jalgo.reclusterer_wta.result(jet)
		self.jet_sd01 = self._jalgo.sd01.result(jet)
		self.jet_sd02 = self._jalgo.sd02.result(jet)
		self.wtastd 	= self.jet_wta.delta_R(jet)
		self.wtasd01 	= self.jet_sd01.delta_R(jet)
		self.wtasd02 	= self.jet_sd02.delta_R(jet)
		self.angk1a1 	= self._jalgo.angularity(jet, 1.0, 1.0, self.jetR)
		self.angk1a2 	= self._jalgo.angularity(jet, 2.0, 1.0, self.jetR)
		self.angk1a3 	= self._jalgo.angularity(jet, 3.0, 1.0, self.jetR)
		self.mjet 		= self._jalgo.mass(jet)
		# self.lund_delta_kt = self._jalgo.lund_delta_kt(jet)
		# self.lund_log = self._jalgo.lund_log(jet)
		self.lunds = self._jalgo.lunds_dict_list(jet)

		self._base_props_list = []
		self._base_props_list = self._gen_base_props_list()

	def _gen_base_props_list(self):
		_g0 = GenericObject()
		l = [a for a in _g0.__dict__]
		_ = [l.append(a) for a in self.__dict__ if a[0] == '_' and a not in l]
		return l

	def to_dict(self):
		# Convert the object to a dictionary
		d = {}
		for key in self.__dict__:
			if key not in self._base_props_list:
				d[key] = getattr(self, key)
		return d

	def to_basic_type_dict(self):
		# Convert the object to a dictionary
		d = {}
		for key in self.__dict__:
			if key not in self._base_props_list:
				o = getattr(self, key)
				if isinstance(o, fj.PseudoJet):
					d[key] = [o.px(), o.py(), o.pz(), o.e()]
				else:
					d[key] = getattr(self, key)
		return d


def merge_events(levents, offsets):
	psjv = vector[fj.PseudoJet]()
	for i, e in enumerate(levents):
		for ip in range(e.size()):
			p = fj.PseudoJet(e[ip].px(), e[ip].py(), e[ip].pz(), e[ip].e())
			p.set_user_index(ip + offsets[i])
			psjv.push_back(p)
	return psjv

def match_z(jet, jets):
	for j in jets:
		if jet.delta_R(j) < 0.1:
			# print(f'[i] jet pT={jet.perp()} eta={jet.eta()} matched with pT={j.perp()} eta={j.eta()}')
			return j
	# print(f'[e] no match found for jet with pT={jet.perp()} eta={jet.eta()}')
	return None

def main():

	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('--jet-pt-min', help='jet pt min', default=100.0, type=float)
	parser.add_argument('--jet-pt-max', help='jet pt max', default=120.0, type=float)
	parser.add_argument('--etadet', help='detector eta', default=2.5, type=float)
	parser.add_argument('--shape', help='fill the jet shape histograms', action='store_true', default=False)
	parser.add_argument('--jetR', help='jet radius', default=0.4, type=float)
	parser.add_argument('--output', '-o', help='output file name', default='pythia_lund_jet.parquet', type=str)
	parser.add_argument('--fixed-label', help='label the jets with a fixed value', default=-1, type=int)
	parser.add_argument('--thermal', help='embed jets into a thermal (Boltzmann) background', action='store_true', default=False)
	args = parser.parse_args()

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_def = fj.JetDefinition(fj.antikt_algorithm, args.jetR)
	jet_selector = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(args.etadet - args.jetR * 1.05)
	if args.jet_pt_max > 0:
		jet_selector *= fj.SelectorPtMax(args.jet_pt_max)
	
	jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
	reclusterer_wta =  fj.contrib.Recluster(jet_def_wta)

	sd01 = fj.contrib.SoftDrop(0, 0.1, 1.0)
	sd02 = fj.contrib.SoftDrop(0, 0.2, 1.0)

	bg_event = None

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	if args.nev < 10:
		args.nev = 10
	count_jets = 0
	count_jets_emb = 0

	pbar = tqdm.tqdm(total=args.nev)
	jets_dicts = []
	jets_dicts_emb = []

	while pbar.n < args.nev:
		if not pythia.next():
			continue
		# parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isVisible()])
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)

		jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		if len(jets) > 0:
			pbar.update(1)
		else:
			pbar.update(0)
			continue
		# print('number of jets in the pythia event:', len(jets), 'leading jet pT:', jets[0].perp())
		for j in jets:
			count_jets += 1
			# add the LundPlane calculation
			# add jet and the LundPlane to the output
			lj = LundJet(jet=j, jetR=args.jetR, label=args.fixed_label)
			# lj_dict = lj.to_dict()
			lj_dict = lj.to_basic_type_dict()
			jets_dicts.append(lj_dict)
			# print(lj_dict)
   
		# idea for photons: make jets single particle and embed them into the event ...
		# write to the separate file...

		if args.thermal:
			dndeta = 2200
			mean_pt = 0.7
			bg_event = generate_boltzmann_pseudojets(dndeta, mean_pt, mass=0.139, eta_range=(-args.etadet, args.etadet))
			merged_event = merge_events([parts, bg_event], [0, 10000])
			jets_emb = fj.sorted_by_pt(jet_def(merged_event))
			# print('number of jets in the embedded event:', len(jets_emb), 'leading jet pT:', jets_emb[0].perp())
			for j in jets:
				jmatched = match_z(j, jets_emb)
				if jmatched is None:
					continue
				# add the LundPlane calculation
				count_jets_emb += 1
				lj_emb = LundJet(jet=jmatched, jetR=args.jetR, label=args.fixed_label)
				lj_dict_emb = lj_emb.to_basic_type_dict()
				jets_dicts_emb.append(lj_dict_emb)

	pythia.stat()

	df = pd.DataFrame(jets_dicts)
	print(f'number of jets: {len(jets_dicts)}')
	df.to_parquet(args.output, engine="pyarrow")

	df_emb = pd.DataFrame(jets_dicts_emb)
	print(f'number of jets embedded: {len(jets_dicts_emb)}')
	df_emb.to_parquet(args.output.replace('.parquet', '_emb.parquet'), engine="pyarrow")

if __name__ == "__main__":
	main()