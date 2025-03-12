#!/usr/bin/env python
# -*- coding: utf-8 -*-

import yasp
yasp.module_load_cppyy('bundle/hepbase')
yasp.module_load_cppyy('heppyy/current')
yasp.module_load_cppyy('alian/current')

import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
Pythia8 = heppyy.load_cppyy('pythia8.Pythia8')
alian_cpp = heppyy.load_cppyy('alian')

import alian

from heppyy.pythia_util import configuration as pyconf

import math
import argparse
import tqdm
import json

class LundJet(object):
	_index = 0

	@classmethod
	def reset_index(cls):
		cls._index = 0

	@classmethod
	def next_index(cls):
		cls._index += 1
		return cls._index
  # here we store the jets cluster sequence

	def __init__(self, jet, parent=None):
		self.jet = jet
		self.parent = parent
		if parent is None:
			self.index = 0
			self.parent_index = -1
			self.reset_index()
		else:
			self.index = self.next_index()
			self.parent_index = parent.index
		self.d1_index = -1
		self.d2_index = -1
		s1 = fj.PseudoJet()
		s2 = fj.PseudoJet()
		if self.jet.has_parents(s1, s2):
			if s1.pt() < s2.pt():
				s1, s2 = s2, s1
			self.d1 = LundJet(s1, self)
			self.d2 = LundJet(s2, self)
			self.d1_index = self.d1.index
			self.d2_index = self.d2.index
		else:
			self.d1 = None
			self.d2 = None

	def __str__(self):
		return self.__repr__()

	def __repr__(self):
		return self.to_string()

	def jet_string(self, indent=0):
		sindent = ''.join(['|-' for j in range(indent)])
		if self.d1 is not None:
			return f'{sindent}[{self.index}->{self.d1_index},{self.d2_index}]<-({self.parent_index}) pt={self.jet.pt():.2f}, eta={self.jet.eta():.2f}, phi={self.jet.phi():.2f}, m={self.jet.m():.2f} '
		else:
			return f'{sindent}[{self.index}]<-({self.parent_index}) pt={self.jet.pt():.2f}, eta={self.jet.eta():.2f}, phi={self.jet.phi():.2f}, m={self.jet.m():.2f} '

	def to_string(self, indent=0):
		if self.d1 is None:
			return self.jet_string(indent)
		else:
			return f'{self.jet_string(indent)}\n{self.d1.to_string(indent+1)}\n{self.d2.to_string(indent+1)}'

	def to_dict(self):
		if self.d1 is None:
			return {'index': self.index, 'parent_index': self.parent_index, 'pt': self.jet.pt(), 'eta': self.jet.eta(), 'phi': self.jet.phi(), 'm': self.jet.m(), 'd1': None, 'd2': None}
		else:
			return {'index': self.index, 'parent_index': self.parent_index, 'pt': self.jet.pt(), 'eta': self.jet.eta(), 'phi': self.jet.phi(), 'm': self.jet.m(), 'd1': self.d1.to_dict(), 'd2': self.d2.to_dict()}

class JetAnalyzer:
	def __init__(self, foutname="pythia_jets_lund.json", pt_min=20.0, jet_R=0.4, abs_rap_max=2.5):
		self.foutname = foutname
		self.pt_min = pt_min
		self.jet_R = jet_R
		self.abs_rap_max = abs_rap_max
		self.jets_dict = {}
		self.jet_index = 0

	def analyze_event(self, pythia):
		"""
		Convert final state Pythia particles to FastJet PseudoJets,
		cluster jets using anti-kt algorithm, and fill histograms.
		"""
		# Collect final state particles.
		py_fj_parts = std.vector(fj.PseudoJet)([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])
		# Define jet clustering: anti-kt algorithm with radius self.jet_R.
		event = pythia.event
		
		# Define jet clustering: anti-kt algorithm with radius self.jet_R.
		jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jet_R)
		cs = fj.ClusterSequence(py_fj_parts, jet_def)
		# Get jets above the pt threshold.
		jets = cs.inclusive_jets(self.pt_min)
		# Apply the rapidity cut
		jet_selector = fj.SelectorAbsRapMax(self.abs_rap_max) * fj.SelectorPtMin(self.pt_min)
		jets = fj.sorted_by_pt(jet_selector(jets))
		jet_def_ca = fj.JetDefinition(fj.cambridge_algorithm, self.jet_R * 2.0)
		# Loop through jets and fill histograms.
		for jet in jets:
			constituents = jet.constituents()
			cs_ca = fj.ClusterSequence(constituents, jet_def_ca)
			jets_ca = fj.sorted_by_pt(jet_selector(cs_ca.inclusive_jets()))
			if jets_ca.size() > 1:
				print('More than one CA jet in the jet ? - skipping')
				continue
			ca_jet = jets_ca[0]
			ljet = LundJet(ca_jet)
			# print(ljet)
			self.jets_dict[self.jet_index] = ljet.to_dict()
			self.jet_index += 1
			break
			
	def save_file(self):
		"""
		Save jet Lund JSON file.
		"""
		print('Saving JSON file to', self.foutname)
		with open(self.foutname, 'w') as fout:
			json.dump(self.jets_dict, fout)

	def __del__(self):
			"""
			Destructor: save histograms before deleting the object.
			"""
			self.save_file()
			

def main():
	argparser = argparse.ArgumentParser(description="Pythia8 jet analysis")
	argparser.add_argument("-o", "--output", default="pythia_jets_lund.json", help="Output ROOT file")
	argparser.add_argument("--ptMin", type=float, default=100.0, help="Minimum jet pT")
	argparser.add_argument("--jetR", type=float, default=0.4, help="Jet radius")
	argparser.add_argument("--nev", type=int, default=1000, help="Number of events")
	argparser.add_argument("--zBias", type=float, default=0.0, help="zBias")
	argparser.add_argument("--coherenceFactor", type=float, default=0.0, help="coherenceFactor")
	argparser.add_argument("--azimuthalBias", type=float, default=0.0, help="azimuthalBias")
	argparser.add_argument("--recoilWeight", type=float, default=0.0, help="recoilWeight")
	argparser.add_argument("--angleBias", type=float, default=0.0, help="angleBias")
	argparser.add_argument("--alphaModifier", type=float, default=0.0, help="alphaModifier")
	
												
	
	args = argparser.parse_args()
	# only in version 8.313
	# current_shower = Pythia8.getShowerModelPtr()
	# print(current_shower)

	myShower = Pythia8.SimpleShowerModelTCustomized()
	myShowerPtr = Pythia8.ShowerModelPtr(myShower)
	#myShower = Pythia8.SimpleTimeShowerCustomized()
	#myShower.getTimeShowerPtr()

	#import pySimpleTimeShower
	#from pySimpleTimeShower import Pythia, SimpleTimeShowerCustomized

	# Create Pythia instance
	pythia = Pythia8.Pythia()
	myShower.addParameters(pythia.settings)
	# Register custom time shower
	# shower = SimpleTimeShowerCustomized()
	pythia.setShowerModelPtr(myShowerPtr)
	extras = ["Next:numberCount = 0", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0", "Stat:showProcessLevel = on"]
	mycfg = [f'PhaseSpace:pThatMin = {args.ptMin}', 'HardQCD:all = on', "Beams:eCM = 13000."]
	mycfg.extend(extras)
	# Initialize with settings
	for c in mycfg:
			pythia.readString(c)
			
	# filepath: /Users/ploskon/devel/alian/alian/sandbox/test_custom_shower.py
	pythia.readString("SimpleTimeShowerCustomized:zBias = 0.0")
	pythia.readString("SimpleTimeShowerCustomized:coherenceFactor = 0.0")
	pythia.readString("SimpleTimeShowerCustomized:azimuthalBias = 0.0")
	pythia.readString("SimpleTimeShowerCustomized:recoilWeight = 0.0")
	pythia.readString("SimpleTimeShowerCustomized:angleBias = 0.0")
	pythia.readString("SimpleTimeShowerCustomized:alphaModifier = 0.0")

	pythia.init()
	if not pythia:
			print("[e] pythia initialization failed.")
			raise RuntimeError("Pythia initialization failed.")

	jana = JetAnalyzer(foutname=args.output, pt_min=args.ptMin, jet_R=args.jetR)
	# Event loop
	# for i in range(args.nev):
	for i in tqdm.tqdm(range(args.nev)):
			if not pythia.next():
					continue
			jana.analyze_event(pythia)

	# Show stats
	pythia.stat()

	pythia.settings.writeFile('test_custom_shower_all.cmnd', True)
	pythia.settings.writeFile('test_custom_shower_changed.cmnd', False)
	
if __name__ == "__main__":
	main()