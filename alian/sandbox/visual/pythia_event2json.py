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

	def accept_event(self, pythia):
		"""
		Convert final state Pythia particles to FastJet PseudoJets,
		cluster jets using anti-kt algorithm - return True if there is a jet above the pt threshold.
		"""
		# Collect final state particles.
		self.py_fj_parts = std.vector(fj.PseudoJet)([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])
		# reset user index according to pythia event
		indexes = [p.index() for idx, p in enumerate(pythia.event) if p.isFinal()]
		# zip the indexes with the particles
		_ = [p.set_user_index(idx) for idx, p in zip(indexes, self.py_fj_parts)]
		# Define jet clustering: anti-kt algorithm with radius self.jet_R.
		event = pythia.event
		
		# Define jet clustering: anti-kt algorithm with radius self.jet_R.
		self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jet_R)
		self.cs = fj.ClusterSequence(self.py_fj_parts, self.jet_def)
		# set the rapidity cut and minimum pt cut
		jet_selector = fj.SelectorAbsRapMax(self.abs_rap_max) * fj.SelectorPtMin(self.pt_min)
		# Get jets above the pt threshold.
		self.jets = fj.sorted_by_pt(jet_selector(self.cs.inclusive_jets(self.pt_min)))
		# Loop over jets
		for jet in self.jets:
			if jet.pt() > self.pt_min:
				return True
		return False	

def build_particle_tree_1(event):
    def add_particle_to_tree(particle_id, tree):
        particle = event[particle_id]
        if particle_id in tree:
            # If particle is already in tree, update its information
            node = tree[particle_id]
        else:
            # Create a new node for this particle
            node = {
                'id': particle_id,
                'pdgId': particle.id(),
                'status': particle.status(),
                'momentum': (particle.px(), particle.py(), particle.pz(), particle.e()),
                'daughters': {}
            }
            tree[particle_id] = node
        
        # Link daughters recursively
        daughter1 = particle.daughter1()
        daughter2 = particle.daughter2()
        if daughter1 > 0 and daughter2 > 0:
            for d in range(daughter1, daughter2 + 1):
                node['daughters'][d] = add_particle_to_tree(d, tree)
        
        return node
    
    tree = {}
    for i in range(event.size()):
        particle = event[i]
        if particle.mother1() == 0 and particle.mother2() == 0:
            # This is a primary particle (i.e., without parents)
            tree[i] = add_particle_to_tree(i, tree)
    
    return tree
			
def build_particle_tree(event):
    # Step 1: Create all particles in the tree
    tree = {}
    for i in range(event.size()):
        particle = event[i]
        tree[i] = {
            'id': i,
            'pdgId': particle.id(),
            'status': particle.status(),
						'isFinal' : particle.isFinal(),
						'isParton' : particle.isParton(),
						'isHadron' : particle.isHadron(),
						'isLepton' : particle.isLepton(),
						'isQuark' : particle.isQuark(),
						'isGluon' : particle.isGluon(),
						'isNeutral' : particle.isNeutral(),
						'isCharged' : particle.isCharged(),
						'name' : particle.name(),
            'momentum': (particle.px(), particle.py(), particle.pz(), particle.e()),
            'mothers': [],
            'daughters': {}
        }
    
    # Step 2: Link mothers and daughters
    for i in range(event.size()):
        particle = event[i]
        mother1 = particle.mother1()
        mother2 = particle.mother2()
        
        if mother1 > 0:
            tree[mother1]['daughters'][i] = tree[i]
            tree[i]['mothers'].append(mother1)
        if mother2 > 0 and mother2 != mother1:
            tree[mother2]['daughters'][i] = tree[i]
            tree[i]['mothers'].append(mother2)
    
    # Step 3: Return the full tree (not just top-level particles)
    return tree
  
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

	# pythia.readString("HadronLevel:all = off")
			
	# filepath: /Users/ploskon/devel/alian/alian/sandbox/test_custom_shower.py
	# pythia.readString("SimpleTimeShowerCustomized:zBias = 0.0")
	# pythia.readString("SimpleTimeShowerCustomized:coherenceFactor = 0.0")
	# pythia.readString("SimpleTimeShowerCustomized:azimuthalBias = 0.0")
	# pythia.readString("SimpleTimeShowerCustomized:recoilWeight = 0.0")
	# pythia.readString("SimpleTimeShowerCustomized:angleBias = 0.0")
	# pythia.readString("SimpleTimeShowerCustomized:alphaModifier = 0.0")

	pythia.init()
	if not pythia:
			print("[e] pythia initialization failed.")
			raise RuntimeError("Pythia initialization failed.")

	jana = JetAnalyzer(foutname=args.output, pt_min=args.ptMin, jet_R=args.jetR)
	# Event loop
	# for i in range(args.nev):
	lund_jets = {}
	for i in tqdm.tqdm(range(args.nev)):
			if not pythia.next():
					continue
			if jana.accept_event(pythia):
				for i, j in enumerate(jana.jets):
					lund_jets[i] = LundJet(j)
					print(lund_jets[i].to_string())
				break

	tree = build_particle_tree(pythia.event)
	with open("pythia_jets_lund.json", "w") as f:
		json.dump(tree, f)

	# Show stats
	pythia.stat()

	pythia.settings.writeFile('test_custom_shower_all.cmnd', True)
	pythia.settings.writeFile('test_custom_shower_changed.cmnd', False)
	
if __name__ == "__main__":
	main()