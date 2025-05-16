#!/usr/bin/env python

import yaml
import os
import heppyy
import yasp
from yasp import GenericObject

fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
# Pythia8 = heppyy.load_cppyy('pythia8.Pythia8')
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy
from cppyy.gbl import Pythia8

if yasp.in_jupyter_notebook():
	from tqdm.notebook import tqdm	
	print("[i-PythiaInput] Running in Jupyter Notebook, using tqdm.notebook")
else:
	from tqdm import tqdm

def psj_from_particle_with_index(particle, index):
	psj = fj.PseudoJet(particle.px(), particle.py(), particle.pz(), particle.e())
	psj.set_user_index(index)
	return psj

def get_pythia_info(pythia):
	_info = Pythia8.getInfo(pythia)
	return _info

class PythiaInput(GenericObject):
	def __init__(self, pythia_cmnd_file = None, user_settings=[], **kwargs):
		super(PythiaInput, self).__init__(**kwargs)
		self.cmnd_files = []
		if pythia_cmnd_file is not None:
			if isinstance(pythia_cmnd_file, str):
				if pythia_cmnd_file.endswith(".txt"):
					with open(pythia_cmnd_file, "r") as file:
						_cmnd_files = file.readlines()
					self.cmnd_files = [x.strip() for x in _cmnd_files]
				else:
					self.cmnd_files = [pythia_cmnd_file]
		self.user_settings = []
		for fn in self.cmnd_files:
			if os.path.exists(fn):
				with open(fn, 'r') as file:
					_ = [self.user_settings.append(line.strip()) for line in file.readlines()]
				break
		for setting in user_settings:
			self.user_settings.append(setting)
		if self.name is None:
			self.name = "PythiaIO"
		self.event = None
		self.event_count = 0
		self.initialize()

	def initialize(self):
		self.pythia = Pythia8.Pythia()
		for setting in self.user_settings:
			self.pythia.readString(setting)
		self.init = False
		if self.pythia.init():
			print("[i-PythiaInput] Pythia initialized successfully.")
			self.init = True
		else:
			raise RuntimeError("[i-PythiaInput] Pythia initialization failed.")
		if self.n_events is None:
			self.n_events = -1
	
	# Efficiently iterate over the tree as a generator
	def next_event_obsolete(self):
		self.event = None
		pbar_total = None
		if self.n_events > 0:
			pbar_total = tqdm(total=self.n_events, desc="Total events")
		while self.event_count < self.n_events or self.n_events < 0:
			n_try = 100
			gen_ok = False
			while n_try > 0:
				if self.pythia.next():
					gen_ok = True
					n_try = 0
				n_try -= 1
			if gen_ok:
				self.event = self.pythia
				self.event_count += 1
				if pbar_total is not None:
					pbar_total.update(1)
				yield self.event
			else:
				break

	def next_event(self):
		ntry = 100
		gen_ok = False
		while ntry > 0:
			if self.pythia.next():
				gen_ok = True
				ntry = 0
				return gen_ok
			ntry -= 1
		return gen_ok

	def __del__(self):
		if self.pythia and self.init:
			self.pythia.stat()