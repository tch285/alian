#!/usr/bin/env python

import yaml
import os
import heppyy
from tqdm import tqdm
from yasp import GenericObject

fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
Pythia8 = heppyy.load_cppyy('pythia8.Pythia8')

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
			print("Pythia initialized successfully.")
			self.init = True
		else:
			raise RuntimeError("Pythia initialization failed.")
	
	# Efficiently iterate over the tree as a generator
	def next_event(self):
		self.event = None
		pbar_total = None
		if self.n_events > 0:
			pbar_total = tqdm(total=self.n_events, desc="Total events")
		while self.event_count < self.n_events:
			if self.pythia.next():
				self.event = self.pythia
				yield self.event
				self.event_count += 1
				if pbar_total is not None:
					pbar_total.update(1)
			else:
				break