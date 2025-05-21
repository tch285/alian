#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
	
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import pandas as pd

def main(args):
	# generate a sample of events with n=1000 particles
	n = args.mult
	dndeta = args.dndeta
	mean_pt = args.mean_pt
	nevents = args.events
	events = {}
	for i in range(nevents):
		# Generate n particles for each event
		if n > 0:
			events[i] = generate_boltzmann_particles_mult(n, mean_pt)
		else:
			# Use dndeta to calculate the number of particles
			# n = int(dndeta * (eta_range[1] - eta_range[0]))
			events[i] = generate_boltzmann_particles_dndeta(dndeta, mean_pt)

	df = pd.DataFrame(events)	

	all_particles = {}
	# make a flat table of all particles
	for i in range(nevents):
		# Flatten the dictionary and append to all_particles
		for key in events[i].keys():
			if key not in all_particles:
				all_particles[key] = []
			all_particles[key].extend(events[i][key])

	all_particles['weights'] = [1/nevents] * len(all_particles['pT'])
	# Convert to DataFrame for easier plotting
	all_particles = pd.DataFrame(all_particles)

	fig, axs = plt.subplots(1, 2, figsize=(12, 5))

	# pT histogram
	sns.histplot(data=all_particles, x='pT', bins=50, kde=True, color='blue', weights='weights', ax=axs[0])
	axs[0].set_title('Transverse Momentum (pT) Distribution')
	axs[0].set_xlabel('Transverse Momentum (pT)')
	axs[0].set_ylabel('Counts')

	# eta histogram
	sns.histplot(all_particles['eta'], bins=50, kde=True, color='green', ax=axs[1])
	axs[1].set_title('Pseudorapidity (η) Distribution')
	axs[1].set_xlabel('Pseudorapidity (η)')
	axs[1].set_ylabel('Counts')

	plt.tight_layout()
	plt.show()
    
if __name__ == "__main__":
	# some numbers in https://arxiv.org/pdf/1612.08966 ~2k dn/dy
	parser = argparse.ArgumentParser(description='Generate Boltzmann distributed particles')
	parser.add_argument('--mult', type=int, default=0, help='Number of particles per event')
	parser.add_argument('--dndeta', type=float, default=2200, help='Number of particles per event')
	parser.add_argument('--mean-pt', type=float, default=0.7, help='Mean transverse momentum (pT)')
	parser.add_argument('--events', type=int, default=1, help='Number of events to generate')
	args = parser.parse_args()
	main(args)
