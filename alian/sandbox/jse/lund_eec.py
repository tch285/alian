#!/usr/bin/env python3
"""
lund_eec.py — Primary Lund plane + subjet-resolved E2C analysis.

Input modes
-----------
  Pythia8  (default, standard pyconf args)
  File     (--input file.root, JEWEL ROOT file via JewelReader)

For each primary Lund splitting that passes a kT and/or kappa threshold,
four E2C distributions are accumulated (binned in ln ΔR):

  eec_AA  : pairs both within harder subjet A
  eec_BB  : pairs both within softer subjet B (the emission)
  eec_AB  : cross pairs — one constituent from A, one from B
  eec_all : all pairs from A ∪ B  (AA + BB + AB, each unordered pair once)

Output
------
  <stem>_eec.parquet         — EEC histograms per kT cut, as a tidy DataFrame
  <stem>_lundplane.npy       — 2-D Lund plane histogram (H[ix, iy])
  <stem>_lundplane_xedges.npy
  <stem>_lundplane_yedges.npy
"""

from __future__ import print_function
import os
import math
import argparse
import numpy as np
import tqdm
import pandas as pd

import heppyy.util.fastjet_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl.std import vector


# ── E2C accumulator ──────────────────────────────────────────────────────────

class EECAccum:
	"""
	Accumulate weighted E2C pairs into ln(ΔR) histograms.

	Convention: each unordered pair (i,j) with i<j is counted once,
	weighted by pT_i * pT_j / pT_jet^2.
	"""

	def __init__(self, nbins=60, lndr_min=-5.0, lndr_max=1.5):
		self.edges   = np.linspace(lndr_min, lndr_max, nbins + 1)
		self.nbins   = nbins
		self.h       = {t: np.zeros(nbins) for t in ('AA', 'BB', 'AB', 'all')}
		self.n_jets  = 0
		self.n_splits = 0

	def fill(self, parts_a, parts_b, jet_pt2):
		"""
		Fill all four E2C types for one Lund splitting.

		Parameters
		----------
		parts_a  : list of fj.PseudoJet  (harder subjet constituents)
		parts_b  : list of fj.PseudoJet  (softer subjet constituents)
		jet_pt2  : float  (jet pT² for normalisation)
		"""
		def _add(p1, p2, key):
			dR = p1.delta_R(p2)
			if dR <= 0:
				return
			idx = np.searchsorted(self.edges, math.log(dR), side='right') - 1
			if 0 <= idx < self.nbins:
				w = p1.pt() * p2.pt() / jet_pt2
				self.h[key][idx]   += w
				self.h['all'][idx] += w

		# AA: pairs within harder subjet
		for i in range(len(parts_a)):
			for j in range(i + 1, len(parts_a)):
				_add(parts_a[i], parts_a[j], 'AA')

		# BB: pairs within softer subjet
		for i in range(len(parts_b)):
			for j in range(i + 1, len(parts_b)):
				_add(parts_b[i], parts_b[j], 'BB')

		# AB: cross pairs
		for pa in parts_a:
			for pb in parts_b:
				_add(pa, pb, 'AB')

		self.n_splits += 1

	def normalise(self):
		if self.n_jets > 0:
			for k in self.h:
				self.h[k] = self.h[k] / self.n_jets

	def to_df(self, label=None, kt_cut=None):
		centers = 0.5 * (self.edges[:-1] + self.edges[1:])
		d = {
			'ln_dR':    centers,
			'bin_lo':   self.edges[:-1],
			'bin_hi':   self.edges[1:],
			'eec_AA':   self.h['AA'],
			'eec_BB':   self.h['BB'],
			'eec_AB':   self.h['AB'],
			'eec_all':  self.h['all'],
			'n_jets':   self.n_jets,
			'n_splits': self.n_splits,
		}
		if label   is not None: d['label']  = label
		if kt_cut  is not None: d['kt_cut'] = kt_cut
		return pd.DataFrame(d)


# ── Lund-plane accumulator ────────────────────────────────────────────────────

class LundPlaneAccum:
	"""2-D histogram of ln(1/Δ) vs ln(kt) for the primary Lund plane."""

	def __init__(self, nbins=50, lninvd_range=(0, 5), lnkt_range=(-5, 3)):
		self.xedges = np.linspace(*lninvd_range, nbins + 1)  # ln(1/Δ)
		self.yedges = np.linspace(*lnkt_range,   nbins + 1)  # ln(kt)
		self.h2     = np.zeros((nbins, nbins))

	def fill(self, lund_seq):
		for l in lund_seq:
			if l.Delta() <= 0 or l.kt() <= 0:
				continue
			xi = np.searchsorted(self.xedges, math.log(1.0 / l.Delta()), side='right') - 1
			yi = np.searchsorted(self.yedges, math.log(l.kt()),           side='right') - 1
			if 0 <= xi < self.h2.shape[0] and 0 <= yi < self.h2.shape[1]:
				self.h2[xi, yi] += 1


# ── per-jet analysis ──────────────────────────────────────────────────────────

def process_jet(jet, lund_gen, eec_accums, lund_accum, kt_cuts, kappa_cuts):
	"""
	Run the Lund + E2C analysis for a single jet.

	Parameters
	----------
	jet        : fj.PseudoJet
	lund_gen   : fj.contrib.LundGenerator
	eec_accums : dict { (kt_cut, kappa_cut) -> EECAccum }
	lund_accum : LundPlaneAccum
	kt_cuts    : list of float   (matched 1-to-1 with eec_accums keys)
	kappa_cuts : list of float
	"""
	lund_seq = lund_gen.result(jet)
	lund_accum.fill(lund_seq)
	jet_pt2 = jet.pt() ** 2

	for (kt_cut, kappa_cut), accum in eec_accums.items():
		accum.n_jets += 1
		for l in lund_seq:
			if l.kt()    < kt_cut:    continue
			if l.kappa() < kappa_cut: continue
			try:
				parts_a = list(l.harder().constituents())
				parts_b = list(l.softer().constituents())
			except Exception:
				continue
			accum.fill(parts_a, parts_b, jet_pt2)


# ── event sources ─────────────────────────────────────────────────────────────

def pythia_events(args):
	"""Yield vector<PseudoJet> per event from Pythia8."""
	import heppyy.util.pythia8_cppyy
	from cppyy.gbl import Pythia8
	from heppyy.pythia_util import configuration as pyconf

	pythia = pyconf.create_and_init_pythia_from_args(args, [])
	if not pythia:
		raise RuntimeError('[e] Pythia initialisation failed.')

	pbar = tqdm.tqdm(total=args.nev)
	while pbar.n < args.nev:
		if not pythia.next():
			continue
		yield vector[fj.PseudoJet]([
			fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
			for p in pythia.event if p.isFinal() and p.isCharged()
		])
		pbar.update(1)
	pythia.stat()


def file_events(file_path, max_events=None):
	"""Yield vector<PseudoJet> per event from a JEWEL ROOT file."""
	import uproot

	with uproot.open(file_path) as f:
		tracks_df     = f['tracks'].arrays(library='pd')
		event_info_df = f['event_info'].arrays(library='pd')
	events_df = pd.merge(tracks_df, event_info_df, on='eventID', how='inner')

	for n, (_, ev) in enumerate(tqdm.tqdm(events_df.groupby('eventID'))):
		if max_events is not None and n >= max_events:
			break
		yield vector[fj.PseudoJet]([
			fj.PseudoJet(px, py, pz, e)
			for px, py, pz, e in zip(ev['px'], ev['py'], ev['pz'], ev['energy'])
		])


# ── main ──────────────────────────────────────────────────────────────────────

def main():
	parser = argparse.ArgumentParser(
		description='Lund plane + subjet E2C analysis',
		prog=os.path.basename(__file__))

	try:
		from heppyy.pythia_util import configuration as pyconf
		pyconf.add_standard_pythia_args(parser)
	except Exception:
		pass

	parser.add_argument('--input',       default=None,
	                    help='JEWEL ROOT file (enables file mode, overrides Pythia)')
	parser.add_argument('--max-events',  default=None,      type=int)
	parser.add_argument('--output',      default='lund_eec_output.parquet', type=str)
	parser.add_argument('--jetR',        default=0.4,       type=float)
	parser.add_argument('--jet-pt-min',  default=100.0,     type=float)
	parser.add_argument('--jet-pt-max',  default=120.0,     type=float)
	parser.add_argument('--etadet',      default=2.5,       type=float)
	parser.add_argument('--kt-cuts',     default='0,1,2,5', type=str,
	                    help='Comma-separated kT minima [GeV] (applied to each Lund splitting)')
	parser.add_argument('--kappa-cuts',  default='0',       type=str,
	                    help='Comma-separated kappa minima (matched to --kt-cuts, padded with 0)')
	parser.add_argument('--label',       default='',        type=str)
	args = parser.parse_args()

	# jet finder + Lund generator
	fj.ClusterSequence.print_banner()
	jetR     = args.jetR
	jet_def  = fj.JetDefinition(fj.antikt_algorithm, jetR)
	jet_sel  = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(args.etadet - jetR * 1.05)
	if args.jet_pt_max > 0:
		jet_sel *= fj.SelectorPtMax(args.jet_pt_max)
	lund_gen = fj.contrib.LundGenerator()
	print('[i] LundGenerator:', lund_gen)

	# build cut pairs
	kt_cuts    = [float(x) for x in args.kt_cuts.split(',')]
	kappa_vals = [float(x) for x in args.kappa_cuts.split(',')]
	kappa_cuts = [kappa_vals[i] if i < len(kappa_vals) else 0.0 for i in range(len(kt_cuts))]

	eec_accums = {(kt, kp): EECAccum() for kt, kp in zip(kt_cuts, kappa_cuts)}
	lund_accum = LundPlaneAccum()

	# event loop
	events = file_events(args.input, max_events=args.max_events) if args.input \
	         else pythia_events(args)

	for parts in events:
		cseq = fj.ClusterSequence(parts, jet_def)
		jets = fj.sorted_by_pt(jet_sel(cseq.inclusive_jets()))
		for jet in jets:
			process_jet(jet, lund_gen, eec_accums, lund_accum, kt_cuts, kappa_cuts)

	# ── output ────────────────────────────────────────────────────────────────
	stem = os.path.splitext(args.output)[0]

	frames = []
	for (kt_cut, kappa_cut), accum in eec_accums.items():
		accum.normalise()
		lbl = f'{args.label}_kt>{kt_cut}_kappa>{kappa_cut}' if args.label \
		      else f'kt>{kt_cut}_kappa>{kappa_cut}'
		df = accum.to_df(label=lbl, kt_cut=kt_cut)
		df['kappa_cut'] = kappa_cut
		frames.append(df)

	eec_df   = pd.concat(frames, ignore_index=True)
	eec_path = stem + '_eec.parquet'
	eec_df.to_parquet(eec_path, engine='pyarrow')
	print(f'[i] EEC histograms  -> {eec_path}')

	lund_path = stem + '_lundplane.npy'
	np.save(lund_path,                       lund_accum.h2)
	np.save(stem + '_lundplane_xedges.npy',  lund_accum.xedges)
	np.save(stem + '_lundplane_yedges.npy',  lund_accum.yedges)
	print(f'[i] Lund plane      -> {lund_path}')


if __name__ == '__main__':
	main()
