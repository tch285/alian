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

Optional single-splitting selections (one splitting chosen per jet):
  --max-kt-eec   : use the primary splitting with the largest kT
  --soft-drop-eec: use the first primary splitting satisfying z > --z-sd (default 0.1)

Output
------
  <stem>_eec.parquet         — EEC histograms (all selections), tidy DataFrame
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


# ── single-splitting selectors ────────────────────────────────────────────────

def select_max_kt(lund_seq):
	"""Return the primary Lund splitting with the largest kT, or None."""
	best = None
	for l in lund_seq:
		if best is None or l.kt() > best.kt():
			best = l
	return best


def select_soft_drop(lund_seq, z_cut):
	"""
	Return the first primary Lund splitting satisfying z > z_cut (soft-drop
	condition with beta=0), or None if no splitting passes.
	"""
	for l in lund_seq:
		if l.z() > z_cut:
			return l
	return None


# ── per-jet analysis ──────────────────────────────────────────────────────────

def _fill_splitting(l, accum, jet_pt2):
	"""Fill accum from a single LundDeclustering object."""
	try:
		parts_a = list(l.harder().constituents())
		parts_b = list(l.softer().constituents())
	except Exception:
		return
	accum.fill(parts_a, parts_b, jet_pt2)


def process_jet(jet, lund_gen, eec_accums, lund_accum, kt_cuts, kappa_cuts,
                eec_maxkt=None, eec_sd=None, z_sd=0.1):
	"""
	Run the Lund + E2C analysis for a single jet.

	Parameters
	----------
	jet        : fj.PseudoJet
	lund_gen   : fj.contrib.LundGenerator
	eec_accums : dict { (kt_cut, kappa_cut) -> EECAccum }  threshold-based
	lund_accum : LundPlaneAccum
	kt_cuts    : list of float
	kappa_cuts : list of float
	eec_maxkt  : EECAccum or None   filled with the max-kT splitting
	eec_sd     : EECAccum or None   filled with the soft-drop splitting
	z_sd       : float              soft-drop z threshold (default 0.1)
	"""
	lund_seq = lund_gen.result(jet)
	lund_accum.fill(lund_seq)
	jet_pt2 = jet.pt() ** 2

	# ── threshold-based (all splittings passing cuts) ──────────────────────
	for (kt_cut, kappa_cut), accum in eec_accums.items():
		accum.n_jets += 1
		for l in lund_seq:
			if l.kt()    < kt_cut:    continue
			if l.kappa() < kappa_cut: continue
			_fill_splitting(l, accum, jet_pt2)

	# ── max-kT splitting ───────────────────────────────────────────────────
	if eec_maxkt is not None:
		eec_maxkt.n_jets += 1
		l = select_max_kt(lund_seq)
		if l is not None:
			_fill_splitting(l, eec_maxkt, jet_pt2)

	# ── soft-drop splitting ────────────────────────────────────────────────
	if eec_sd is not None:
		eec_sd.n_jets += 1
		l = select_soft_drop(lund_seq, z_sd)
		if l is not None:
			_fill_splitting(l, eec_sd, jet_pt2)


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
	parser.add_argument('--max-kt-eec',  action='store_true', default=False,
	                    help='Also compute E2C for the max-kT primary splitting per jet')
	parser.add_argument('--soft-drop-eec', action='store_true', default=False,
	                    help='Also compute E2C for the soft-drop selected splitting per jet')
	parser.add_argument('--z-sd',        default=0.1,       type=float,
	                    help='Soft-drop z threshold (default 0.1, used with --soft-drop-eec)')
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
	eec_maxkt  = EECAccum() if args.max_kt_eec    else None
	eec_sd     = EECAccum() if args.soft_drop_eec else None

	if eec_maxkt: print('[i] max-kT EEC enabled')
	if eec_sd:    print(f'[i] soft-drop EEC enabled  (z_sd > {args.z_sd})')

	# event loop
	events = file_events(args.input, max_events=args.max_events) if args.input \
	         else pythia_events(args)

	for parts in events:
		cseq = fj.ClusterSequence(parts, jet_def)
		jets = fj.sorted_by_pt(jet_sel(cseq.inclusive_jets()))
		for jet in jets:
			process_jet(jet, lund_gen, eec_accums, lund_accum, kt_cuts, kappa_cuts,
			            eec_maxkt=eec_maxkt, eec_sd=eec_sd, z_sd=args.z_sd)

	# ── output ────────────────────────────────────────────────────────────────
	stem = os.path.splitext(args.output)[0]

	frames = []
	for (kt_cut, kappa_cut), accum in eec_accums.items():
		accum.normalise()
		lbl = f'{args.label}_kt>{kt_cut}_kappa>{kappa_cut}' if args.label \
		      else f'kt>{kt_cut}_kappa>{kappa_cut}'
		df = accum.to_df(label=lbl, kt_cut=kt_cut)
		df['kappa_cut']  = kappa_cut
		df['selection']  = 'threshold'
		frames.append(df)

	if eec_maxkt is not None:
		eec_maxkt.normalise()
		lbl = f'{args.label}_max_kt' if args.label else 'max_kt'
		df  = eec_maxkt.to_df(label=lbl)
		df['kt_cut']    = float('nan')
		df['kappa_cut'] = float('nan')
		df['selection'] = 'max_kt'
		frames.append(df)

	if eec_sd is not None:
		eec_sd.normalise()
		lbl = f'{args.label}_soft_drop_z>{args.z_sd}' if args.label \
		      else f'soft_drop_z>{args.z_sd}'
		df  = eec_sd.to_df(label=lbl)
		df['kt_cut']    = float('nan')
		df['kappa_cut'] = float('nan')
		df['selection'] = f'soft_drop_z>{args.z_sd}'
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
