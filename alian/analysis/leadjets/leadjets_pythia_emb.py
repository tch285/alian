#!/usr/bin/env python

import os
import argparse
import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')

from alian.io.root_io import SingleRootFile
from alian.io import data_io
from alian.utils import data_fj
from alian.steer.glob import globals
from alian.utils.treewriter import RTreeWriter
from alian.analysis.base import BaseAnalysis
from alian.analysis.base import CEventSubtractor
from alian.io.pythia_io import PythiaInput, get_pythia_info, psj_from_particle_with_index
from alian.analysis.leadjets.leadjets_data import JetFinding, gDefaultGridSpacing, gDefaultGridBGyMax

import ROOT

class LeadingJetAnalysis(BaseAnalysis):
	_defaults = {
								'name' : 'leadjets',
              	'jet_Rs': [0.2, 0.4, 0.6],
								'nleading_write': -1,
								'jet_algorithm': fj.antikt_algorithm,
								'part_eta_max': 0.9,
								'bg_y_max': gDefaultGridBGyMax,
								'bg_grid_spacing': gDefaultGridSpacing,
								'write_constituents': False,
        				'pythia_offset_index': 100000,
								'min_jet_pt_signal': 20,
								}

	def __init__(self, **kwargs):
		super(LeadingJetAnalysis, self).__init__(**kwargs)
		# open an output rootfile and create a tree
		_rf = SingleRootFile()
		self.tne = ROOT.TNtuple(f'tne_{self.name}', f'tne_{self.name}', "mult:track_count:centr:bgrho:bgsigma")
		# declare jet analysis for each jet radius
		self.jfs_emb = []
		self.jfs_sig = []
		self.tnj_sig = {}
		self.tnj_mis = {}
		self.tnj_emb = {}
		for R in self.jet_Rs:
			self.jfs_sig.append(JetFinding(jet_R=R, jet_algorithm=self.jet_algorithm, part_eta_max=self.part_eta_max, bg_y_max=self.bg_y_max, bg_grid_spacing=self.bg_grid_spacing))
			self.jfs_emb.append(JetFinding(jet_R=R, jet_algorithm=self.jet_algorithm, part_eta_max=self.part_eta_max, bg_y_max=self.bg_y_max, bg_grid_spacing=self.bg_grid_spacing))
			self.tnj_sig[R] = ROOT.TNtuple(f"tnj_sig_{self.name}_R{R}".replace('.', ''), f"tnj_sig_{self.name}_R{R}".replace('.', ''),  "mult:track_count:centr:bgrho:bgsigma:pt:eta:phi:xsec:xsecErr")
			self.tnj_mis[R] = ROOT.TNtuple(f"tnj_mis_{self.name}_R{R}".replace('.', ''), f"tnj_mis_{self.name}_R{R}".replace('.', ''),  "mult:track_count:centr:bgrho:bgsigma:pt:eta:phi:xsec:xsecErr")
			self.tnj_emb[R] = ROOT.TNtuple(f"tnj_emb_{self.name}_R{R}".replace('.', ''), f"tnj_emb_{self.name}_R{R}".replace('.', ''),  "mult:track_count:centr:bgrho:bgsigma:pt:eta:phi:a:deltaR:deltapt:deltapt_rho:xsec:xsecErr:pt_sig")

		self.part_selector = fj.SelectorAbsEtaMax(self.part_eta_max)
		self.jet_selector = fj.SelectorPtMin(self.min_jet_pt_signal)

	def embed_jet(self, j):
			# print('Embedding jets')
			_rv = std.vector[fj.PseudoJet]()
			for jconst in j.constituents():
					_rv.push_back(jconst)
			for bg in self.psjv:
					_rv.push_back(bg)
			return _rv            

	def pt_match_jets(self, jsig, jemb):
			ptmatch = 0
			for jsigc in jsig.constituents():
					for jembc in jemb.constituents():
							if jsigc.user_index() == jembc.user_index():
									ptmatch += jsigc.pt()
			return ptmatch / jsig.pt()

	def gen_pythia_event(self, e):
			_pythia = None
			ntries = 10000
			# generate a pythia event
			while ntries > 0:
				ntries -= 1
				_pythia = next(e.pythiaio.next_event())
				if _pythia.event.size() > 0:
					self.pythia_parts = self.part_selector(std.vector[fj.PseudoJet]([psj_from_particle_with_index(p, i + self.pythia_offset_index) 
																																					for i, p in enumerate(_pythia.event) if p.isFinal() and p.isVisible() and p.isCharged()]))
					if self.pythia_parts.size() > 0:						
						break
			if ntries == 0:
				RuntimeError('Failed to gen valid pythia event after {} tries'.format(ntries))

	def analyze(self, e):
			self.gen_pythia_event(e)
			_info = get_pythia_info(e.pythiaio.pythia)
			self.parts_emb = std.vector[fj.PseudoJet](e.psjv)
			_ = [self.parts_emb.push_back(p) for p in self.pythia_parts]
			for jf_sig, jf_emb in zip(self.jfs_sig, self.jfs_emb):
				jf_sig.analyze(self.pythia_parts)
				jf_emb.analyze(self.parts_emb)
				e.rho = jf_emb.rho
				e.sigma = jf_emb.sigma
				for j in self.jet_selector(jf_sig.jets):
					miss = True
					for jemb in jf_emb.jets:
						deltapt = self.pt_match_jets(j, jemb)
						if deltapt < 0.5:
							continue
						miss = False
						deltaR = j.delta_R(jemb)
						deltapt_rho = j.perp() - (jemb.perp() - jf_emb.rho * jemb.area())
						self.tnj_emb[jf_emb.jet_R].Fill(e.multiplicity, e.track_count, e.centrality, jf_emb.rho, jf_emb.sigma, jemb.pt(), jemb.eta(), jemb.phi(), jemb.area(), deltaR, deltapt, deltapt_rho, _info.sigmaGen(), _info.sigmaErr(), j.pt())
					if miss:
						self.tnj_mis[jf_sig.jet_R].Fill(e.multiplicity, e.track_count, e.centrality, jf_sig.rho, jf_sig.sigma, j.pt(), j.eta(), j.phi(), _info.sigmaGen(), _info.sigmaErr())
					self.tnj_sig[jf_sig.jet_R].Fill(e.multiplicity, e.track_count, e.centrality, jf_sig.rho, jf_sig.sigma, j.pt(), j.eta(), j.phi(), _info.sigmaGen(), _info.sigmaErr())
			self.tne.Fill(e.multiplicity, e.track_count, e.centrality, e.rho, e.sigma)

	def finalize(self):
		self.rt.write()
		# self.rt.write_and_close()
		# self.root_output.close()


def main():
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results_emb.root")
	parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=data_io.get_default_tree_structure())
	parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)
	parser.add_argument('--nlead', type=int, default=-1, help='how many lead jets to write to the output file - <= 0 is all')
	parser.add_argument('--save-tracks', action='store_true', help='Save track information in the output file')
	parser.add_argument('--cent-min', type=int, help='Minimum centrality', default=-1)
	parser.add_argument('--cent-max', type=int, help='Maximum centrality', default=101)
	parser.add_argument('--no-tqdm', action='store_true', help='Disable tqdm progress bars', default=False)
	parser.add_argument('--cs-dRmax', type=float, help='Max distance for constituent subtraction', default=-1)
	parser.add_argument('--cs-alpha', default=0, type=float)
	parser.add_argument('--part-eta-max', default=0.9, type=float)
	parser.add_argument('--pythia-cmnd', type=str, help='Pythia configuration file', default='$ALIAN_DEV/alian/config/pythia-pp-hardQCD-5TeV-Monash.cmnd')   
	parser.add_argument('--pt-hat-min', type=float, help='Minimum pT hat', default=20.0)

	args = parser.parse_args()
	print(args)

	globals.tqdm_silent = args.no_tqdm

	# open the output file
	root_file = SingleRootFile(fname=args.output)
	print(root_file.root_file.GetName())
 
	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# initialize pythia
	pythia_settings = ['PhaseSpace:pThatMin = {}'.format(args.pt_hat_min), 'PhaseSpace:bias2Selection = on'] #, 'PhaseSpace:bias2SelectionPow = 4', 'PhaseSpace:bias2SelectionRef=50']
	_extra_settings = ['Random:setSeed = on', 'Random:seed = {}'.format(-1), "Next:numberCount = 0", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0", "Stat:showProcessLevel = on"]
	_ = [pythia_settings.append(setting) for setting in _extra_settings]
	pythia_cmnd_file = os.path.expandvars(args.pythia_cmnd)
	pythia_input = PythiaInput(pythia_cmnd_file=pythia_cmnd_file, user_settings=pythia_settings, auto_next=False)

	# get the fj banner out of the way
	fj.ClusterSequence().print_banner()

	cs = None
	if args.cs_dRmax > 0:
		cs = CEventSubtractor(alpha=args.cs_alpha, max_distance=args.cs_dRmax, max_eta=args.part_eta_max, bge_rho_grid_size=gDefaultGridSpacing, max_pt_correct=100)
		print(cs)
	parts_selector = fj.SelectorAbsEtaMax(args.part_eta_max)


	# open the output file
	an = LeadingJetAnalysis(name='ljetana_std', part_eta_max=args.part_eta_max, bg_y_max=gDefaultGridBGyMax, bg_grid_spacing=gDefaultGridSpacing, save_tracks=args.save_tracks, write_constituents=args.save_tracks, nleading_write=args.nlead, min_jet_pt_signal=args.pt_hat_min)
	if args.cs_dRmax > 0:
		an_cs = LeadingJetAnalysis(name='ljetana_cs', part_eta_max=args.part_eta_max, bg_y_max=gDefaultGridBGyMax, bg_grid_spacing=gDefaultGridSpacing, save_tracks=args.save_tracks, write_constituents=args.save_tracks, nleading_write=args.nlead)

	# event loop using the data source directly
	for i,e in enumerate(data_source.next_event()):
		if e.centrality < args.cent_min or e.centrality > args.cent_max:
			continue
		psjv = data_fj.data_tracks_to_pseudojets(e)
		e.psjv = parts_selector(psjv)
		e.pythiaio = pythia_input
		an.analyze(e)
		if cs:
			psjv_cs = cs.process_event(psjv)
			# e.psjv = parts_selector(psjv_cs)
			e.psjv = psjv_cs
			an_cs.analyze(e)
		if i > args.entries and args.entries > 0:
			break

  # finalize the analyses - write the output, close the file, etc.
	# an.finalize()
	# if cs:
	# 	an_cs.analyze(e)
	root_file.close()

if __name__ == '__main__':
		main()
