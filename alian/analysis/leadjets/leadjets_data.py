#!/usr/bin/env python

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

import ROOT

# this is finding jets - not doing anything with them
class JetFinding(BaseAnalysis):
	_defaults = {
              	'jet_R': 0.4,
								'jet_algorithm': fj.antikt_algorithm,
								'part_eta_max': 0.9,
								'bg_y_max': 1.5,
								'bg_grid_spacing': 0.2,
								}
	def __init__(self, **kwargs):
		super(JetFinding, self).__init__(**kwargs)
		self.jet_eta_max = self.part_eta_max - self.jet_R
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.jet_eta_max + self.jet_R))
		self.jet_selector = fj.SelectorAbsEtaMax(self.jet_eta_max)
		self.bg_estimator = fj.GridMedianBackgroundEstimator(self.bg_y_max, self.bg_grid_spacing)

	def analyze(self, psjv):
		# estimate event background rho with grid estimator
		self.bg_estimator.set_particles(psjv)
		self.rho = self.bg_estimator.rho()
		self.sigma = self.bg_estimator.sigma()
		self.ca = fj.ClusterSequenceArea(psjv, self.jet_def, self.area_def)
		self.jets = fj.sorted_by_pt(self.jet_selector(self.ca.inclusive_jets()))


class LeadingJetAnalysis(BaseAnalysis):
	_defaults = {
								'name' : 'leadjets',
              	'jet_Rs': [0.2, 0.4, 0.6],
								'nleading_write': -1,
								'jet_algorithm': fj.antikt_algorithm,
								'part_eta_max': 0.9,
								'bg_y_max': 1.5,
								'bg_grid_spacing': 0.1,
								'write_constituents': False,
								}

	def __init__(self, **kwargs):
		super(LeadingJetAnalysis, self).__init__(**kwargs)
		# open an output rootfile and create a tree
		_rf = SingleRootFile()
		# print(_rf.root_file.GetName())
		# self.rt = RTreeWriter(name=self.name, file_name=f'{self.name}.root')
		self.rt_e = RTreeWriter(name=self.name+'_ev', fout=_rf.root_file)
		self.rt_j = RTreeWriter(name=self.name+'_js', fout=_rf.root_file)
		self.tne = ROOT.TNtuple(f'tne_{self.name}', f'tne_{self.name}', "mult:track_count:centr:jet_count:bgrho:bgsigma")
		# declare jet analysis for each jet radius
		self.jet_findings = []
		self.tnj = {}
		for R in self.jet_Rs:
			self.jet_findings.append(JetFinding(jet_R=R, jet_algorithm=self.jet_algorithm, part_eta_max=self.part_eta_max, bg_y_max=self.bg_y_max, bg_grid_spacing=self.bg_grid_spacing))
			self.tnj[R] = ROOT.TNtuple(f"tnj_{self.name}_R{R}".replace('.', ''), f"e_{self.name}_R{R}".replace('.', ''),  "mult:track_count:centr:jet_count:bgrho:bgsigma:pt:eta:phi")

	def add_constituents_to_write(self, j, nlead, jetR):
		for c in fj.sorted_by_pt(j.constituents()):
			idx = c.user_index()
			_cdict = self.cs_to_write.get(idx, {})
			# new constituent
			if _cdict == {}: 
				_cdict['c'] = c
				_cdict['ic'] = c.user_index()
				self.cs_to_write[idx] = _cdict
			_cdict[jetR] = nlead

	def analyze(self, e):
		data_ev = {}
		data_ev['cent'] = float(e.centrality)
		data_ev['ntracks'] = e.track_count
		data_ev['npsjv'] = len(e.psjv)
		data_ev['mult'] = float(e.multiplicity)
		data = {}
		data['cent'] = float(e.centrality)
		data['ntracks'] = e.track_count
		data['npsjv'] = len(e.psjv)
		data['mult'] = float(e.multiplicity)
		if self.write_constituents:
			self.cs_to_write = {}
		for jf in self.jet_findings:
			jf.analyze(e.psjv)
			data_ev['rho_R{}'.format(jf.jet_R).replace('.', '')] = jf.rho
			data_ev['sigma_R{}'.format(jf.jet_R).replace('.', '')] = jf.sigma
			if jf == self.jet_findings[0]:
				self.tne.Fill(e.multiplicity, e.track_count, e.centrality, len(jf.jets), jf.rho, jf.sigma)
			if len(jf.jets) > 0:
					sjetR = 'j_R{}'.format(jf.jet_R).replace('.', '')
					if self.nleading_write > 0:
						data['{}'.format(sjetR)] = [j for ij, j in enumerate(jf.jets) if ij < self.nleading_write]
						data['{}_n'.format(sjetR)] = [ij for ij, j in enumerate(jf.jets) if ij < self.nleading_write]
						_ = [self.tnj[jf.jet_R].Fill(e.multiplicity, e.track_count, e.centrality, len(jf.jets), jf.rho, jf.sigma, j.pt(), j.eta(), j.phi()) for ij, j in enumerate(jf.jets) if ij < self.nleading_write]
						if self.write_constituents:
							_ = [self.add_constituents_to_write(j, ij, jf.jet_R) for ij, j in enumerate(jf.jets) if ij < self.nleading_write]
					else:
						data['{}'.format(sjetR)] = jf.jets
						data['{}_n'.format(sjetR)] = [ij for ij, j in enumerate(jf.jets)]
						_ = [self.tnj[jf.jet_R].Fill(e.multiplicity, e.track_count, e.centrality, len(jf.jets), jf.rho, jf.sigma, j.pt(), j.eta(), j.phi()) for j in jf.jets]
						if self.write_constituents:
							_ = [self.add_constituents_to_write(j, ij, jf.jet_R) for ij, j in enumerate(jf.jets)]
					data['{}_rho'.format(sjetR)] = jf.rho
		# data['constituents'] = [c for c in e.psjv if c.user_index() in self.cs_indexes_to_write]
		# print('-'*20)
		# print(self.cs_to_write)
		if self.write_constituents:
			constit_data = {'c' : []}
			for jf in self.jet_findings:
				sjetR = '{}'.format(jf.jet_R).replace('.', '')
				constit_data[sjetR] = []
			for k, v in self.cs_to_write.items():
				v = self.cs_to_write[k]
				constit_data['c'].append(v['c'])
				for jf in self.jet_findings:
					sjetR = '{}'.format(jf.jet_R).replace('.', '')
					constit_data[sjetR].append(v.get(jf.jet_R, -1))
				# print(v)
			self.rt_j.fill_branches(b = data, c = constit_data)
			# print(constit_data)
		else:
			self.rt_j.fill_branches(b = data)
		self.rt_j.fill_tree()
		self.rt_e.fill_branches(e = data_ev)
		self.rt_e.fill_tree()

	def finalize(self):
		self.rt.write()
		# self.rt.write_and_close()
		# self.root_output.close()


def main():
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
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

	args = parser.parse_args()
	print(args)

	globals.tqdm_silent = args.no_tqdm

	# open the output file
	root_file = SingleRootFile(fname=args.output)
	print(root_file.root_file.GetName())
 
	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# get the fj banner out of the way
	fj.ClusterSequence().print_banner()

	cs = None
	if args.cs_dRmax > 0:
		cs = CEventSubtractor(alpha=args.cs_alpha, max_distance=args.cs_dRmax, max_eta=args.part_eta_max, bge_rho_grid_size=0.1, max_pt_correct=100)
		print(cs)
	parts_selector = fj.SelectorAbsEtaMax(args.part_eta_max)

	# open the output file
	an = LeadingJetAnalysis(name='ljetanastd', part_eta_max=args.part_eta_max, bg_y_max=1.5, bg_grid_spacing=0.1, save_tracks=args.save_tracks, write_constituents=args.save_tracks, nleading_write=args.nlead)
	if args.cs_dRmax > 0:
		an_cs = LeadingJetAnalysis(name='ljetanacs', part_eta_max=args.part_eta_max, bg_y_max=1.5, bg_grid_spacing=0.1, save_tracks=args.save_tracks, write_constituents=args.save_tracks, nleading_write=args.nlead)

	# event loop using the data source directly
	for i,e in enumerate(data_source.next_event()):
		if e.centrality < args.cent_min or e.centrality > args.cent_max:
			continue
		psjv = data_fj.data_tracks_to_pseudojets(e)
		e.psjv = parts_selector(psjv)
		an.analyze(e)
		if cs:
			psjv_cs = cs.process_event(psjv)
			e.psjv = parts_selector(psjv_cs)
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
