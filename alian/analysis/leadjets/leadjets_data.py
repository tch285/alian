#!/usr/bin/env python

import argparse
import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
from alian.io import data_io
from alian.utils import data_fj
from alian.steer.glob import globals
from alian.utils.treewriter import RTreeWriter

import ROOT


class BaseAnalysis(heppyy.GenericObject):
    _defaults = {}

    def __init__(self, **kwargs):
        super(BaseAnalysis, self).__init__(**kwargs)
        self.results = []
        for k, val in self.__class__._defaults.items():
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, val)


# this is finding jets - not doing anything with them
class JetFinding(BaseAnalysis):
	_defaults = {
              	'jet_R': 0.4,
								'jet_algorithm': fj.antikt_algorithm,
								'part_eta_max': 0.9,
								'bg_y_max': 1.5,
								'bg_grid_spacing': 0.1,
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
		self.rt = RTreeWriter(name=self.name, file_name=f'{self.name}.root')
		# declare jet analysis for each jet radius
		self.jet_findings = []
		for R in self.jet_Rs:
			self.jet_findings.append(JetFinding(jet_R=R, jet_algorithm=self.jet_algorithm, part_eta_max=self.part_eta_max, bg_y_max=self.bg_y_max, bg_grid_spacing=self.bg_grid_spacing))

	def add_constituents_to_write(self, jets):
		for j in jets:
			for c in j.constituents():
				print(c.user_index())
				self.cs_indexes_to_write.append(c.user_index())
		self.cs_indexes_to_write = list(set(self.cs_indexes_to_write))

	def analyze(self, e):
		data = {}
		data['cent'] = e.centrality
		data['ntracks'] = e.track_count
		if self.write_constituents:
			self.cs_to_write = std.vector[fj.PseudoJet]()
			self.cs_indexes_to_write = []
		for jf in self.jet_findings:
			jf.analyze(e.psjv)
			if len(jf.jets) > 0:
					sjetR = 'j_R{}'.format(jf.jet_R).replace('.', '')
					if self.nleading_write > 0:
						data['{}'.format(sjetR)] = [j for ij, j in enumerate(jf.jets) if ij < self.nleading_write]
						data['{}_n'.format(sjetR)] = [ij for ij, j in enumerate(jf.jets) if ij < self.nleading_write]
						if self.write_constituents:
							self.add_constituents_to_write(data['{}'.format(sjetR)])
					else:
						data['{}'.format(sjetR)] = jf.jets
						data['{}_n'.format(sjetR)] = [ij for ij, j in enumerate(jf.jets)]
						if self.write_constituents:
							self.add_constituents_to_write(data['{}'.format(sjetR)])
					data['{}_rho'.format(sjetR)] = jf.rho
		self.rt.fill_branches(b = data)
		self.rt.fill_tree()

	def finalize(self):
		self.rt.write_and_close()
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

	args = parser.parse_args()
	print(args)

	globals.tqdm_silent = args.no_tqdm

	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# get the fj banner out of the way
	fj.ClusterSequence().print_banner()

	# open the output file
	an = LeadingJetAnalysis(part_eta_max=0.9, bg_y_max=1.5, bg_grid_spacing=0.1, save_tracks=args.save_tracks, write_constituents=args.save_tracks, nleading_write=args.nlead)

	# event loop using the data source directly
	for i,e in enumerate(data_source.next_event()):
		if e.centrality < args.cent_min or e.centrality > args.cent_max:
			continue
		# psjv = data_fj.data_tracks_to_pseudojets(e)
		e.psjv = data_fj.data_tracks_to_pseudojets(e)
		an.analyze(e)

  # finalize the analyses - write the output, close the file, etc.
	an.finalize()


if __name__ == '__main__':
		main()
