#!/usr/bin/env python

import argparse
import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
from alian.io import data_io
from alian.utils import data_fj
from alian.io.root_io import SingleRootFile

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
								'jet_eta_max': 0.5,
								'bg_y_max': 1.5,
								'bg_grid_spacing': 0.1,
								}
	def __init__(self, **kwargs):
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.jet_eta_max + self.jet_R))
		self.jet_selector = fj.SelectorAbsEtaMax(self.jet_eta_max)
		self.bg_estimator = fj.GridMedianBackgroundEstimator(self.bg_y_max, self.bg_grid_spacing)

	def analyze(self, e):
		# estimate event background rho with grid estimator
		self.bg_estimator.set_particles(e.psjv)
		self.rho = self.bg_estimator.rho()
		self.sigma = self.bg_estimator.sigma()
		# print('rho:', rho, 'centrality:', self.centrality, 'track_count:', self.track_count, 'e.psjv.size():', e.psjv.size())
		self.ca = fj.ClusterSequenceArea(e.psjv, self.jet_def, self.area_def)
		self.jets = fj.sorted_by_pt(self.jet_selector(self.ca.inclusive_jets()))


# this is finding jets and doing something with them - storing them in a TNtuple
class SimpleJetAnalysis(BaseAnalysis):
	_defaults = { 
								'name' : 'sja',
              	'jet_R': 0.4, 
								'jet_algorithm': fj.antikt_algorithm, 
								'jet_eta_max': 0.5,
								'bg_y_max': 1.5,
								'bg_grid_spacing': 0.1,
								'n_accepted_jets': 0,
								'n_accepted_events': 0
								}
	def __init__(self, **kwargs):
		super(SimpleJetAnalysis, self).__init__(**kwargs)
		self.root_output = SingleRootFile()
		self.tn_mult = ROOT.TNtuple(f"e_{self.name}", f"e_{self.name}",  "mult:track_count:centr:jet_count:bgrho:bgsigma")
		self.tn_jets = ROOT.TNtuple(f"j_{self.name}", f"j_{self.name}", "emult:track_count:centr:pt:eta:phi:m:e:jmult:nlead:leadpt:area:rho")
		if self.save_tracks:
			self.tn_tracks = ROOT.TNtuple(f't_{self.name}', f't_{self.name}', 'pt:eta:phi:m:jpt:jeta:jphi:jarea:jrho:jmult:jleadpt')
		self.root_output.add(self.tn_mult)
		self.root_output.add(self.tn_jets)
		self.results.append(self.root_output)

		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.jet_eta_max + self.jet_R))
		self.jet_selector = fj.SelectorAbsEtaMax(self.jet_eta_max)
		self.bg_estimator = fj.GridMedianBackgroundEstimator(self.bg_y_max, self.bg_grid_spacing)

	def analyze(self, e):
		self.n_accepted_events += 1
		# estimate event background rho with grid estimator
		self.bg_estimator.set_particles(e.psjv)
		self.rho = self.bg_estimator.rho()
		self.sigma = self.bg_estimator.sigma()
		# print('rho:', rho, 'centrality:', self.centrality, 'track_count:', self.track_count, 'e.psjv.size():', e.psjv.size())
		self.ca = fj.ClusterSequenceArea(e.psjv, self.jet_def, self.area_def)
		self.jets = fj.sorted_by_pt(self.jet_selector(self.ca.inclusive_jets()))
		jet_count = 0
		for ij, j in enumerate(self.jets):
				jconstits = [jconst for jconst in j.constituents() if jconst.perp() > 0.1]
				if len(jconstits) == 0:
						continue
				leadpt = fj.sorted_by_pt(jconstits)[0].perp()
				njconstits = len(jconstits)
				self.tn_jets.Fill(e.multiplicity, e.track_count, e.centrality, j.perp(), j.eta(), j.phi(), j.m(), j.e(), njconstits, ij, leadpt, j.area(), self.rho)
				jet_count += 1
				self.n_accepted_jets += 1
				if self.save_tracks:
					_ = [self.tn_tracks.Fill(t.perp(), t.eta(), t.phi(), t.m(), j.perp(), j.eta(), j.phi(), j.area(), self.rho, njconstits, leadpt) for t in jconstits]
		self.tn_mult.Fill(e.multiplicity, e.track_count, e.centrality, jet_count, self.rho, self.sigma)

	def finalize(self):
		print('[i] Accepted events:', self.n_accepted_events, 'Accepted jets:', self.n_accepted_jets)


def main():
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument("-o", "--output", type=str, help="Output file name.", default="analysis_results.root")
	parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=data_io.get_default_tree_structure())
	parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)
	parser.add_argument('--save-tracks', action='store_true', help='Save track information in the output file')
	parser.add_argument('--cent-min', type=int, help='Minimum centrality', default=-1)
	parser.add_argument('--cent-max', type=int, help='Maximum centrality', default=101)
 
	args = parser.parse_args()
	print(args)

	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# get the fj banner out of the way
	fj.ClusterSequence().print_banner()

	# open the output file
	root_output = SingleRootFile()
	analyses = []
	for R in [0.2, 0.4, 0.6]:
		analyses.append(SimpleJetAnalysis(name=f'sja_R{R}', jet_R=R, svae_tracks=args.save_tracks))

	# event loop using the data source directly
	for i,e in enumerate(data_source.next_event()):
		if e.centrality < args.cent_min or e.centrality > args.cent_max:
			continue
		psjv = data_fj.data_tracks_to_pseudojets(e)
		e.psjv = data_fj.data_tracks_to_pseudojets(e)
		for a in analyses:
			a.analyze(e)
  
  # finalize the analyses - write the output, close the file, etc.
	for a in analyses:
		a.finalize()
	root_output.close()
  
if __name__ == '__main__':
		main()