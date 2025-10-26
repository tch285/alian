# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python

import argparse
import heppyy
fj = heppyy.load_cppyy('fastjet')
from alian.io import data_io
from alian.utils import data_fj
import pandas as pd
import numpy as np
import ROOT
import math
from typing import List
import yaml

#### DATA STRUCTURES ####

class Cluster:
	def __str__(self):
		return f"Cluster(energy={self.energy}, eta={self.eta}, phi={self.phi}, m02={self.m02}, m20={self.m20}, ncells={self.ncells})"

	def __init__(self):
		self.energy = 0
		self.eta = 0
		self.phi = 0
		self.m02 = 0
		self.m20 = 0
		self.ncells = 0
		self.time = 0
		self.isexotic = 0
		self.distancebadchannel = 0
		self.nlm = 0
		self.clusterdef = 0
		self.matchedTrackIndex = 0

class Track:
	def __str__(self):
		return f"Track(pT={self.pt}, eta={self.eta}, phi={self.phi})"

	def __init__(self):
		self.pt = 0
		self.eta = 0
		self.phi = 0
		self.label = 0
		self.tracksel = 0

	def total_p(self):
		return self.pt * math.cosh(self.eta)

class Event:
	def __init__(self, eventStruct):
		self.run_number = eventStruct.data['run_number']
		self.event_selection = eventStruct.data['event_selection']
		self.triggersel = eventStruct.data['triggersel']
		self.centrality = eventStruct.data['centrality']
		self.multiplicity = eventStruct.data['multiplicity']
		self.tracks = []
		for i in range(len(eventStruct.data['track_data_pt'])):
			track = Track()
			track.pt = eventStruct.data['track_data_pt'][i]
			track.eta = eventStruct.data['track_data_eta'][i]
			track.phi = eventStruct.data['track_data_phi'][i]
			track.label = eventStruct.data['track_data_label'][i]
			track.tracksel = eventStruct.data['track_data_tracksel'][i]
			self.tracks.append(track)
		self.clusters = []
		for i in range(len(eventStruct.data['cluster_data_energy'])):
			cluster = Cluster()
			cluster.energy = eventStruct.data['cluster_data_energy'][i]
			cluster.eta = eventStruct.data['cluster_data_eta'][i]
			cluster.phi = eventStruct.data['cluster_data_phi'][i]
			cluster.m02 = eventStruct.data['cluster_data_m02'][i]
			cluster.m20 = eventStruct.data['cluster_data_m20'][i]
			cluster.ncells = eventStruct.data['cluster_data_ncells'][i]
			cluster.time = eventStruct.data['cluster_data_time'][i]
			cluster.isexotic = eventStruct.data['cluster_data_isexotic'][i]
			cluster.distancebadchannel = eventStruct.data['cluster_data_distancebadchannel'][i]
			cluster.nlm = eventStruct.data['cluster_data_nlm'][i]
			cluster.clusterdef = eventStruct.data['cluster_data_clusterdef'][i]
			cluster.matchedTrackIndex = eventStruct.data['cluster_data_matchedTrackIndex'][i]
			self.clusters.append(cluster)

#### DATA STRUCTURES ####

#### JET ANALYSIS ####

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
		super(JetFinding, self).__init__(**kwargs)
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

#### JET ANALYSIS ####

#### ENUMS ####

class Enums:
	EventTriggerSel = {
		"noTrigSel": 0,
		"JetChLowPt": 1,
		"JetChHighPt": 2,
		"TrackLowPt": 3,
		"TrackHighPt": 4,
		"JetD0ChLowPt": 5,
		"JetD0ChHighPt": 6,
		"JetLcChLowPt": 7,
		"JetLcChHighPt": 8,
		"EMCALReadout": 9,
		"JetFullHighPt": 10,
		"JetFullLowPt": 11,
		"JetNeutralHighPt": 12,
		"JetNeutralLowPt": 13,
		"GammaVeryHighPtEMCAL": 14,
		"GammaVeryHighPtDCAL": 15,
		"GammaHighPtEMCAL": 16,
		"GammaHighPtDCAL": 17,
		"GammaLowPtEMCAL": 18,
		"GammaLowPtDCAL": 19,
		"GammaVeryLowPtEMCAL": 20,
		"GammaVeryLowPtDCAL": 21
	}

	EventCollisionSel = {
		"sel8": 0,
		"sel8Full": 1,
		"sel7": 2,
		"selMC": 3,
		"selUnanchoredMC": 4,
		"sel7KINT7": 5
	}

	TrackSel = {
		"trackSign": 0,
		"globalTrack": 1,
		"qualityTrack": 2,
		"hybridTrack": 3,
		"uniformTrack": 4,
		"uniformTrackWoDCA": 5
	}

#### ENUMS ####

#### STATISTICS ####

class Histograms:
	def __init__(self):
		self.QAStatisticsHists = ROOT.TH1F("QA Statistics Histograms", "Counts", 1, 0, 1)
		self.QAStatistics = {}
		# Histogram Limits
		self.min_photon_energy_range = 0
		self.max_photon_energy_range = 50
		self.photon_energy_nbins = 200
		self.min_track_pt_range = 0
		self.max_track_pt_range = 50
		self.track_pt_nbins = 200
		self.min_jet_pt_range = 0
		self.max_jet_pt_range = 50
		self.jet_pt_nbins = 200
		# Photo QA Histograms
		self.h_cluster_energy = {}
		self.h_cluster_phi_vs_energy = {}
		self.h_cluster_eta_vs_energy = {}
		self.h_cluster_m02_vs_energy = {}
		self.h_cluster_m20_vs_energy = {}
		self.h_cluster_ncells_vs_energy = {}
		self.h_cluster_time_vs_energy = {}
		self.h_cluster_isexotic_vs_energy = {}
		# TODO
		# self.h_cluster_nlm_vs_energy = {}
		# Track QA Histograms
		self.h_track_pt = {}
		self.h_track_phi_vs_pt = {}
		self.h_track_eta_vs_pt = {}
		# Jet QA Histograms
		self.h_jet_energy = {}
		self.h_jet_energy_t = {}
		self.h_jet_pt = {}
		self.h_jet_phi_vs_pt = {}
		self.h_jet_eta_vs_pt = {}
		self.h_jet_rapidity_vs_pt = {}
		# Other QA Histograms
		self.h_total_track_pt_in_isolation_cone = ROOT.TH1F("Total photon candidate pT within isolation cone", "Photon pT; GeV; counts", self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
		self.h_photon_jet_phi = ROOT.TH3F("Delta Phi separation between photons and jets", "Delta Phi separation between photons and jets", 50, 0, math.pi, 50, 0, 80, 200, 0, 80)
		self.h_photon_jet_pt = ROOT.TH2F("Jet-Photon pair pT", "Jet-Photon pair pT", 100, 0, 80, 100, 0, 80)

	def fillQAStatistics(self, name: str, count: int):
		self.QAStatisticsHists.Fill(name, count)
		if name not in self.QAStatistics:
			self.QAStatistics[name] = 0
		self.QAStatistics[name] += count

	def fillClusterQA(self, cluster: Cluster, label: str):
		self.fillClustersQA([cluster], label)

	def fillClustersQA(self, clusters: List[Cluster], label: str):
		hist_label = label.lower().replace(" ", "_")
		root_directory_name = "Cluster Hists (" + label + ")"
		tfile.mkdir(root_directory_name, "", True)
		tfile.cd(root_directory_name)

		if label not in self.h_cluster_energy:
			self.h_cluster_energy[label] = ROOT.TH1F("cluster_energy_" + hist_label, "Cluster Energy (" + label + "); GeV; counts", self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_phi_vs_energy[label] = ROOT.TH2F("cluster_phi_vs_energy_" + hist_label, "Cluster Phi (" + label + "); counts", 360, 0, 2*math.pi, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_eta_vs_energy[label] = ROOT.TH2F("cluster_eta_vs_energy_" + hist_label, "Cluster Eta (" + label + "); counts", 100, -1.5, 1.5, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_m02_vs_energy[label] = ROOT.TH2F("cluster_m02_vs_energy_" + hist_label, "Cluster M02 (" + label + "); counts", 100, 0, 3, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_m20_vs_energy[label] = ROOT.TH2F("cluster_m20_vs_energy_" + hist_label, "Cluster M20 (" + label + "); counts", 100, 0, 1, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_ncells_vs_energy[label] = ROOT.TH2F("cluster_ncells_vs_energy_" + hist_label, "Cluster Number of Cells (" + label + "); counts", 10, 0, 10, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_time_vs_energy[label] = ROOT.TH2F("cluster_time_vs_energy_" + hist_label, "Cluster Time (" + label + "); ns", 100, -40, 40, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			self.h_cluster_isexotic_vs_energy[label] = ROOT.TH2F("cluster_isexotic_vs_energy_" + hist_label, "Cluster Number of Cells (" + label + "); counts", 2, 0, 1, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
			# TODO
			# self.h_cluster_nlm_vs_energy[label] = ROOT.TH2F("cluster_nlm_vs_energy_" + hist_label, "Cluster NLM (" + label + "); counts", 100, 0, 1, energy_nbins, 0, max_energy_range)
		for cluster in clusters:
			self.h_cluster_energy[label].Fill(cluster.energy)
			self.h_cluster_phi_vs_energy[label].Fill(cluster.phi, cluster.energy)
			self.h_cluster_eta_vs_energy[label].Fill(cluster.eta, cluster.energy)
			self.h_cluster_m02_vs_energy[label].Fill(cluster.m02, cluster.energy)
			self.h_cluster_m20_vs_energy[label].Fill(cluster.m20, cluster.energy)
			self.h_cluster_ncells_vs_energy[label].Fill(cluster.ncells, cluster.energy)
			self.h_cluster_time_vs_energy[label].Fill(cluster.time, cluster.energy)
			self.h_cluster_isexotic_vs_energy[label].Fill(cluster.isexotic, cluster.energy)
			# TODO
			# self.h_cluster_nlm_vs_energy[label].Fill(cluster.m20, cluster.energy)
		tfile.cd()


	def fillTrackQA(self, track: Track, label: str):
		self.fillTracksQA([track], label)

	def fillTracksQA(self, tracks: List[Track], label: str):
		hist_label = label.lower().replace(" ", "_")
		root_directory_name = "Track Hists (" + label + ")"
		tfile.mkdir(root_directory_name, "", True)
		tfile.cd(root_directory_name)

		if label not in self.h_track_pt:
			self.h_track_pt[label] = ROOT.TH1F("track_pt_" + hist_label, "Track pT (" + label + "); GeV; counts", self.track_pt_nbins, self.min_track_pt_range, self.max_track_pt_range)
			self.h_track_phi_vs_pt[label] = ROOT.TH2F("track_phi_vs_pt_" + hist_label, "Track Phi (" + label + "); counts", 360, 0, 2*math.pi, self.track_pt_nbins, self.min_track_pt_range, self.max_track_pt_range)
			self.h_track_eta_vs_pt[label] = ROOT.TH2F("track_eta_vs_pt_" + hist_label, "Track Eta (" + label + "); counts", 100, -1.5, 1.5, self.track_pt_nbins, self.min_track_pt_range, self.max_track_pt_range)
		for track in tracks:
			self.h_track_pt[label].Fill(track.pt)
			self.h_track_phi_vs_pt[label].Fill(track.phi, track.pt)
			self.h_track_eta_vs_pt[label].Fill(track.eta, track.pt)
		tfile.cd()

	def fillJetsQA(self, jets: List[fj.PseudoJet], label: str):
		hist_label = label.lower().replace(" ", "_")
		root_directory_name = "Jet Hists (" + label + ")"
		tfile.mkdir(root_directory_name, "", True)
		tfile.cd(root_directory_name)

		if label not in self.h_jet_energy:
			self.h_jet_energy[label] = ROOT.TH1F("jet_energy_" + hist_label, "Jet Energy (" + label + "); GeV; counts", 200, 0, 50)
			self.h_jet_energy_t[label] = ROOT.TH1F("jet_energy_t_" + hist_label, "Jet Transverse Energy (" + label + "); GeV; counts", 200, 0, 50)
			self.h_jet_pt[label] = ROOT.TH1F("jet_pt_" + hist_label, "Jet Transverse Momentum (" + label + "); GeV; counts", self.jet_pt_nbins, self.min_jet_pt_range, self.max_jet_pt_range)
			self.h_jet_phi_vs_pt[label] = ROOT.TH2F("jet_phi_" + hist_label, "Jet Phi (" + label + "); GeV; counts", 360, 0, 2*math.pi, self.jet_pt_nbins, self.min_jet_pt_range, self.max_jet_pt_range)
			self.h_jet_eta_vs_pt[label] = ROOT.TH2F("jet_eta_" + hist_label, "Jet Eta (" + label + "); GeV; counts", 100, -1.5, 1.5, self.jet_pt_nbins, self.min_jet_pt_range, self.max_jet_pt_range)
			self.h_jet_rapidity_vs_pt[label] = ROOT.TH2F("jet_rapidity_" + hist_label, "Jet Rapidity (" + label + "); counts", 200, 0, 10, self.jet_pt_nbins, self.min_jet_pt_range, self.max_jet_pt_range)
		for jet in jets:
			self.h_jet_energy[label].Fill(jet.E())
			self.h_jet_energy_t[label].Fill(jet.Et())
			self.h_jet_pt[label].Fill(jet.perp())
			self.h_jet_phi_vs_pt[label].Fill(jet.phi(), jet.perp())
			self.h_jet_eta_vs_pt[label].Fill(jet.eta(), jet.perp())
			self.h_jet_rapidity_vs_pt[label].Fill(jet.rapidity(), jet.perp())
		tfile.cd()

	def dump(self):
		for name in self.QAStatistics:
			print("{}: {}".format(name, self.QAStatistics[name]))

#### STATISTICS ####

#### CUTS ####

class Cuts:
	class EventCuts:
		def __init__(self):
			self.event_selection = [ Enums.EventCollisionSel["sel8"] ]
			self.triggersel = [ Enums.EventTriggerSel["GammaHighPtEMCAL"] ]

		def setCuts(self, cuts: dict):
			if "event_selection" in cuts:
				self.event_selection = []
				for elem in cuts["event_selection"].split(","):
					self.event_selection.append(Enums.EventCollisionSel[elem])
			if "triggersel" in cuts:
				self.triggersel = []
				for elem in cuts["triggersel"].split(","):
					self.triggersel.append(Enums.EventTriggerSel[elem])

		def dump(self):
			print("Event cuts:")
			print("\tevent_selection:", self.event_selection)
			print("\ttriggersel:", self.triggersel)

	class ClusterCuts:
		def __init__(self):
			self.min_E = 0.7 # GeV
			self.min_m02 = 0.1
			self.max_m02 = 0.3
			self.min_time = -30 # ns
			self.max_time = 35 # ns
			self.min_ncells = 2
			self.delta_phi = 0.05
			self.delta_eta = 0.05
			self.energy_pt_ratio_threshold = 1.75
			self.isolation_cone_radius = 0.4
			self.isolation_pt_threshold = 1.5 # GeV

		def setCuts(self, cuts: dict):
			if "min_E" in cuts:
				self.min_E = cuts["min_E"]
			if "min_m02" in cuts:
				self.min_m02 = cuts["min_m02"]
			if "max_m02" in cuts:
				self.max_m02 = cuts["max_m02"]
			if "min_time" in cuts:
				self.min_time = cuts["min_time"]
			if "max_time" in cuts:
				self.max_time = cuts["max_time"]
			if "min_ncells" in cuts:
				self.min_ncells = cuts["min_ncells"]
			if "delta_phi" in cuts:
				self.delta_phi = cuts["delta_phi"]
			if "delta_eta" in cuts:
				self.delta_eta = cuts["delta_eta"]
			if "energy_pt_ratio_threshold" in cuts:
				self.energy_pt_ratio_threshold = cuts["energy_pt_ratio_threshold"]
			if "isolation_cone_radius" in cuts:
				self.isolation_cone_radius = cuts["isolation_cone_radius"]
			if "isolation_pt_threshold" in cuts:
				self.isolation_pt_threshold = cuts["isolation_pt_threshold"]

		def dump(self):
			print("Cluster cuts:")
			print("\tmin_E:", self.min_E, "GeV")
			print("\tmin_m02:", self.min_m02)
			print("\tmax_m02:", self.max_m02)
			print("\tmin_time:", self.min_time)
			print("\tmax_time:", self.max_time)
			print("\tmin_ncells:", self.min_ncells)
			print("\tdelta_phi:", self.delta_phi)
			print("\tdelta_eta:", self.delta_eta)
			print("\tenergy_pt_ratio_threshold:", self.energy_pt_ratio_threshold)
			print("\tisolation_cone_radius:", self.isolation_cone_radius)
			print("\tisolation_pt_threshold:", self.isolation_pt_threshold)

	class TrackCuts:
		def __init__(self):
			self.min_pt = 0.150
			self.tracksel = Enums.TrackSel["globalTrack"]

		def setCuts(self, cuts: dict):
			if "min_pt" in cuts:
				self.min_pt = cuts["min_pt"]
			if "tracksel" in cuts:
				self.tracksel = []
				for elem in cuts["tracksel"].split(","):
					self.tracksel.append(Enums.TrackSel[elem])

		def dump(self):
			print("Track cuts:")
			print("\ttracksel:", self.tracksel)
			print("\tmin_pt:", self.min_pt, "GeV")

	def __init__(self):
		self.event = Cuts.EventCuts()
		self.track = Cuts.TrackCuts()
		self.cluster = Cuts.ClusterCuts()

	def extractFromFile(self, file: str):
		with open(file) as stream:
			try:
				yaml_struct = yaml.safe_load(stream)
				if "event_cuts" in yaml_struct:
					self.event.setCuts(yaml_struct["event_cuts"])
				if "track_cuts" in yaml_struct:
					self.track.setCuts(yaml_struct["track_cuts"])
				if "cluster_cuts" in yaml_struct:
					self.cluster.setCuts(yaml_struct["cluster_cuts"])
			except yaml.YAMLError as exc:
				print("Error processing {} file".format(file))
				exit(0)

	def dump(self):
		self.event.dump()
		self.track.dump()
		self.cluster.dump()

#### CUTS ####

def continuousAngleDiff(theta_1, theta_2):
  theta_diff = abs(theta_1 - theta_2)
  return theta_diff if theta_diff < math.pi else 2*math.pi - theta_diff

def isBitSet(mask: int, bit: int):
	return mask & (1 << bit)

def isAnyBitSet(mask: int, bits: List[int]):
	for bit in bits:
		if isBitSet(mask, bit):
			return True
	return False

def shouldSelectEvent(event: Event) -> bool:
	if not isAnyBitSet(event.event_selection, cuts.event.event_selection):
		hists.fillQAStatistics("Number of Rejected Events due to event selection (event_selection)", 1)
		return False
	if not isAnyBitSet(event.triggersel, cuts.event.triggersel):
		hists.fillQAStatistics("Number of Rejected Events due to event selection (triggersel)", 1)
		return False
	return True

def doTrackClusterMatch(cluster: Cluster, track: Track) -> bool:
	return	abs(cluster.eta - track.eta) <= cuts.cluster.delta_eta and \
			abs(cluster.phi - track.phi) <= cuts.cluster.delta_phi

def shouldSelectClusterAsPhotonCandidate(cluster: Cluster, tracks: List[Track]) -> bool:
	# Cluster shape
	if not (cluster.m02 > cuts.cluster.min_m02 and cluster.m02 < cuts.cluster.max_m02):
		hists.fillQAStatistics("Number of clusters exluded due to shape", 1)
		return False
	hists.fillClusterQA(cluster, "After Cluster Shape Cut")
	# Cluster time
	if not (cluster.time > cuts.cluster.min_time and cluster.time < cuts.cluster.max_time):
		hists.fillQAStatistics("Number of clusters exluded due to time", 1)
		return False
	hists.fillClusterQA(cluster, "After Cluster Time Cut")
	# Checking the number of cells in the cluster
	if cluster.ncells < cuts.cluster.min_ncells:
		hists.fillQAStatistics("Number of clusters exluded due to number of cell", 1)
		return False
	hists.fillClusterQA(cluster, "After Cluster Cell Number Cut")
	# Cluster energy
	if cluster.energy < cuts.cluster.min_E:
		hists.fillQAStatistics("Number of clusters exluded due to low energy", 1)
		return False
	hists.fillClusterQA(cluster, "After Cluster Energy Cut")
	# Is cluster exotic?
	if cluster.isexotic:
		hists.fillQAStatistics("Number of clusters exluded due to being exotic", 1)
		return False
	hists.fillClusterQA(cluster, "After Cluster Exoticity Cut")
	# TODO: nlm
	# Exclude charged particles by finding the tracks leading to the cluster
	for track in tracks:
		if doTrackClusterMatch(cluster, track):
			if cluster.energy / track.total_p() <= cuts.cluster.energy_pt_ratio_threshold:
				hists.fillQAStatistics("Number of clusters identified as charged particles", 1)
				return False
	hists.fillClusterQA(cluster, "After Cluster Charged Particle Cut")
	return True

def findPhotonCandidates(clusters: List[Cluster], tracks: List[Track]) -> List[Cluster]:
	hists.fillClustersQA(clusters, "Raw Clusters")
	hists.fillQAStatistics("Total number of clusters candidates", len(clusters))
	candidates = [cluster for cluster in clusters if shouldSelectClusterAsPhotonCandidate(cluster, tracks)]
	hists.fillQAStatistics("Number of potential photon candidates", len(candidates))
	hists.fillClustersQA(candidates, "Photon Candidates")
	return candidates

def shouldSelectTrack(track: Track) -> bool:
	if not isAnyBitSet(track.tracksel, cuts.track.tracksel):
		hists.fillQAStatistics("Number of tracks excluded due to track selection bit", 1)
		return False
	hists.fillTrackQA(track, "After Track Selection Cut")
	if track.pt < cuts.track.min_pt:
		hists.fillQAStatistics("Number of tracks excluded due to pT cutoff", 1)
		return False
	hists.fillTrackQA(track, "After Track pT Cut")
	return True

def findTracks(tracks: List[Track]) -> List[Track]:
	hists.fillTracksQA(tracks, "Raw Tracks")
	hists.fillQAStatistics("Total number of track candidates", len(tracks))
	candidates = [track for track in tracks if shouldSelectTrack(track)]
	hists.fillQAStatistics("Number of selected tracks", len(candidates))
	hists.fillTracksQA(candidates, "Selected Tracks")
	return candidates

def isIsolatedPhotonCandidate(candidate: Cluster, tracks: List[Track]) -> bool:
	total_pt = 0
	for track in tracks:
		delta_phi = continuousAngleDiff(candidate.phi, track.phi)
		if np.sqrt((candidate.eta - track.eta) ** 2 + (delta_phi ** 2)) < cuts.cluster.isolation_cone_radius:
			total_pt += track.pt
	# Fill corresponding histogram
	hists.h_total_track_pt_in_isolation_cone.Fill(total_pt)
	if total_pt >= cuts.cluster.isolation_pt_threshold:
		return False
	return True

def excludeNonIsolatedPhotons(candidates: List[Cluster], tracks: List[Track]) -> int:
	candidates = [candidate for candidate in candidates if isIsolatedPhotonCandidate(candidate, tracks)]
	hists.fillQAStatistics("Number of isolated photon candidates", len(candidates))
	hists.fillClustersQA(candidates, "Isolated Photon Candidates")
	return candidates

def findJets(e) -> List[fj.PseudoJet]:
	e.psjv = data_fj.data_tracks_to_pseudojets(e)
	finder = JetFinding()
	finder.analyze(e)
	return finder.jets

def computePhotonJetCorrelations(photon_candidates: List[Cluster], jets: List[fj.PseudoJet]):
	for photon in photon_candidates:
		for jet in jets:
			# x: jet pT (100 bins, from 0 to 80)
			# y: photon pT (100 bins, from 0 to 80)
			hists.h_photon_jet_pt.Fill(jet.perp(), photon.energy)
			delta_phi = continuousAngleDiff(photon.phi, jet.phi())
			# x: delta phi (50 bins, from -pi to pi)
			# y: jet pT (50 bins, from 0 to 80)
			# z: photon pT (200 bins, from 0 to 80)
			hists.h_photon_jet_phi.Fill(delta_phi, jet.perp(), photon.energy)
	return

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
	parser.add_argument('input_file', type=str, help="Input file or file containing a list of input files.")
	parser.add_argument("-e", "--entries", type=int, help="Number of entries to process.", default=-1) #-1 means all entries
	parser.add_argument("-o", "--output", type=str, help="File to write the analysis.", default="analysis.root")
	parser.add_argument('-t', '--tree-struct', type=str, help="Path to a YAML file describing the tree structure.", default=None)
	parser.add_argument('-c', '--cuts', type=str, help="Path to a YAML file describing the cuts to be applied to the data.", default=None)
	parser.add_argument('--lhc-run', type=int, help='LHC Run', default=3)

	args = parser.parse_args()

	# initialize an output root file
	tfile = ROOT.TFile(args.output, "recreate")
	tfile.cd()  # making sure we are in the right 'directory' to create some histograms

	# initialize histogram collection
	hists = Histograms()

	# initialize cuts from the YAML file
	cuts = Cuts()
	cuts.extractFromFile(args.cuts)

	# Print cuts
	cuts.dump()

	# initialize the data input
	data_source = data_io.DataInput(args.input_file, lhc_run=args.lhc_run, yaml_file=args.tree_struct, n_events=args.entries)

	# get the fj banner out of the way
	fj.ClusterSequence().print_banner()

	# TODO: Remove idx, needed only for test runs
	idx = 0
	# event loop using the data source directly
	for e in data_source.next_event():
		idx += 1
		# if idx > 100000:
		# 	break
		event = Event(e)
		hists.fillQAStatistics("Total number of events", 1)
		# Remove events based on the triggers
		if not shouldSelectEvent(event):
			hists.fillQAStatistics("Total number of rejected events", 1)
			continue
		# Find tracks satisfying the cuts
		tracks = findTracks(event.tracks)
		# Identify potential photon candidates from the clusters
		photon_candidates = findPhotonCandidates(event.clusters, tracks)
		# Exclude candidates that are not isolated
		isolated_photon_candidates = excludeNonIsolatedPhotons(photon_candidates, tracks)
		# Skip the event if there are no isolated photon candidates
		if len(isolated_photon_candidates) == 0:
			hists.fillQAStatistics("Number of Events with no isolated photons", 1)
			continue
		# Perform jet finding
		jets = findJets(e)
		hists.fillQAStatistics("Number of Jets", len(jets))
		# Calculate EECs
		# TODO
		hists.fillJetsQA(jets, "All Jets")
		# Look at correlations between jets and photon candidates
		computePhotonJetCorrelations(isolated_photon_candidates, jets)
		# TODO: projections needs to be done in a different program


	hists.dump()
	tfile.Write()
	# close the file
	tfile.Close()
