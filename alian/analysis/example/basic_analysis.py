#!/usr/bin/env python

"""
This script is used on the ROOT file (provided as an argument) in
the BerkeleyTree format (YAML file with structure must be provided
with -t) to perform the analysis of the events containing isolated
photons and jets.
The script performs necessary cuts specified in YAML file (enabled
with -c), creates QA histograms for the cleaned data, and stores them
in the output ROOT file (specified with -o)
"""


import argparse
import heppyy
fj = heppyy.load_cppyy('fastjet')
from alian.io import data_io
from alian.utils import data_fj
from alian.analysis.base import AnalysisSelector, Event, Track, Cluster, BaseAnalysis, HistogramCollection, setup_logger
import numpy as np
import ROOT
from typing import List

#### DATA STRUCTURES ####

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
        self.ca = fj.ClusterSequenceArea(e.psjv, self.jet_def, self.area_def)
        self.jets = fj.sorted_by_pt(self.jet_selector(self.ca.inclusive_jets()))

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
        self.h_photon_jet_phi = ROOT.TH3F("Delta Phi separation between photons and jets", "Delta Phi separation between photons and jets", 50, 0, np.pi, 50, 0, 80, 200, 0, 80)
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
            self.h_cluster_phi_vs_energy[label] = ROOT.TH2F("cluster_phi_vs_energy_" + hist_label, "Cluster Phi (" + label + "); counts", 360, 0, 2*np.pi, self.photon_energy_nbins, self.min_photon_energy_range, self.max_photon_energy_range)
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
            self.h_track_phi_vs_pt[label] = ROOT.TH2F("track_phi_vs_pt_" + hist_label, "Track Phi (" + label + "); counts", 360, 0, 2*np.pi, self.track_pt_nbins, self.min_track_pt_range, self.max_track_pt_range)
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
            self.h_jet_phi_vs_pt[label] = ROOT.TH2F("jet_phi_" + hist_label, "Jet Phi (" + label + "); GeV; counts", 360, 0, 2*np.pi, self.jet_pt_nbins, self.min_jet_pt_range, self.max_jet_pt_range)
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
            logger.info("{}: {}".format(name, self.QAStatistics[name]), stacklevel = 2)

#### STATISTICS ####

#### CUTS ####

def continuousAngleDiff(theta_1, theta_2):
    """Calculate phi angle difference."""
    theta_diff = abs(theta_1 - theta_2)
    return theta_diff if theta_diff < np.pi else 2*np.pi - theta_diff

def shouldSelectEvent(event: Event) -> bool:
    ev_pass, evsel_pass, trgsel_pass = cuts.event.selects(event, verbose = True)
    retval = True
    if not evsel_pass:
        hists.fillQAStatistics("Number of Rejected Events due to event selection (event_selection)", 1)
        retval = False
    if not trgsel_pass:
        hists.fillQAStatistics("Number of Rejected Events due to event selection (triggersel)", 1)
        retval = False
    return retval

def shouldSelectClusterAsPhotonCandidate(cluster: Cluster, tracks: List[Track]) -> bool:
    # Cluster shape
    if not cuts.cluster.pass_shape(cluster):
        hists.fillQAStatistics("Number of clusters excluded due to shape", 1)
        return False
    hists.fillClusterQA(cluster, "After Cluster Shape Cut")
    # Cluster time
    if not cuts.cluster.pass_time(cluster):
        hists.fillQAStatistics("Number of clusters excluded due to time", 1)
        return False
    hists.fillClusterQA(cluster, "After Cluster Time Cut")
    # Checking the number of cells in the cluster
    if not cuts.cluster.pass_ncells(cluster):
        hists.fillQAStatistics("Number of clusters excluded due to number of cell", 1)
        return False
    hists.fillClusterQA(cluster, "After Cluster Cell Number Cut")
    # Cluster energy
    if not cuts.cluster.pass_E(cluster):
        hists.fillQAStatistics("Number of clusters excluded due to low energy", 1)
        return False
    hists.fillClusterQA(cluster, "After Cluster Energy Cut")
    # Is cluster exotic?
    if not cuts.cluster.pass_exotic(cluster):
        hists.fillQAStatistics("Number of clusters excluded due to being exotic", 1)
        return False
    hists.fillClusterQA(cluster, "After Cluster Exoticity Cut")
    # TODO: nlm
    # Exclude charged particles by finding the tracks leading to the cluster
    for track in tracks:
        if cuts.cluster.is_geo_matched(cluster, track):
            if cuts.cluster.is_e_pt_matched(cluster, track):
                hists.fillQAStatistics("Number of clusters identified as charged particles", 1)
                return False
    hists.fillClusterQA(cluster, "After Cluster Charged Particle Cut")
    return True

def findPhotonCandidates(clusters: List[Cluster], tracks: List[Track]) -> List[Cluster]:
    hists.fillClustersQA(clusters, "Raw Clusters")
    hists.fillQAStatistics("Total number of clusters candidates", len(clusters))
    candidates = [cluster for cluster in clusters if cuts.cluster.selects(cluster, tracks)]
    hists.fillQAStatistics("Number of potential photon candidates", len(candidates))
    hists.fillClustersQA(candidates, "Photon Candidates")
    return candidates

def shouldSelectTrack(track: Track) -> bool:
    track_pass, tracksel_pass, pt_min_pass = cuts.track.selects(track, verbose = True)
    if not tracksel_pass:
        hists.fillQAStatistics("Number of tracks excluded due to track selection bit", 1)
        return False
    hists.fillTrackQA(track, "After Track Selection Cut")
    if not pt_min_pass:
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
    is_iso_cluster, total_pt = cuts.cluster.is_iso(candidate, tracks, verbose = True)
    # Fill corresponding histogram
    hists.h_total_track_pt_in_isolation_cone.Fill(total_pt)
    if not is_iso_cluster:
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
    logger = setup_logger(__name__)

    parser = argparse.ArgumentParser(description="Run analysis on ROOT file using YAML configuration.")
    parser.add_argument('-i', '--input-file', type=str, help="Input file or file containing a list of input files.", required = True)
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
    cuts = AnalysisSelector.from_file(args.cuts)
    # hc = HistogramCollection(args.cuts)

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
    logger.info(f"Saved histograms to {args.output}")
