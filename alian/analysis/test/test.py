#!/usr/bin/env python


# import heppyy
# fj = heppyy.load_cppyy('fastjet')

# import alian.analysis.base
# import heppyy.util.fastjet_cppyy # noqa: F401
# # import heppyy.util.pythia8_cppyy
# # import heppyy.util.heppyy_cppyy

# from cppyy.gbl import fastjet as fj

# print(fj.antikt_algorithm)


import cppyy
# import sys
import numpy as np

# Load your C++ code
# cppyy.add_include_path(np.get_include())
# cppyy.load_library("libalian_fjutil.so")
# cppyy.include("fjutil/fjutil.hh")
# cppyy.include("fastjet/PseudoJet.hh")
# cppyy.include("fastjet/ClusterSequence.hh")
# from cppyy.gbl import alian

# from alian.analysis.base.selection import TrackSel
# from cppyy.gbl import fastjet as fj
# from cppyy.gbl import std
import heppyy
import heppyy.util.fastjet_cppyy  # noqa: F401
alian = heppyy.load_cppyy("alian")
from cppyy.gbl import fastjet as fj






# E  = np.array([9.0], dtype = np.float32)
# pt = np.array([3.0, 4.0], np.float32)
# eta = np.array([2.0, 1.0], dtype = np.float32)
# phi = np.array([3.0, 0.2], dtype = np.float32)

# m02 = np.array([0.3], dtype = np.float32)
# m20 = np.array([0.9], dtype = np.float32)
# ncells = np.array([0.9], dtype = np.uint16)
# time = np.array([0.9], dtype = np.float32)
# isexotic = np.array([0.9], dtype = np.bool)
# distancebadchannel = np.array([0.9], dtype = np.uint16)
# nlm = np.array([0.9], dtype = np.uint16)
# clusterdef = np.array([0.9], dtype = np.uint16)
# matchedTrackIndex = np.array([0.9], dtype = np.uint16)
# fj.ClusterSequence.print_banner()
# px = np.array([1.0, 3.0], dtype = np.float32)
# py = np.array([2.0, 4.0], dtype = np.float32)
# pz = np.array([3.0, 2.0], dtype = np.float32)
# tracksel = np.array([6, 9], dtype=np.uint8)
# arr1 = np.array([], dtype = np.float64)
# print(tracksel.dtype)
# tracks = alian.numpy_pxpypz_to_tracks(px, py, pz, tracksel, 0)
# tracks = alian.numpy_ptetaphi_to_tracks(px, py, pz, tracksel, 0)
# psj1 = alian.numpy_pxpypz_to_pseudojets(px, py, pz, 0.1359, 0)
# psj1 = alian.numpy_ptetaphi_to_pseudojets(pt, eta, phi, 0.1359, 0)
# psj1 = alian.numpy_ptetaphi_to_tracks(pt, eta, phi, tracksel, 0)
# psj1 = alian.numpy_pxpypz_to_tracks(px, py, pz, tracksel, 0)
# psj2 = alian.numpy_ptetaphi_to_pseudojets(pt, eta, phi, 0.1359, 0)
# psjlist = [psj for psj in psj1]
# psj1[0].set_user_info(alian.TrackInfo())
# print(type(psj1))
# print(type(psj1[0]))
# print(psj1[0].has_user_info())
# print(psj1[0].has_user_info[alian.TrackInfo]())
# # print(type(psj1[0].user_info_shared_ptr()))
# print(type(psj1[0].user_info[alian.TrackInfo]()))
# print(psj1[0].user_info[alian.TrackInfo]().q())
# print(psj1[0].user_info[alian.TrackInfo]().tracksel())
# print(psj1[1].user_info[alian.TrackInfo]().q())
# print(psj1[1].user_info[alian.TrackInfo]().tracksel())
# print(type(psjlist))
# print(type(psjlist[0]))

# jet_def = fj.JetDefinition(fj.antikt_algorithm, 0.4)
# cs = fj.ClusterSequence(psj1, jet_def)
# jets = fj.sorted_by_pt(cs.inclusive_jets())
# print(len(jets))
# print(type(jets))
# print(type(jets[0].constituents()[0]))
# print(jets[0].constituents()[0].user_info[alian.TrackInfo]().tracksel())
# print(jets[0].constituents()[0].user_info[alian.TrackInfo]().q())
# clusters = alian.numpy_energyetaphi_to_clusters(E, eta, phi, m02, m20, ncells, time, isexotic, distancebadchannel, nlm, clusterdef, matchedTrackIndex, 0)
# print(clusters)
# print(type(clusters))
# print(clusters[0].E())
# print(clusters[0].m())
# print(clusters[0].pt())
# print(clusters[0].eta())
# print(clusters[0].phi())
# print(clusters[0].delta_R(tracks[0]))
# print(TrackSel.globalTrack & tracks[0].tracksel())
# sys.exit(0)
# print(tracks[0].delta_eta(tracks2[0]))