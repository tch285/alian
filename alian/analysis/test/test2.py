#!/usr/bin/env python


# import heppyy
# fj = heppyy.load_cppyy('fastjet')

# import alian.analysis.base
# import heppyy.util.fastjet_cppyy # noqa: F401
# # import heppyy.util.pythia8_cppyy
# # import heppyy.util.heppyy_cppyy

# from cppyy.gbl import fastjet as fj

# print(fj.antikt_algorithm)


# import cppyy
# import sys
import numpy as np

# Load your C++ code
# cppyy.add_include_path(np.get_include())
# cppyy.load_library("libalian_track.so")
# cppyy.include("track/track.hh")
# cppyy.load_library("libalian_fjutil.so")
# cppyy.include("fjutil/fjutil.hh")
# from cppyy.gbl import alian
# import heppyy
# alian = heppyy.load_cppyy("alian")

from alian.analysis.base import JetFinder
# import heppyy.util.fastjet_cppyy  # noqa: F401
# from cppyy.gbl import fastjet as fj
# from cppyy.gbl import std
import heppyy
alian = heppyy.load_cppyy("alian")


# px = np.array([1.0])
# py = np.array([2.0])
# pz = np.array([3.0])
# E  = np.array([9.0])
# tracksel = np.array([15], dtype=np.uint8)
# arr1 = np.array([], dtype = np.float64)
# tracks = alian.numpy_pxpypz_to_tracks(px, py, pz, tracksel, 0.1359)
# tracks = alian.numpy_pxpypz_to_pseudojets(px, py, pz, 0.1359, 0)
pt = np.array([2.0], dtype = np.float32)
eta = np.array([0.3], dtype = np.float32)
phi = np.array([0.5], dtype = np.float32)
print(pt.dtype)
tracks = alian.numpy_ptetaphi_to_pseudojets(pt, eta, phi, 0.1359, 0)
print(tracks)
print(tracks[0].pt())
print(tracks[0].m())
print(tracks[0].E())

# sys.exit(0)


# Create a NumPy array
# arr1 = np.array([1.0, 2.0, 3.0, 4.0])
# print(arr1.shape)
# arr = np.array([65535, 65535, 65535, 65535, 65535], dtype=np.uint16)

# Pass to C++ using __array_interface__
# alian.process_array(arr1, *arr1.shape, arr, *arr.shape)
# print(arr1)  # [2., 4., 6., 8.]
# print(arr)  # [2., 4., 6., 8.]

# For 2D arrays
# arr_2d = np.array([[1.0, 2.0], [3.0, 4.0, 5.0]])
# alian.process_2d_array(
#     arr_2d, *arr_2d.shape
# )

# print(arr_2d)


# import cppyy
# import numpy as np

# cppyy.cppdef("""
# void process_array(double* data, int size) {
#     for (int i = 0; i < size; ++i) {
#         data[i] *= 2.0;
#     }
# }
# """)

# arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)

# # Direct call - cppyy should convert automatically
# cppyy.gbl.process_array(arr, len(arr))
# print(arr)