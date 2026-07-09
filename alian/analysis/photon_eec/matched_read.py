# This script is used on the ROOT file (provided as an argument) in
# the BerkeleyTree format (YAML file with structure must be provided
# with -t) to perform the analysis of the events containing isolated
# photons and jets.
# The script performs necessary cuts specified in YAML file (enabled
# with -c), creates QA histograms for the cleaned data, and stores them
# in the output ROOT file (specified with -o)

#!/usr/bin/env python

# import argparse
# import itertools
# import heppyy
# fj = heppyy.load_cppyy('fastjet')
from alian.analysis.base.selection import TrackSel
from alian.analysis.base import Event
# import numpy as np
# from alian.analysis.base.utils import delta_phi_limit, delta_R
from alian.io import data_io
# import numpy as np
import heppyy
# import heppyy.util.fastjet_cppyy  # noqa: F401
alian = heppyy.load_cppyy("alian")
# from cppyy.gbl import fastjet as fj

# input_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/run3/data/staging/mhwang/convert_test/BerkeleyTree.root"
input_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/run3/data/LHC24ao/BerkeleyTrees/1/BerkeleyTree.root"
data_source = data_io.DataInput(input_file,
                lhc_run = 3,
                yaml_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/run3/data/LHC24ao/BerkeleyTrees/tstruct.yaml",
                n_events = 1000)


for e in data_source.next_event(disable_bar = True):
    # build the event only
    event = Event(e)
    if len(e.data["cluster_time"]) > 0 and len(e.data["cluster_matched_track_delta_eta"]) > 1:
        print("event:")
        print(e.data["cluster_time"])
        # print(e.data["cluster_ncells"])
        print(e.data["cluster_matched_track_n"])
        print(e.data["cluster_matched_track_delta_eta"])
        # print(e.data["cluster_matchedTrackPhi"])
        # print(e.data["cluster_matchedTrackP"])
        # print(type(e.data["cluster_matchedTrackDeltaEta"]))
        # print(e)
        clus = alian.numpy_energyetaphi_to_clusters(
            e.data["cluster_energy"],
            e.data["cluster_eta"],
            e.data["cluster_phi"],
            e.data["cluster_m02"],
            e.data["cluster_m20"],
            e.data["cluster_ncells"],
            e.data["cluster_time"],
            e.data["cluster_exoticity"],
            e.data["cluster_dbc"],
            e.data["cluster_nlm"],
            e.data["cluster_defn"],
            e.data["cluster_matched_track_n"],
            e.data["cluster_matched_track_delta_eta"],
            e.data["cluster_matched_track_delta_phi"],
            e.data["cluster_matched_track_p"],
            e.data["cluster_matched_track_pt"],
            e.data["cluster_matched_track_sel"],
            0
        )
        # print(type(clus))
        print("nclus:", clus.size())
        for cluster in clus:
            N = cluster.matchedTrackN()
            print("matched eta clusterN:", N)
            print(cluster.matchedTrackDeltaEta())
            # for i in range(N):
            #     print(cluster.matchedTrackDeltaEta()[i])
            #     print(cluster.matchedTrackDeltaPhi()[i])
            #     print(cluster.matchedTrackSel()[i])
            #     print(TrackSel(cluster.matchedTrackSel()[i]))

        # break
        # print('enxt')
        # print(type(e.data["cluster_matchedTrackEta"].__array__()))
        # print(type(e.data["cluster_matchedTrackEta"][0]))
        # print(type(e.data["cluster_matchedTrackEta"][0].__array__()))
        # print(type(e.data["cluster_matchedTrackEta"][0][0].__array__()))
        # print(type(e.data["cluster_matchedTrackEta"]._values))
        # print(type(e.data["cluster_matchedTrackEta"].tolist()))
        # print(type(e.data["cluster_matchedTrackEta"][0]))
        # print(dir(e.data["cluster_matchedTrackEta"]))
        # newarr = np.asarray(e.data["cluster_matchedTrackEta"].__array__())
        # newarr = np.array([np.asarray(row, dtype=np.float64) for row in e.data["cluster_matchedTrackEta"].__array__()])
        # print(newarr)
        # print(newarr.dtype)
        # print(newarr[0])
        # print(newarr[0].dtype)
        # alian.test(newarr)
        # newarr = e.data["cluster_matchedTrackEta"].tolist()
        # print(type(newarr))
        # print(type(newarr[0]))
        # print(type(newarr[0][0]))
        # alian.testvec(e.data["cluster_matchedTrackEta"].tolist())
    # print(event)
    # skip if the event doesn't pass the selection
#     if not self.selector.event.selects(self.event):
#         continue
#     # build the rest of the event and analyze
#     self.build_event_objs(e)
#     self.analyze_event()
# self.note_time("Analyzed all events")
# import cppyy
# import cppyy.gbl as gbl

# cppyy.cppdef("""
# #include <iostream>
# #include <vector>
# #include <string>

# // Function that accepts a list of lists (represented as a vector of vectors in C++)
# void process_nested_list(const std::vector<std::vector<int>>& data) {
#     std::cout << "Received nested list with " << data.size() << " inner lists." << std::endl;
#     for (size_t i = 0; i < data.size(); ++i) {
#         std::cout << "  Inner list " << i << " has " << data[i].size() << " elements: ";
#         for (int value : data[i]) {
#             std::cout << value << " ";
#         }
#         std::cout << std::endl;
#     }
# }
# """)

# Access the C++ function from the global namespace
# process_nested_list = gbl.process_nested_list

# # A Python list
# python_nested_list = [
#     [1],
#     [4, 5],
#     [6, 7, 8, 9]
# ]

# # Pass the Python list directly to the C++ function
# process_nested_list(python_nested_list)