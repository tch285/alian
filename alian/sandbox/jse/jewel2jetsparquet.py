#!/usr/bin/env python

import uproot
import pandas as pd

from yasp import GenericObject
import heppyy.util.fastjet_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl.std import vector

from heppyy.util.mputils import logbins
from lund_jet import LundJet

# now we want to use the reader class but we want to use the function apply to build events on the fly without loading the entire file into memory
class JewelReader:
		def __init__(self, file_path):
				self.file_path = file_path

		def read_file(self):
				with uproot.open(self.file_path) as file:
						# Read tracks tree into a DataFrame
						tracks_df = file["tracks"].arrays(library="pd")
						# Read event_info tree into a DataFrame
						event_info_df = file["event_info"].arrays(library="pd")
				
				# Merge tracks and event_info on eventID
				events_df = pd.merge(tracks_df, event_info_df, on="eventID", how="inner")
				return events_df
			
		def build_events(self):
				events_df = self.read_file()
				for n, (event_id, event) in enumerate(events_df.groupby("eventID")):
						yield {
								"event_id": event_id,
								"num_tracks": len(event),
								"weight": event["weight"].iloc[0],
								"xsec": event["xsec"].iloc[0],
								"tracks": event[["label", "px", "py", "pz", "energy"]],
								# vector of pseudo-jets for fastjet
								"pseudo_jets": [fj.PseudoJet(px, py, pz, energy) for px, py, pz, energy in zip(event["px"], event["py"], event["pz"], event["energy"])],
								"vector_pseudo_jets": vector[fj.PseudoJet]([fj.PseudoJet(px, py, pz, energy) for px, py, pz, energy in zip(event["px"], event["py"], event["pz"], event["energy"])])
						}

def jewel_to_parquet(file_path, output=None, jetR=0.4, label=None, max_events=None):
		"""
		Read a JEWEL file, cluster jets, build LundJet dicts and write to parquet.

		Parameters
		----------
		file_path  : str   input JEWEL .root file
		output     : str   output parquet path (default: input basename + .parquet)
		jetR       : float jet radius (default 0.4)
		label      : str   jet label stored in LundJet
		max_events : int   stop after this many events (None = all)
		"""
		if output is None:
				import os
				output = os.path.splitext(os.path.basename(file_path))[0] + '.parquet'
				
		jet_selector = fj.SelectorPtMin(20.0) * fj.SelectorAbsRapMax(2.5)  # example: select jets with pt > 20 GeV and |rapidity| < 2.5
		jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
		jets_dicts = []

		reader = JewelReader(file_path)
		i_current = 0
		for event in reader.build_events():
			if max_events is not None and i_current >= max_events:
				break
			cseq = fj.ClusterSequence(event['vector_pseudo_jets'], jet_def)
			jets = fj.sorted_by_pt(jet_selector(cseq.inclusive_jets()))
			for j in jets:
				lj = LundJet(jet=j, jetR=jetR, label=label)
				jets_dicts.append(lj.to_basic_type_dict())
			i_current += 1
			if i_current % 1000 == 0:
				print(f"Processed {i_current} events...")
		df = pd.DataFrame(jets_dicts)
		print(f'number of jets: {len(df)}  ->  {output}')
		df.to_parquet(output, engine='pyarrow')
		return df

# example
# jewel_to_parquet(pp_100_GeV_file_path, output='pp_lund_jets.parquet', label='pp', max_events=100)

if __name__ == "__main__":
		import argparse
		parser = argparse.ArgumentParser(description='Convert JEWEL .root file to parquet with LundJet info')
		parser.add_argument('input', help='Input JEWEL .root file path')
		parser.add_argument('--output', help='Output parquet file path (default: input basename + .parquet)')
		parser.add_argument('--jetR', type=float, default=0.4, help='Jet radius for clustering (default: 0.4)')
		parser.add_argument('--label', default="", type=str, help='Label to store in LundJet (default: None)')
		parser.add_argument('--max-events', type=int, default=None, help='Maximum number of events to process (default: all)')
		args = parser.parse_args()

		jewel_to_parquet(args.input, output=args.output, jetR=args.jetR, label=args.label, max_events=args.max_events)