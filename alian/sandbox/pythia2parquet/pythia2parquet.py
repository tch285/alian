#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import json
import yasp
import cppyy


import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

# from cppyy.gbl import pythiaext

from heppyy.pythia_util import configuration as pyconf
import logging
from heppyy.util.logger import Logger
log = Logger()

import ROOT
import math
import array
import itertools


sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from alian.io.root_io import SingleRootFile


def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr


def idx_match_to_out_parton(jet, out_parton1, out_parton2, jet_R0=0.4):
    # print(jet, out_parton1, out_parton2)
    # print(jet.delta_R(out_parton1), jet.delta_R(out_parton2))
    # if jet.delta_R(out_parton1) < jet_R0:
    # 		return out_parton1.user_index()
    # if jet.delta_R(out_parton2) < jet_R0:
    # 		return out_parton2.user_index()
    if jet.delta_R(out_parton1) < jet.delta_R(out_parton2):
        return out_parton1.user_index()
    else:
        return out_parton2.user_index()
    return None

# we measured for k=1 and a>0
def angularity(jet, a, k, jetR):
	ang = 0.0
	for p in jet.constituents():
		dr = jet.delta_R(p) / jetR
		pt = p.perp() / jet.perp()
		ang += ((dr)**a) * ((pt)**k)
	return ang


def main():
  parser = argparse.ArgumentParser(description='pythia 2 parquet', prog=os.path.basename(__file__))
  pyconf.add_standard_pythia_args(parser)
  parser.add_argument('-o','--output', help='root output filename', default='pythia8_simple_eec_output.root', type=str)
  parser.add_argument('--charged', help='charged particles only', action='store_true')
  parser.add_argument('--jet-pt-min', help='jet pT min', type=float, default=20.0)
  parser.add_argument('--jet-pt-max', help='jet pT min', type=float, default=1000.0)
  parser.add_argument('--part-eta-max', help='part eta max', type=float, default=1.0)
  parser.add_argument('--debug', help='debug', action='store_true')
  parser.add_argument('--write-config', help='write config to yaml and quit', type=str, default='')
  parser.add_argument('--jet-R', help='jet R', type=float, default=0.4)
  args = parser.parse_args()	
 
  if args.debug:
    log.set_level(logging.DEBUG)

  if args.output == 'pythia8_simple_eec_output.root':
    if args.py_vincia:
      args.output = args.output.replace('.root', '_vincia.root')
    if args.py_dire:
      args.output = args.output.replace('.root', '_dire.root')
    print("[w] using [modified] default output file:", args.output)
  else:
    print("[w] using specified output file:", args.output)

  rf = SingleRootFile(fname=args.output)
  rf.root_file.cd()
  h_jet_pt  = ROOT.TH1F('h_jet_pt', 'h_jet_pt', 100, 0, 100)
  tn_norm   = ROOT.TNtuple('tn_norm', 'tn_norm', 'nev:xsec:xsec_err:sum_of_weights')
  tn_hard 	= ROOT.TNtuple('tn_hard', 'tn_hard', 'nev:xsec:ev_weight:x1:x2:QFac:id1:id2:id3:pt3:eta3:id4:pt4:eta4')
  tn_events = ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts:x1:x2:QFac:ncoll')
  tn_jet 		= ROOT.TNtuple(f'tn_jet', 'tn_jet', 'nev:xsec:ev_weight:nj:ij:pt:eta:phi:m:ptlead:pid')

  # jet finder
  # print the banner first
  fj.ClusterSequence.print_banner()
  print()
  jet_R0 = args.jet_R
  jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
  jet_selector = fj.SelectorPtMin(args.jet_pt_min)
  jet_selector = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_max) * fj.SelectorAbsEtaMax(args.part_eta_max - jet_R0 * 1.05)
 
  mycfg = []
  pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
  if not pythia:
    print("[e] pythia initialization failed.")
    return
  
  jets_out = []

  sum_weights = 0.0
  _stop = False 
  pbar = tqdm.tqdm(range(args.nev))
  njets = 0
  iev = 0
  while not _stop:
    if not pythia.next():
      continue
    # _info = Pythia8.pythia.info - defunct
    _info = Pythia8.getInfo(pythia)
    sigmaGen = _info.sigmaGen()
    ev_weight = _info.weight()
    sum_weights += ev_weight
    iev += 1

    if args.charged:
      fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
    else:
      fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])

    jets = fj.sorted_by_pt(jet_selector(jet_def(fjparts)))
    njets += len(jets)
    # see https://pythia.org/latest-manual/HeavyIons.html 
    ncoll = 1
    hiinfo = Pythia8.getHIInfo(pythia)
    if hiinfo:
      ncoll = hiinfo.nCollTot()   
    if jets.size() > 0:
      # tn_events 	= ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts')
      tn_events.Fill(iev, sigmaGen, ev_weight, pythia.event.size(), _info.x1(), _info.x2(), _info.QFac(), ncoll)
      in_parton1 = pythia.event[3]
      in_parton2 = pythia.event[4]
      out_parton1 = pythia.event[5]
      out_parton2 = pythia.event[6]
      fj_out_partons = vector[fj.PseudoJet]()
      for p in [out_parton1, out_parton2]:
        fjp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
        fjp.set_user_index(p.index())
        fj_out_partons.push_back(fjp)
      # tn_hard = ROOT.TNtuple('tn_hard', 'tn_hard', 'nev:xsec:ev_weight:x1:x2:QFac:id1:id2:id3:id4:pt3:eta3:pt4:eta4')
      tn_hard.Fill(iev, sigmaGen, ev_weight, _info.x1(), _info.x2(), _info.QFac(), in_parton1.id(), in_parton2.id(), out_parton1.id(), out_parton1.pT(), out_parton1.eta(), out_parton2.id(), out_parton2.pT(), out_parton2.eta())
      for ij, j in enumerate(jets):
        h_jet_pt.Fill(j.perp(), ev_weight)
        pid_idx = idx_match_to_out_parton(j, fj_out_partons[0], fj_out_partons[1], jet_R0)
        pid = pythia.event[pid_idx].id() if pid_idx is not None else 0
        constituents = fj.sorted_by_pt(j.constituents())
        ptlead = constituents[0].perp()
        # print(j.perp())
        # tn_jet = ROOT.TNtuple(f'tn_jet', 'tn_jet', 'nev:xsec:ev_weight:nj:ij:pt:eta:phi:m:ptlead')				
        tn_jet.Fill(iev, sigmaGen, ev_weight, njets, ij, j.perp(), j.eta(), j.phi(), j.m(), ptlead, pid)
        part_dicts = []
        for c in constituents:
            part_dicts.append({'pt': float(np.float32(c.perp())), 'eta': float(np.float32(c.eta())), 'phi': float(np.float32(c.phi())), 'm': float(np.float32(c.m()))})
        j_dict = {'nev': iev, 'xsec': float(np.float32(sigmaGen)), 'ev_weight': float(np.float32(ev_weight)), 'nj': njets, 'ij': ij, 'pt': float(np.float32(j.perp())), 'eta': float(np.float32(j.eta())), 'phi': float(np.float32(j.phi())), 'm': float(np.float32(j.m())), 'ptlead': float(np.float32(ptlead)), 'pid': pid, 'constituents': part_dicts,
                  'angk1a1': float(np.float32(angularity(j, 1.0, 1.0, jet_R0))),
									'angk1a2': float(np.float32(angularity(j, 2.0, 1.0, jet_R0))),
									'angk1a3': float(np.float32(angularity(j, 3.0, 1.0, jet_R0)))
									}
        jets_out.append(j_dict)
    else:
      continue
    pbar.update(jets.size())
    if pbar.n >= args.nev:
      _stop = True

  pythia.stat()

  # Get the total cross section and weight sum
  sigma_gen = _info.sigmaGen()
  sigma_gen_err = _info.sigmaErr()
  weight_sum = _info.weightSum()  # Same as sum_weights
  nAccepted = _info.nAccepted()

  # tn_norm =   ROOT.TNtuple('tn_norm', 'tn_norm', 'nev:xsec:xsec_err:sum_of_weights')
  tn_norm.Fill(nAccepted, sigma_gen, sigma_gen_err, weight_sum)

  # Save to JSON summary
  json_file = args.output.replace('.root', '.json')
  with open(json_file, "w") as f:
      json.dump({
          "n_accepted": nAccepted,
          "sigma_gen": sigma_gen,
          "sum_weights": weight_sum
      }, f, indent=2)

  # Output for verification
  print(f"sigma_gen = {sigma_gen:.3f} +- {sigma_gen_err:.3f} mb")
  print(f"Sum of weights = {weight_sum:.3f} [check: {sum_weights:.3f}]")
 
  print('[i] number of jets:', njets)
  
  h_jet_pt.Scale(sigma_gen / sum_weights)
  rf.close()
  
  #now save also to parquet
  import pandas as pd
  df = pd.DataFrame(jets_out)
  
  # convert to more efficient types
  # # Option 1: Convert float64 columns to float32 (half precision)
  # float_columns = ['xsec', 'ev_weight', 'pt', 'eta', 'phi', 'm', 'ptlead', 'angk1a1', 'angk1a2', 'angk1a3']
  # for col in float_columns:
  #   if col in df.columns:
  #     df[col] = df[col].astype('float32')
  # 
  # # Also convert constituent data to float32
  # if 'constituents' in df.columns:
  #   for idx, constituents in enumerate(df['constituents']):
  #     if constituents:
  #       for part in constituents:
  #         for key in ['pt', 'eta', 'phi', 'm']:
  #           if key in part:
  #             part[key] = float(np.float32(part[key]))
  
  pq_file = args.output.replace('.root', '.parquet')
  # Save with compression and additional options for better storage efficiency
  df.to_parquet(pq_file, compression='snappy', index=False)

if __name__ == '__main__':
  main()
