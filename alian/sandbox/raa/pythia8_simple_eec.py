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


def perpendicular_cone_particles(jet, particles, R0=0.4):
  # get phi of a rotated jet in phi by pi/2
  jet_phi = jet.phi()
  jet_phi_perp_plus_min = jet_phi + math.pi/2.
  jet_phi_perp_plus_max = jet_phi + math.pi/2. + R0
  jet_phi_perp_minus_min = jet_phi - math.pi/2. - R0
  jet_phi_perp_minus_max = jet_phi - math.pi/2.
  for p in particles:
    #delta_phi returns -pi to pi
    if p.delta_phi_to(jet) > jet_phi_perp_minus_min and p.delta_phi_to(jet) < jet_phi_perp_minus_max:
      p.set_user_index(-1)
      yield p
    if p.delta_phi_to(jet) > jet_phi_perp_plus_min and p.delta_phi_to(jet) < jet_phi_perp_plus_max:
      p.set_user_index(-2)
      yield p

def perpendicular_part_sets(jet, particles, R0, ptcut):
  # jet_phi = jet.phi()
  jet_dphi_perp_plus_min   = 0 + math.pi/2. - R0/2.
  jet_dphi_perp_plus_max   = 0 + math.pi/2. + R0/2.
  jet_dphi_perp_minus_min  = 0 - math.pi/2. - R0/2.
  jet_dphi_perp_minus_max  = 0 - math.pi/2. + R0/2.
  jconsts = [p for p in jet.constituents() if p.perp() >= ptcut]
  pcone_plus = vector[fj.PseudoJet]()
  pcone_minus = vector[fj.PseudoJet]()
  for p in jconsts:
    p.set_user_index(0)
    pcone_minus.push_back(p)
    pcone_plus.push_back(p)
  for p in particles:
    if p.perp() >= ptcut:
      ptmp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
      if p.delta_phi_to(jet) > jet_dphi_perp_plus_min and p.delta_phi_to(jet) < jet_dphi_perp_plus_max:
        # reset_PtYPhiM (double pt, double y, double phi, double m=0.0)
        ptmp.reset_PtYPhiM(p.perp(), p.rap(), p.phi() - math.pi/2., p.m())
        ptmp.set_user_index(-1)
        pcone_plus.push_back(ptmp)
      if p.delta_phi_to(jet) > jet_dphi_perp_minus_min and p.delta_phi_to(jet) < jet_dphi_perp_minus_max:
        ptmp.reset_PtYPhiM(p.perp(), p.rap(), p.phi() + math.pi/2., p.m())
        ptmp.set_user_index(-2)
        pcone_minus.push_back(ptmp)
  # print(' ->', ptcut, len(jconsts), len(pcone_plus), len(pcone_minus))
  return pcone_plus, pcone_minus


def main():
  parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
  pyconf.add_standard_pythia_args(parser)
  parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
  parser.add_argument('-o','--output', help='root output filename', default='pythia8_simple_eec_output.root', type=str)
  parser.add_argument('--charged', help='charged particles only', action='store_true')
  parser.add_argument('--jet-pt-min', help='jet pT min', type=float, default=20.0)
  parser.add_argument('--jet-pt-max', help='jet pT min', type=float, default=1000.0)
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
  # tn_pi0 = ROOT.TNtuple(f'tn_pi0', 'tn_pi0', 'nev:xsec:ev_weight:pt:eta:phi:m')
  # tn_pich = ROOT.TNtuple(f'tn_pich', 'tn_pich', 'nev:xsec:ev_weight:pt:eta:phi:m')
  tn_pich_jev = ROOT.TNtuple(f'tn_pich_jev', 'tn_pich_jev', 'nev:xsec:ev_weight:pt:eta:phi:m:ptjet')
  # pt_cuts = [0.15, 1.0]
  pt_cuts = [1.0]
  tn_eec = {}
  tn_eec_format_str = 'nev:xsec:ev_weight:ij:dr:pt1:pt2:eec:ptjet:ptcut:iidx:jidx:ptlead:pid:dp'
  for ptcut in pt_cuts:
    tn_eec[ptcut] = ROOT.TNtuple(f'tn_eec_ptcut{ptcut}', 'tn_eec_ptcut{ptcut}', tn_eec_format_str)
  tn_eec_pcone = {}
  for ptcut in pt_cuts:
    tn_eec_pcone[ptcut] = ROOT.TNtuple(f'tn_eec_pcone_ptcut{ptcut}', 'tn_eec_pcone_ptcut{ptcut}', tn_eec_format_str)

  pythia = Pythia8.Pythia()

  # jet finder
  # print the banner first
  fj.ClusterSequence.print_banner()
  print()
  jet_R0 = args.jet_R
  jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
  jet_selector = fj.SelectorPtMin(args.jet_pt_min)
  jet_selector = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_max) * fj.SelectorAbsEtaMax(1 - jet_R0 * 1.05)
 
  mycfg = []
  pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
  if not pythia:
    print("[e] pythia initialization failed.")
    return

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
    # tn_pi0 = ROOT.TNtuple(f'tn_pi0', 'tn_pi0', 'nev:xsec:ev_weight:pt:eta:phi:m')
    # _ = [tn_pi0.Fill(iev, sigmaGen, ev_weight, p.pT(), p.eta(), p.phi(), p.m()) for p in pythia.event if p.id() == 111]
    # _ = [tn_pich.Fill(iev, sigmaGen, ev_weight, p.pT(), p.eta(), p.phi(), p.m()) for p in pythia.event if abs(p.id()) == 211 and p.isFinal()]
    if jets.size() > 0:
      # tn_events 	= ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts')
      tn_events.Fill(iev, sigmaGen, ev_weight, pythia.event.size(), _info.x1(), _info.x2(), _info.QFac(), ncoll)
      # tn_pich_jev = ROOT.TNtuple(f'tn_pich_jev', 'tn_pich_jev', 'nev:xsec:ev_weight:pt:eta:phi:m:ptjet')
      _ = [tn_pich_jev.Fill(iev, sigmaGen, ev_weight, p.pT(), p.eta(), p.phi(), p.m(), jets[0].perp()) for p in pythia.event if abs(p.id()) == 211 and p.isFinal()]
      # _ = [print(p.index(), p.id(), p.status(), p.pT(), p.eta(), p.phi(), p.m()) for p in pythia.event]
      in_parton1 = pythia.event[3]
      in_parton2 = pythia.event[4]
      out_parton1 = pythia.event[5]
      out_parton2 = pythia.event[6]
      fj_out_partons = vector[fj.PseudoJet]()
      for p in [out_parton1, out_parton2]:
        fjp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
        fjp.set_user_index(p.index())
        fj_out_partons.push_back(fjp)
      for ij, j in enumerate(jets):
        h_jet_pt.Fill(j.perp(), ev_weight)
        pid_idx = idx_match_to_out_parton(j, fj_out_partons[0], fj_out_partons[1], jet_R0)
        pid = pythia.event[pid_idx].id() if pid_idx is not None else 0
        ptlead = fj.sorted_by_pt(j.constituents())[0].perp()
        # print(j.perp())
        # tn_jet = ROOT.TNtuple(f'tn_jet', 'tn_jet', 'nev:xsec:ev_weight:nj:ij:pt:eta:phi:m:ptlead')				
        tn_jet.Fill(iev, sigmaGen, ev_weight, njets, ij, j.perp(), j.eta(), j.phi(), j.m(), fj.sorted_by_pt(j.constituents())[0].perp(), pid)
        # tn_hard = ROOT.TNtuple('tn_hard', 'tn_hard', 'nev:xsec:ev_weight:x1:x2:QFac:id1:id2:id3:id4:pt3:eta3:pt4:eta4')
        tn_hard.Fill(iev, sigmaGen, ev_weight, _info.x1(), _info.x2(), _info.QFac(), in_parton1.id(), in_parton2.id(), out_parton1.id(), out_parton1.pT(), out_parton1.eta(), out_parton2.id(), out_parton2.pT(), out_parton2.eta())
        for ptcut in pt_cuts:
            _parts_cut = [p for p in j.constituents() if p.perp() >= ptcut]
            _pairs = list(itertools.product(_parts_cut, repeat=2))
            if len(_pairs) < 1:
                continue
            for first, second in _pairs:
                dr = first.delta_R(second)
                eec = first.perp() * second.perp() / pow(j.perp(), 2.)
                dp = math.sqrt((first - second).modp2())
                # tn_eec[ptcut] = ROOT.TNtuple(f'tn_eec_ptcut{ptcut}', 'tn_eec_ptcut{ptcut}', 'nev:xsec:ev_weight:ij:dr:pt1:pt2:eec:ptjet:ptcut:iidx:jidx:ptlead')
                tn_eec[ptcut].Fill(iev, sigmaGen, ev_weight, ij, dr, first.perp(), second.perp(), eec, j.perp(), ptcut, first.user_index(), second.user_index(), ptlead, pid, dp)

            # perpendicular cones
            _parts_cut_pcone_plus, _parts_cut_pcone_minus = perpendicular_part_sets(j, fjparts, jet_R0, ptcut)
            for pcone in [_parts_cut_pcone_plus, _parts_cut_pcone_minus]:
              _pairs_pcone = list(itertools.product(pcone, repeat=2))
              if len(_pairs_pcone) < 1:
                  continue
              for first, second in _pairs_pcone:
                  if first.user_index() >= 0 and second.user_index() >= 0:
                      continue
                  # print('accepting pair:', first.user_index(), second.user_index())
                  dr = first.delta_R(second)
                  eec = first.perp() * second.perp() / pow(j.perp(), 2.)
                  dp = math.sqrt((first - second).modp2())
                  # tn_eec[ptcut] = ROOT.TNtuple(f'tn_eec_ptcut{ptcut}', 'tn_eec_ptcut{ptcut}', 'nev:xsec:ev_weight:ij:dr:pt1:pt2:eec:ptjet:ptcut:iidx:jidx:ptlead')
                  tn_eec_pcone[ptcut].Fill(iev, sigmaGen, ev_weight, ij, dr, first.perp(), second.perp(), eec, j.perp(), ptcut, first.user_index(), second.user_index(), ptlead, pid, dp)
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

if __name__ == '__main__':
  main()
