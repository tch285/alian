#!/usr/bin/env python

import yasp
yasp.module_load_cppyy('bundle/hepbase')
yasp.module_load_cppyy('heppyy/current')
yasp.module_load_cppyy('alian/current')

import heppyy
fj = heppyy.load_cppyy('fastjet')
std = heppyy.load_cppyy('std')
Pythia8 = heppyy.load_cppyy('pythia8.Pythia8')
alian_cpp = heppyy.load_cppyy('alian')

import alian

from heppyy.pythia_util import configuration as pyconf

import math
import ROOT
import argparse
import tqdm

class JetAnalyzer:
    def __init__(self, histFile="jetAnalysis.root", ptMin=20.0, jetR=0.4):
        # Book ROOT histograms.
        self.histFile = histFile
        self.ptMin = ptMin
        self.jetR = jetR
        self.hFrag = ROOT.TH1F("hFrag", "Fragmentation Function;z;dP/dx", 100, 0.0, 1.0)
        self.hAng = ROOT.TH1F("hAng", "Jet Angularity;Angularity;dP/dx", 100, 0.0, 1.0)
        self.hDr = ROOT.TH1F("hDr", "Constituent angles;#Delta R / R_{jet};dP/dx", 100, 0.0, 1.0)

    def analyze_event(self, pythia):
        """
        Convert final state Pythia particles to FastJet PseudoJets,
        cluster jets using anti-kt algorithm, and fill histograms.
        """
        # Collect final state particles.
        fjInputs = std.vector(fj.PseudoJet)([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])
        # Define jet clustering: anti-kt algorithm with radius self.jetR.
        event = pythia.event
        
        # Define jet clustering: anti-kt algorithm with radius self.jetR.
        jetDef = fj.JetDefinition(fj.antikt_algorithm, self.jetR)
        cs = fj.ClusterSequence(fjInputs, jetDef)
        # Get jets above the pt threshold.
        jets = cs.inclusive_jets(self.ptMin)
        # Sort jets in descending order of pT.
        jets = sorted(jets, key=lambda jet: jet.pt(), reverse=True)

        # Loop through jets and fill histograms.
        for jet in jets:
            jetPt = jet.pt()
            constituents = jet.constituents()

            # Fill fragmentation histogram: z = constituent pT / jet pT.
            for p in constituents:
                z = p.pt() / jetPt
                self.hFrag.Fill(z)

            # Compute a simple jet angularity:
            # Sum_i (pT_i * DeltaR(i, jet)) divided by jet pT.
            ang = 0.0
            for p in constituents:
                dr = jet.delta_R(p)
                self.hDr.Fill(dr/self.jetR)
                ang += p.pt() * dr
            if jetPt > 0:
                ang /= jetPt
            self.hAng.Fill(ang)

    def save_histograms(self):
        """
        Write histograms to a ROOT file.
        """
        print('Saving histograms to', self.histFile)
        outFile = ROOT.TFile(self.histFile, "RECREATE")
        norm = self.hAng.GetEntries()
        if norm > 0:
            self.hAng.Scale(1.0 / norm / self.hAng.GetBinWidth(1))
            self.hDr.Scale(1.0 / norm / self.hDr.GetBinWidth(1))
            self.hFrag.Scale(1.0 / norm / self.hFrag.GetBinWidth(1))
        self.hFrag.Write()
        self.hAng.Write()
        self.hDr.Write()
        outFile.Close()
        print("Histograms saved to", self.histFile)
        
    def __del__(self):
        """
        Destructor: save histograms before deleting the object.
        """
        self.save_histograms()
        

def main():
    argparser = argparse.ArgumentParser(description="Pythia8 jet analysis")
    argparser.add_argument("-o", "--output", default="h_test_custom_shower.root", help="Output ROOT file")
    argparser.add_argument("--ptMin", type=float, default=100.0, help="Minimum jet pT")
    argparser.add_argument("--jetR", type=float, default=0.4, help="Jet radius")
    argparser.add_argument("--nev", type=int, default=1000, help="Number of events")
    argparser.add_argument("--zBias", type=float, default=0.0, help="zBias")
    argparser.add_argument("--coherenceFactor", type=float, default=0.0, help="coherenceFactor")
    argparser.add_argument("--azimuthalBias", type=float, default=0.0, help="azimuthalBias")
    argparser.add_argument("--recoilWeight", type=float, default=0.0, help="recoilWeight")
    argparser.add_argument("--angleBias", type=float, default=0.0, help="angleBias")
    argparser.add_argument("--alphaModifier", type=float, default=0.0, help="alphaModifier")
    
                          
    
    args = argparser.parse_args()
    # only in version 8.313
    # current_shower = Pythia8.getShowerModelPtr()
    # print(current_shower)

    myShower = Pythia8.SimpleShowerModelTCustomized()
    myShowerPtr = Pythia8.ShowerModelPtr(myShower)
    #myShower = Pythia8.SimpleTimeShowerCustomized()
    #myShower.getTimeShowerPtr()

    #import pySimpleTimeShower
    #from pySimpleTimeShower import Pythia, SimpleTimeShowerCustomized

    # Create Pythia instance
    pythia = Pythia8.Pythia()
    myShower.addParameters(pythia.settings)
    # Register custom time shower
    # shower = SimpleTimeShowerCustomized()
    pythia.setShowerModelPtr(myShowerPtr)
    extras = ["Next:numberCount = 0", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0", "Stat:showProcessLevel = on"]
    mycfg = ['PhaseSpace:pThatMin = 100', 'HardQCD:all = on', "Beams:eCM = 13000."]
    mycfg.extend(extras)
    # Initialize with settings
    for c in mycfg:
        pythia.readString(c)
        
    # filepath: /Users/ploskon/devel/alian/alian/sandbox/test_custom_shower.py
    pythia.readString("SimpleTimeShowerCustomized:zBias = 0.0")
    pythia.readString("SimpleTimeShowerCustomized:coherenceFactor = 0.0")
    pythia.readString("SimpleTimeShowerCustomized:azimuthalBias = 0.0")
    pythia.readString("SimpleTimeShowerCustomized:recoilWeight = 0.0")
    pythia.readString("SimpleTimeShowerCustomized:angleBias = 0.0")
    pythia.readString("SimpleTimeShowerCustomized:alphaModifier = 0.0")

    pythia.init()
    if not pythia:
        print("[e] pythia initialization failed.")
        raise RuntimeError("Pythia initialization failed.")

    jana = JetAnalyzer(histFile=args.output, ptMin=args.ptMin, jetR=args.jetR)
    # Event loop
    # for i in range(args.nev):
    for i in tqdm.tqdm(range(args.nev)):
        if not pythia.next():
            continue
        jana.analyze_event(pythia)

    # Show stats
    pythia.stat()

    pythia.settings.writeFile('test_custom_shower_all.cmnd', True)
    pythia.settings.writeFile('test_custom_shower_changed.cmnd', False)
    
if __name__ == "__main__":
    main()