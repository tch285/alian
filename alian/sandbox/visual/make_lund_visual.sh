#!/bin/bash

./pythia_jets_lund2json.py 
./plot_jet_json.py

./pythia_jets_lund2json.py --ptMin 1000
./plot_jet_json.py
