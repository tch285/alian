#!/bin/bash

nev=100000
jet_pt_min=20

# ./pythia8_simple_eec.py --jet-R 0.2 --output pythia8_simple_eec_output_nobias.root 	--py-cmnd pythia8_pp_nobias.cmnd	--nev ${nev} --jet-pt-min ${jet_pt_min} --py-pthatmin ${jet_pt_min}

./pythia8_simple_eec.py --jet-R 0.2 --output pythia8_simple_eec_output.root --py-cmnd pythia8_pp.cmnd	--nev ${nev} --jet-pt-min ${jet_pt_min} --py-pthatmin ${jet_pt_min}

#  --jet-pt-min ${jet_pt_min} --py-biasref 20
