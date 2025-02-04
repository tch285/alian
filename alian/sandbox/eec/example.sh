#!/bin/bash

./pythia8_simple_eec.py --nev 1000 --py-cmnd ./pythia8_pp.cmnd --jet-pt-min 100 --jet-pt-max 120 --py-pthatmin 100
../draw_from_yaml.py -c tdraw_eec_dp.yaml -d input=pythia8_simple_eec_output.root output=hist.root jetptmin=100 jetptmax=120
