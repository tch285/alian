#!/bin/bash

inputfile=$1
jetptmin=$2
jetptmax=$3
outputfile=$4

if [ -z "$inputfile" ]; then
	inputfile=pythia8_simple_eec_output.root
fi

if [ -z "$jetptmin" ]; then
	jetptmin=20
fi

if [ -z "$jetptmax" ]; then
	jetptmax=150
fi

if [ -z "$outputfile" ]; then
	outputfile="${inputfile%.root}_${jetptmin}_${jetptmax}_h.root"
fi

draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.0 fquenchconst=0.0
outputfile_q="${outputfile%.root}_q0p075.root"
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile_q} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.075 fquenchconst=7
outputfile_q="${outputfile%.root}_q0p1.root"
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile_q} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.1 fquenchconst=10
outputfile_q="${outputfile%.root}_q0p125.root"
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile_q} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.125 fquenchconst=12
outputfile_q="${outputfile%.root}_q0p2.root"
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile_q} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.2 fquenchconst=14
outputfile_q="${outputfile%.root}_q0p4.root"
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile_q} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.4 fquenchconst=20

inputfile=${inputfile%.root}_nobias.root
outputfile="${inputfile%.root}_h.root"
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax} fquench=0.0 fquenchconst=0.0
