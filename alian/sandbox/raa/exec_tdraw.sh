#!/bin/bash

inputfile=$1
outputfile=$2
jetptmin=$3
jetptmax=$4

if [ -z "$inputfile" ]; then
	inputfile=pythia8_simple_eec_output.root
fi

if [ -z "$outputfile" ]; then
	outputfile=pythia_raa.root
fi

if [ -z "$jetptmin" ]; then
	jetptmin=20
fi

if [ -z "$jetptmax" ]; then
	jetptmax=100
fi

draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax}
inputfile=${inputfile%.root}_nobias.root
outputfile=${outputfile%.root}_nobias.root
draw_from_yaml.py -c tdraw_eec_eloss.yaml -d output=${outputfile} input=${inputfile} jetptmin=${jetptmin} jetptmax=${jetptmax}
