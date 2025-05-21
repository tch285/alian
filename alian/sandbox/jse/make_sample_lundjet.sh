#!/bin/bash

nev=$1
if [ -z "$1" ]; then
	nev=10000
fi
# ./lund_jet.py --nev ${nev} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120

thermal_bg=0
if [ "$2" == "bg" ]; then
	thermal_bg=1
fi

function make_lund_jet() 
{
	thermal_bg_flag="$3"
	outdir="./sample_lundjet_${1}"
	mkdir -p ${outdir}
	# make the lund jet
	if [ "$2" == "0" ]; then
		./lund_jet.py ${thermal_bg_flag} --nev ${1} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120 --fixed-label ${2} --py-hardQCD --output "${outdir}/lund_jet_hardQCDany.parquet"
	elif [ "$2" == "1" ]; then
		./lund_jet.py ${thermal_bg_flag} --nev ${1} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120 --fixed-label ${2} --py-hardQCDquarks --output ${outdir}/lund_jet_hardQCDquarks.parquet
	elif [ "$2" == "2" ]; then
		./lund_jet.py ${thermal_bg_flag} --nev ${1} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120 --fixed-label ${2} --py-hardQCDgluons --output ${outdir}/lund_jet_hardQCDgluons.parquet
	elif [ "$2" == "3" ]; then
		./lund_jet.py ${thermal_bg_flag} --nev ${1} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120 --fixed-label ${2} --py-hardQCDlf --output ${outdir}/lund_jet_hardQCDlf.parquet
	elif [ "$2" == "4" ]; then
		./lund_jet.py ${thermal_bg_flag} --nev ${1} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120 --fixed-label ${2} --py-hardQCDcharm --output ${outdir}/lund_jet_hardQCDcharm.parquet
	elif [ "$2" == "5" ]; then
		./lund_jet.py ${thermal_bg_flag} --nev ${1} --py-pthatmin 100 --jet-pt-min 100 --jet-pt-max 120 --fixed-label ${2} --py-hardQCDbeauty --output ${outdir}/lund_jet_hardQCDbeauty.parquet
	fi
}
export -f make_lund_jet

if [ "$thermal_bg" == "1" ]; then
	thermal_bg_flag="--thermal"
	parallel --jobs 20 --progress make_lund_jet ::: ${nev} ::: 0 1 2 3 4 5 ::: ${thermal_bg_flag}
else
	thermal_bg_flag=""
	parallel --jobs 20 --progress make_lund_jet ::: ${nev} ::: 0 1 2 3 4 5
fi
