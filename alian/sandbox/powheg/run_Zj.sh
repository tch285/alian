#!/bin/bash

# -*- coding: utf-8 -*-

wdir=$(dirname $(readlink -f $0))
wdir=$wdir/Zj
echo $wdir

mkdir -p $wdir

cp Zj.input $wdir/powheg.input

cd $wdir

# run the generator
# pwhg_Zj
# run the generator with explicit stages
pwhg_Zj -xgrid 1     # stage 1: grid generation
pwhg_Zj -xgrid 2     # stage 2: grid refinement
pwhg_Zj -xgrid 3     # stage 3: upper bound optimization
pwhg_Zj -xgrid 4     # stage 4: event generation

# run the hadronization
if [ -e "pwgevents.lhe" ]; then
	$HEPPYY_DEV/heppyy/example/powheg_pythia/pythia8_powheg.py pwgevents.lhe --nev 1000 --verbose
fi
cd -