#!/bin/bash

# -*- coding: utf-8 -*-

wdir=$(dirname $(readlink -f $0))
wdir=$wdir/hf_b
echo $wdir

mkdir -p $wdir

cp hf_b.input $wdir/powheg.input

cd $wdir

# run the generator
pwhg_hvq

# run the hadronization
$HEPPYY_DEV/heppyy/example/powheg_pythia/pythia8_powheg.py pwgevents.lhe --nev 1000 --verbose

cd -