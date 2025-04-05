#!/bin/bash

# -*- coding: utf-8 -*-

this_dir=$(dirname $(readlink -f $0))
wdir=$this_dir/Zj
echo $wdir

mkdir -p $wdir

cp Zj.input $wdir/powheg.input

if [ ! -e "pwgseeds.dat" ]; then
		echo "[w] pwgseeds.dat not found, creating a new one"
		./gen_seeds.py
	if [ ! -e "pwgseeds.dat" ]; then
		echo "[e] pwgseeds.dat not found, exiting"
		exit 1
	fi
fi
cp pwgseeds.dat $wdir/pwgseeds.dat

cd $wdir

# Function to run POWHEG with a specific seed
run_powheg() {
    seed=$1
    mkdir -p run_$seed
    cp powheg.input run_$seed/
		cp pwgseeds.dat run_$seed/
    cd run_$seed
    # pwhg_Zj -xgrid 1     # stage 1: grid generation
    # pwhg_Zj -xgrid 2     # stage 2: grid refinement
    # pwhg_Zj -xgrid 3     # stage 3: upper bound optimization
    # pwhg_Zj -xgrid 4     # stage 4: event generation
    echo ${seed} | pwhg_Zj
    mv pwgevents.lhe ../pwgevents-${seed}.lhe
    cd ..
    # rm -rf run_$seed
}

export -f run_powheg

ncpu=$(nproc --all)
njobs=$((ncpu - 1))
echo "Number of jobs: $njobs"

# Run multiple seeds in parallel
parallel -j ${njobs} --progress run_powheg ::: {1..10}  # Adjust -j for the number of parallel jobs

cd -