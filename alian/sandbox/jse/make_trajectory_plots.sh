#!/bin/bash

njets=1000

function make_trajectory() 
{
	fname="pythia_lund_jet.parquet"
	njets=1000
	# make the trajectory
	./trajectory_current.py --kappa --y-high 3.5 --njets $1 -f ${fname} --png --njet-scale ${njets}
}
export -f make_trajectory

parallel --jobs 20 --progress make_trajectory ::: $(seq 1 ${njets})

rm -rf ./plots
mkdir -p ./plots
mv trajectory_current_*.png ./plots/

# make a gif from the png files
# convert -delay 10 -loop 0 plots/trajectory_current_*.png plots/trajectory_current.gif

#make a movie from the png files
ffmpeg -framerate 30 -i plots/trajectory_current_%d.png -c:v libx264 -pix_fmt yuv420p plots/trajectory_current.mp4
