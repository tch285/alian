#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import argparse
import os

def trajectory_current(fname, s_kt_or_kappa='kt', njets=0, 
                      y_limit_low=None, y_limit_high=None,
                      n_bins_x=None, n_bins_y=None, savefig=False,
                      fixed_njet_scale=0):
	"""
	Function to compute plot of the trajectory current from a given file.
	
	Parameters:
	fname (str): Path to the input file.
	s_kt_or_kappa (str): Specify 'kt' or 'kappa' for the analysis.
	njets (int): Number of jets to process.
	
	Returns:
	None
	"""
	# Load your parquet file
	df_jets = pd.read_parquet('pythia_lund_jet.parquet')

	# limit the number of jets based on the njets parameter
	if njets > 0:
		df_jets = df_jets[:njets]		

	# Extract lund records from all jets
	records = []
	for _, row in df_jets.iterrows():
			lunds = row['lunds']
			for i in range(len(lunds) - 1):
					kt_i, delta_i = lunds[i][s_kt_or_kappa], lunds[i]['delta']
					kt_j, delta_j = lunds[i+1][s_kt_or_kappa], lunds[i+1]['delta']
					x_i, y_i = np.log(1/delta_i), np.log(kt_i)
					x_j, y_j = np.log(1/delta_j), np.log(kt_j)
					records.append({'x': x_i, 'y': y_i, 'dx': x_j - x_i, 'dy': y_j - y_i})

	df = pd.DataFrame(records)

	# Define grid and bin the (x,y) positions
	nx, ny = 25, 25
	if n_bins_x is not None:
		nx = n_bins_x
	if n_bins_y is not None:
		ny = n_bins_y
	x_bins = np.linspace(df['x'].min(), df['x'].max(), nx + 1)
	y_bins = np.linspace(df['y'].min(), df['y'].max(), ny + 1)
	df['ix'] = pd.cut(df['x'], bins=x_bins, labels=False, include_lowest=True)
	df['iy'] = pd.cut(df['y'], bins=y_bins, labels=False, include_lowest=True)

	# Compute average displacement per bin
	group = (
			df.groupby(['ix', 'iy'])
				.agg(x=('x', 'mean'), y=('y', 'mean'),
						dx=('dx', 'mean'), dy=('dy', 'mean'))
				.dropna()
				.reset_index()
	)

	if y_limit_low is None:
		y_limit_low = df['y'].min()
		y_limit_low = -8
		if s_kt_or_kappa == 'kt':
				# Apply cut for plotting only (ln(kT) > -4)
				y_limit_low = -4

	if y_limit_high is None:
		y_limit_high = df['y'].max()

	# df_plot = df[df['y'] > y_limit_low]
	# group_plot = group[group['y'] > y_limit_low]

	df_plot = df
	group_plot = group

	plt.figure(figsize=(8, 6))

	# 1) density background (only for ln(kT) > -4)
	weights_plot = np.ones_like(df_plot['x']) / len(df_jets)
	if fixed_njet_scale > 0:
		weights_plot = np.ones_like(df_plot['x']) / fixed_njet_scale
	plt.hist2d(
			df_plot['x'], df_plot['y'],
			bins=[x_bins, y_bins],
			cmap='Blues',
			norm=mcolors.LogNorm(),
			weights=weights_plot
	)
	plt.clim(0.001, 1.0)
	plt.colorbar(label='Splitting count per jet (log scale)')

	# 2) trajectory current (quiver only for arrows starting at ln(kT) > -4)
	# Calculate scale depending on njets (scale increases from 1 to 5 for up to 1000 jets)
	max_scale = 5
	min_scale = 1
	scale_njets = 1000
	scale = min_scale + (max_scale - min_scale) * min(len(df_jets), scale_njets) / scale_njets
	plt.quiver(group_plot['x'], group_plot['y'], group_plot['dx'], group_plot['dy'], 
						angles='xy', scale_units='xy', scale=scale, width=0.003, color='red', alpha=1.0)

	plt.xlabel(r'$\ln(1/\Delta)$')
	ylabel_latex = r'$\ln(k_T)$' if s_kt_or_kappa == 'kt' else r'$\ln(\kappa)$'
	plt.ylabel(ylabel_latex)
	plt.ylim(y_limit_low, y_limit_high)
	plt.xlim(0.9, 5.5)
	plt.title(f'Splitting Density with Trajectory Current - {len(df_jets)} Jets')
	plt.grid(True, linestyle='--', linewidth=0.5)
	plt.tight_layout()
	if savefig:
		plt.savefig(f"trajectory_current_{njets}.png", dpi=300)
	else:
		plt.show()
 
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='just plot the lund plane', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--file', help='input file name', default='pythia_lund_jet.parquet', type=str)
	parser.add_argument('--kappa', help='use kappa instead of kt', action='store_true')
	parser.add_argument('--njets', help='number of jets to process', default=0, type=int)
	parser.add_argument('--y-low', help='y-axis lower limit', default=None, type=float)
	parser.add_argument('--y-high', help='y-axis upper limit', default=None, type=float)
	parser.add_argument('--n-bins-x', help='n bins on x == ln(1/Delta)', default=None, type=int)
	parser.add_argument('--n-bins-y', help='n bins on y == ln(kT|kappa)', default=None, type=int)
	parser.add_argument('--png', help='output file name', action='store_true', default=False)
	parser.add_argument('--njet-scale', help='normalize to a fixed number of jets', default=0, type=int)
	args = parser.parse_args()

	s_kt_or_kappa = 'kappa' if args.kappa else 'kt'
	trajectory_current('pythia_lund_jet.parquet', 
                    	s_kt_or_kappa='kt', njets=args.njets, 
                    	y_limit_high=args.y_high, y_limit_low=args.y_low,
											n_bins_x=args.n_bins_x, n_bins_y=args.n_bins_y, savefig=args.png, 
           						fixed_njet_scale=args.njet_scale)
 
 