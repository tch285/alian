#!/usr/bin/env python3
"""
Create seaborn pair plot for lunds variables from parquet file
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json
from pathlib import Path

def load_lunds_data(parquet_file):
    """
    Load lunds data from parquet file and extract variables for pair plot
    """
    # Read the parquet file
    df = pd.read_parquet(parquet_file)
    
    # Extract all lunds data into a flat DataFrame
    all_lunds = []
    
    for idx, row in df.iterrows():
        lunds = row['lunds']
        for lund in lunds:
            # Add jet-level information to each lund entry
            lund_entry = lund.copy()
            lund_entry['jet_pt'] = row['pt']
            lund_entry['jet_eta'] = row['eta']
            lund_entry['jet_phi'] = row['phi']
            lund_entry['jet_m'] = row['m']
            lund_entry['label'] = row['label']
            all_lunds.append(lund_entry)
    
    return pd.DataFrame(all_lunds)

def create_lunds_pairplot(lunds_df, output_file='lunds_pairplot.png'):
    """
    Create seaborn pair plot for lunds variables
    """
    # Select key variables for the pair plot
    # These are the most physically meaningful variables from the lunds splitting
    key_vars = ['delta', 'eta', 'kappa', 'kt', 'm', 'pt', 'z', 'psi']
    
    # Filter the dataframe to only include these variables
    plot_df = lunds_df[key_vars + ['label']].copy()
    
    # Remove any infinite or NaN values
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()
    
    # Apply log scale to some variables that have wide dynamic range
    plot_df['log_kt'] = np.log10(plot_df['kt'] + 1e-10)  # Add small constant to avoid log(0)
    plot_df['log_kappa'] = np.log10(plot_df['kappa'] + 1e-10)
    
    # Select final variables for plotting
    final_vars = ['delta', 'eta', 'log_kappa', 'log_kt', 'm', 'pt', 'z']
    
    # Create the pair plot
    plt.figure(figsize=(12, 10))
    
    # Use hue to color by label if available
    if 'label' in plot_df.columns and plot_df['label'].nunique() > 1:
        g = sns.pairplot(plot_df[final_vars + ['label']], 
                        hue='label',
                        diag_kind='hist',
                        plot_kws={'alpha': 0.6, 's': 2},
                        diag_kws={'alpha': 0.7})
    else:
        g = sns.pairplot(plot_df[final_vars], 
                        diag_kind='hist',
                        plot_kws={'alpha': 0.6, 's': 2},
                        diag_kws={'alpha': 0.7})
    
    # Customize the plot
    g.fig.suptitle('Lund Plane Variables Pair Plot', y=1.02, fontsize=16)
    
    # Add variable descriptions as axis labels
    var_labels = {
        'delta': 'Δ (opening angle)',
        'eta': 'η (pseudorapidity)', 
        'log_kappa': 'log₁₀(κ)',
        'log_kt': 'log₁₀(kₜ)',
        'm': 'mass',
        'pt': 'pₜ',
        'z': 'z (momentum fraction)'
    }
    
    # Update axis labels
    for i, ax_row in enumerate(g.axes):
        for j, ax in enumerate(ax_row):
            if i == len(g.axes) - 1:  # Bottom row
                ax.set_xlabel(var_labels.get(final_vars[j], final_vars[j]))
            if j == 0:  # Left column
                ax.set_ylabel(var_labels.get(final_vars[i], final_vars[i]))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Pair plot saved to: {output_file}")
    
    return g

def print_data_summary(lunds_df):
    """
    Print summary statistics of the lunds data
    """
    print("Lunds Data Summary:")
    print(f"Total number of lund splittings: {len(lunds_df)}")
    print(f"Number of jets: {lunds_df['jet_pt'].nunique() if 'jet_pt' in lunds_df.columns else 'Unknown'}")
    
    if 'label' in lunds_df.columns:
        print(f"Labels distribution:")
        print(lunds_df['label'].value_counts().sort_index())
    
    print(f"\nVariable ranges:")
    key_vars = ['delta', 'eta', 'kappa', 'kt', 'm', 'pt', 'z']
    for var in key_vars:
        if var in lunds_df.columns:
            print(f"  {var}: [{lunds_df[var].min():.3e}, {lunds_df[var].max():.3e}]")

def main():
    """
    Main function to create the pair plot
    """
    # Define the parquet file path
    parquet_file = "/Users/ploskon/devel/alian/alian/sandbox/jse/sample_lundjet_1000/lund_jet_hardQCDgluons.parquet"
    
    # Check if file exists
    if not Path(parquet_file).exists():
        print(f"Error: File {parquet_file} not found!")
        return
    
    print("Loading lunds data from parquet file...")
    lunds_df = load_lunds_data(parquet_file)
    
    print_data_summary(lunds_df)
    
    print("\nCreating pair plot...")
    g = create_lunds_pairplot(lunds_df)
    
    # Also create a correlation matrix plot
    key_vars = ['delta', 'eta', 'kappa', 'kt', 'm', 'pt', 'z']
    corr_data = lunds_df[key_vars].copy()
    corr_data['log_kt'] = np.log10(corr_data['kt'] + 1e-10)
    corr_data['log_kappa'] = np.log10(corr_data['kappa'] + 1e-10)
    
    plt.figure(figsize=(10, 8))
    correlation_matrix = corr_data[['delta', 'eta', 'log_kappa', 'log_kt', 'm', 'pt', 'z']].corr()
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, 
                square=True, fmt='.2f')
    plt.title('Lund Variables Correlation Matrix')
    plt.tight_layout()
    plt.savefig('lunds_correlation_matrix.png', dpi=300, bbox_inches='tight')
    print("Correlation matrix saved to: lunds_correlation_matrix.png")
    
    plt.show()

if __name__ == "__main__":
    main()
