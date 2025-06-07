#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import plotly.graph_objects as go

def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def plot_event(tree):
    nodes_x = []
    nodes_y = []  # will store pz
    nodes_z = []  # will store py
    labels = []
    nodes_color = []  # use RGBA strings for alpha control
    edges_x = []
    edges_y = []
    edges_z = []
    
    # Create nodes
    for particle_id, particle in tree.items():
        # Unpack momentum; swap py and pz for plotting
        px, py, pz, _ = particle['momentum']
        nodes_x.append(px)
        nodes_y.append(pz)  # plot pz along y-axis
        nodes_z.append(py)  # plot py along z-axis
        labels.append(f"ID: {particle_id}, PDG: {particle['pdgId']}, Status: {particle['status']}, Name: {particle['name']}, Energy: {particle['momentum'][3]}")
        
        # Set color based on whether the particle is final:
        # Final particles: opaque green; Non-final: semi-transparent red.
        if particle.get('isFinal', False):
            nodes_color.append("rgba(0,128,0,1)")
        else:
            nodes_color.append("rgba(255,0,0,0.5)")
        
        # Create edges (connections between mothers and daughters)
        for daughter_id, daughter in particle['daughters'].items():
            daughter_px, daughter_py, daughter_pz, _ = daughter['momentum']
            # Start of edge (mother)
            edges_x.append(px)
            edges_y.append(pz)  # use pz
            edges_z.append(py)  # use py
            # End of edge (daughter)
            edges_x.append(daughter_px)
            edges_y.append(daughter_pz)
            edges_z.append(daughter_py)
            # Break between edges
            edges_x.append(None)
            edges_y.append(None)
            edges_z.append(None)
    
    # Create edge trace with alpha-blended edges
    edge_trace = go.Scatter3d(
        x=edges_x, y=edges_y, z=edges_z,
        mode='lines',
        line=dict(color='blue', width=1),
        opacity=0.5,  # overall edge transparency
        hoverinfo='none'
    )
    
    # Create node trace with variable colors (using the RGBA values)
    node_trace = go.Scatter3d(
        x=nodes_x, y=nodes_y, z=nodes_z,
        mode='markers',
        marker=dict(size=5, color=nodes_color),
        text=labels,
        hoverinfo='text'
    )
    
    # Create and display the figure
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title="Pythia8 Event Tree",
        scene=dict(
            xaxis_title="px",
            yaxis_title="pz",
            zaxis_title="py"
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    fig.write_html("pythia_jets_json2plotly.html")
    fig.show()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Plot particle event tree from JSON file.")
    parser.add_argument('file', type=str, help="Path to JSON file")
    
    args = parser.parse_args()
    
    # Load JSON and plot event
    event_tree = load_json(args.file)
    plot_event(event_tree)