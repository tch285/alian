#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import plotly.graph_objects as go
import numpy as np

global args

def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def bezier_curve_constant_pull(p0, p2, num_points=20, pull=0.5):
    """
    Generate a quadratic Bézier curve between p0 and p2.
    The control point p1 is set to be the midpoint scaled by 'pull'.
    p0 and p2 are 3-tuples.
    """
    # Compute midpoint:
    mid = ((p0[0] + p2[0]) / 2, (p0[1] + p2[1]) / 2, (p0[2] + p2[2]) / 2)
    # Pull the midpoint toward the origin:
    p1 = (mid[0] * pull, mid[1] * pull, mid[2] * pull)
    xs, ys, zs = [], [], []
    for t in [i/(num_points-1) for i in range(num_points)]:
        x = (1-t)**2 * p0[0] + 2*(1-t)*t * p1[0] + t**2 * p2[0]
        y = (1-t)**2 * p0[1] + 2*(1-t)*t * p1[1] + t**2 * p2[1]
        z = (1-t)**2 * p0[2] + 2*(1-t)*t * p1[2] + t**2 * p2[2]
        xs.append(x)
        ys.append(y)
        zs.append(z)
    return xs, ys, zs

import math


def bezier_curve(p0, p2, num_points=20, pull=0.5, threshold=1.0):
    """
    Generate a quadratic Bézier curve between p0 and p2.
    The control point p1 is set to be the midpoint scaled by an adjusted pull.
    For particles where the distance between p0 and p2 is less than 'threshold',
    the pull is softened (raised closer to 1, meaning less pull toward the origin).
    
    p0 and p2 are 3-tuples.
    """
    # Compute Euclidean distance between p0 and p2
    dist = math.sqrt((p0[0]-p2[0])**2 + (p0[1]-p2[1])**2 + (p0[2]-p2[2])**2)
    
    # If points are close, adjust pull toward 1 (no pull)
    if dist < threshold:
        # adjusted_pull goes linearly from 1 (when dist==0) to pull (when dist==threshold)
        adjusted_pull = pull + (1 - pull) * (1 - dist/threshold)
    else:
        adjusted_pull = pull

    # Compute midpoint and then pull control point toward origin by adjusted_pull factor:
    mid = ((p0[0] + p2[0]) / 2, (p0[1] + p2[1]) / 2, (p0[2] + p2[2]) / 2)
    p1 = (mid[0] * adjusted_pull, mid[1] * adjusted_pull, mid[2] * adjusted_pull)
    
    xs, ys, zs = [], [], []
    for t in [i/(num_points-1) for i in range(num_points)]:
        x = (1-t)**2 * p0[0] + 2*(1-t)*t * p1[0] + t**2 * p2[0]
        y = (1-t)**2 * p0[1] + 2*(1-t)*t * p1[1] + t**2 * p2[1]
        z = (1-t)**2 * p0[2] + 2*(1-t)*t * p1[2] + t**2 * p2[2]
        xs.append(x)
        ys.append(y)
        zs.append(z)
    return xs, ys, zs

def log_rescale_x(x, max_x):
    if max_x == 0:
        return x
    min_accept = 1e-6
    if np.abs(x) < min_accept:
        return x
    else:
        # return np.log(np.abs(x / max_x)) * np.sign(x)
        return np.log(np.abs(x)) / np.log(max_x) * np.sign(x)

def rescale_x(x, max_x):
    return x / max_x

def rescale(px, py, pz, max_E):
    if args.rescale or args.log:
        if args.log:
            return log_rescale_x(px, max_E), log_rescale_x(py, max_E), log_rescale_x(pz, max_E)
        return rescale_x(px, max_E), rescale_x(py, max_E), rescale_x(pz, max_E)
    return px, py, pz

def plot_event(args):

    tree = load_json(args.file)

    nodes_x = []
    nodes_y = []  # will store pz
    nodes_z = []  # will store py
    labels = []
    nodes_color = []  # use RGBA strings for alpha control
    nodes_size = []
    edges_x = []
    edges_y = []
    edges_z = []
    edges_labels = []

    max_E = np.max([particle['momentum'][3] for particle in tree.values()])
    print('max_E:', max_E)
    
    # Create nodes
    for particle_id, particle in tree.items():
        # Unpack momentum; swap py and pz for plotting
        px, py, pz, _ = particle['momentum']
        px, py, pz = rescale(px, py, pz, max_E)
        nodes_x.append(px)
        nodes_y.append(pz)  # plot pz along y-axis
        nodes_z.append(py)  # plot py along z-axis
        labels.append(f"ID: {particle_id}, PDG: {particle['pdgId']}, Status: {particle['status']}, Name: {particle['name']}, Energy: {particle['momentum'][3]} isHardon: {particle['isHadron']}")
        
        # Set color based on whether the particle is final:
        # Final particles: opaque green; Non-final: semi-transparent red.
        if (particle_id == "1") or (particle_id == "2"):
            nodes_color.append("rgba(255,0,0,0.5)")
            nodes_size.append(10)
        else:
            if particle.get('isFinal', False):
                nodes_color.append("rgba(0,128,0,1)")
                nodes_size.append(10)
            else:
                nodes_color.append("rgba(0,0,128,0.5)")
                nodes_size.append(5)

        # Create edges (connections between mothers and daughters) with curved lines
        for daughter_id, daughter in particle['daughters'].items():
            # For plotting, we use: 
            # x: px, y: pz, z: py
            px, py, pz, _ = particle['momentum']
            px, py, pz = rescale(px, py, pz, max_E)
            daughter_px, daughter_py, daughter_pz, _ = daughter['momentum']
            daughter_px, daughter_py, daughter_pz = rescale(daughter_px, daughter_py, daughter_pz, max_E)
            p0 = (px, pz, py)
            p2 = (daughter_px, daughter_pz, daughter_py)
            xs, ys, zs = bezier_curve(p0, p2, num_points=20, pull=0.7, threshold=0.0)
            edges_x.extend(xs + [None])
            edges_y.extend(ys + [None])
            edges_z.extend(zs + [None])
            edges_labels.append(f"Mother ID: {particle_id} {particle['name']} -> Daughter ID: {daughter_id} {daughter['name']}")
                    
        # # Create edges (connections between mothers and daughters)
        # for daughter_id, daughter in particle['daughters'].items():
        #     daughter_px, daughter_py, daughter_pz, _ = daughter['momentum']
        #     # Start of edge (mother)
        #     edges_x.append(px)
        #     edges_y.append(pz)  # use pz
        #     edges_z.append(py)  # use py
        #     # End of edge (daughter)
        #     edges_x.append(daughter_px)
        #     edges_y.append(daughter_pz)
        #     edges_z.append(daughter_py)
        #     # Break between edges
        #     edges_x.append(None)
        #     edges_y.append(None)
        #     edges_z.append(None)
    
    # Create edge trace with alpha-blended edges
    edge_trace = go.Scatter3d(
        x=edges_x, y=edges_y, z=edges_z,
        mode='lines',
        line=dict(color='blue', width=5),
        opacity=0.1,  # overall edge transparency
        hoverinfo='text',
        text=edges_labels,
        name='mother-daughter links'
    )
    
    # Create node trace with variable colors (using the RGBA values)
    node_trace = go.Scatter3d(
        x=nodes_x, y=nodes_y, z=nodes_z,
        mode='markers',
        marker=dict(size=nodes_size, color=nodes_color),
        text=labels,
        hoverinfo='text',
        name='particles'
    )
    
    x_axis_title = "px (GeV/c)"
    y_axis_title = "pz (GeV/c)"
    z_axis_title = "py (GeV/c)"
    if args.log:
        x_axis_title = "log(px)/log(Ebeam)"
        y_axis_title = "log(pz)/log(Ebeam)"
        z_axis_title = "log(py)/log(Ebeam)"
    elif args.rescale:
        x_axis_title = "px / Ebeam"
        y_axis_title = "pz / Ebeam"
        z_axis_title = "py / Ebeam"
    # Create and display the figure
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
            title=dict(
                # text="Pythia8 Event Tree",
                x=0.5,
                y=0.98,         # Adjust title vertical position lower
                xanchor="center",
                yanchor="top"
            ),
            legend=dict(
                orientation="h",
                xanchor="left",
                yanchor="top"
            ),        
            scene=dict(
            xaxis=dict(
                title=x_axis_title,
                showspikes=False
            ),
            yaxis=dict(
                title=y_axis_title,
                showspikes=False
            ),
            zaxis=dict(
                title=z_axis_title,
                showspikes=False
            )
        ),
        margin=dict(l=0, r=0, b=0, t=0)
    )    
    foutname = 'pythia_jets_json2plotly_bezier.html'
    if args.log:
        foutname = foutname.replace('.html', '_log.html')
    if args.rescale:
        foutname = foutname.replace('.html', '_rescaled.html')
    fig.write_html(foutname)
    fig.show()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot particle event tree from JSON file.")
    parser.add_argument('file', type=str, help="Path to JSON file")
    parser.add_argument('--pull', type=float, default=0.5, help="Pull factor for Bézier curves")
    parser.add_argument('--threshold', type=float, default=1.0, help="Distance threshold for pull adjustment")
    parser.add_argument('--log', action='store_true', help="Use log scaling for momentum")
    parser.add_argument('--rescale', action='store_true', help="Use beam E scaling for momentum")
    
    args = parser.parse_args()
    
    # Load JSON and plot event
    plot_event(args)