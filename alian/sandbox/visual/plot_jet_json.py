#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import plotly.graph_objects as go

def traverse(node, parent_coord=None, points=None, edges=None):
    if points is None:
        # list of tuples (phi, eta, pt, is_leaf, hover_text)
        points = []
    if edges is None:
        # list of tuples ((phi_p, eta_p, pt_p), (phi_c, eta_c, pt_c))
        edges = []
    
    # Extract coordinates from the current node
    phi = node.get("phi")
    eta = node.get("eta")
    pt = node.get("pt")
    current_coord = (phi, eta, pt)
    
    # Determine if this node is a leaf (i.e. no children in d1 or d2)
    has_child = node.get("d1") or node.get("d2")
    is_leaf = not bool(has_child)
    
    # Build hover text
    hover_text = f"Index: {node.get('index')}"
    hover_text += f"<br>pₜ,η,φ : ({pt:.2f}, {eta:.2f}, {phi:.2f})"
    if parent_coord is not None:
        parent_index = node.get("parent_index")
        hover_text += f"<br>Parent: {parent_index} ({parent_coord[2]:.2f}, {parent_coord[1]:.2f}, {parent_coord[0]:.2f})"
    if has_child:
        d1_index = node.get("d1").get("index")
        d2_index = node.get("d2").get("index")
        hover_text += f"<br>Daughters: present<br>D1: {d1_index}<br>D2: {d2_index}"
    else:
        hover_text += "<br>No daughters (leaf)"
    
    # Append the tuple including the hover text
    points.append((phi, eta, pt, is_leaf, hover_text))
    
    # If there is a parent coordinate, add an edge from parent to current
    if parent_coord is not None:
        edges.append((parent_coord, current_coord))
    
    # Recursively process children (keys "d1", "d2")
    for key in ["d1", "d2"]:
        child = node.get(key)
        if child:
            traverse(child, current_coord, points, edges)
    return points, edges

def main():
    # Load the JSON data. Adjust the path if necessary.
    with open("pythia_jets_lund.json", "r") as f:
        data = json.load(f)

    # Assuming the JSON has a single root with key "0"
    root = data["0"]

    # Traverse the JSON structure
    points, edges = traverse(root)

    # Prepare lists for plotting and collect hover text
    leaf_x, leaf_y, leaf_z, leaf_text = [], [], [], []
    nonleaf_x, nonleaf_y, nonleaf_z, nonleaf_text = [], [], [], []
    pts = []  # collect all pₜ for later axis range determination

    for phi, eta, pt, is_leaf, hover in points:
        pts.append(pt)
        if is_leaf:
            leaf_x.append(phi)
            leaf_y.append(eta)
            leaf_z.append(pt)
            leaf_text.append(hover)
        else:
            nonleaf_x.append(phi)
            nonleaf_y.append(eta)
            nonleaf_z.append(pt)
            nonleaf_text.append(hover)

    # Create the Plotly 3D scatter for the points
    fig = go.Figure()

    # Plot non-leaf nodes (splittings) in red with tooltips
    fig.add_trace(go.Scatter3d(
        x=nonleaf_x, y=nonleaf_y, z=nonleaf_z,
        mode='markers',
        marker=dict(size=5, color='red'),
        text=nonleaf_text,
        hoverinfo="text",
        name='splittings'
    ))

    # Plot leaf nodes (final state particles) in green with tooltips
    fig.add_trace(go.Scatter3d(
        x=leaf_x, y=leaf_y, z=leaf_z,
        mode='markers',
        marker=dict(size=5, color='green'),
        text=leaf_text,
        hoverinfo="text",
        name='final state particles'
    ))

    # Add connecting edges as individual line segments
    for parent, child in edges:
        fig.add_trace(go.Scatter3d(
            x=[parent[0], child[0]],
            y=[parent[1], child[1]],
            z=[parent[2], child[2]],
            mode='lines',
            line=dict(color='blue'),
            showlegend=False
        ))

    # Adjust the maximum of the pₜ (z) axis to be 1.1 * pt of index 0
    pt_root = root.get("pt")
    z_max = 1.1 * pt_root
    z_min = 0.0 # min(pts)

    fig.update_layout(
        scene=dict(
            xaxis_title="φ (rad)",
            yaxis_title="η",
            zaxis=dict(
                title="pₜ (GeV/c)",
                range=[z_min, z_max]
            )
        ),
        title="Pythia8 jet - declustered with C/A algorithm using FastJet",
        margin=dict(l=0, r=0, b=0, t=40)
    )

    # Save to HTML and display the plot
    fig.write_html("plot_jet_json.html", include_plotlyjs="cdn")
    fig.show()

if __name__ == "__main__":
    main()