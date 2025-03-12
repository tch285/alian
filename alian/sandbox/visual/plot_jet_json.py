#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import plotly.graph_objects as go
import numpy as np

def create_cylinder(x, y, z_bottom, z_top, radius=0.05, resolution=20):
    """
    Create a vertical cylinder mesh (as a dict for go.Mesh3d) that extends from z_bottom to z_top.
    """
    theta = np.linspace(0, 2*np.pi, resolution, endpoint=False)
    # Bottom circle
    xb = x + radius * np.cos(theta)
    yb = y + radius * np.sin(theta)
    zb = np.full(theta.shape, z_bottom)
    # Top circle
    xt = x + radius * np.cos(theta)
    yt = y + radius * np.sin(theta)
    zt = np.full(theta.shape, z_top)
    
    # Combine vertices: bottom then top
    vertices_x = list(xb) + list(xt)
    vertices_y = list(yb) + list(yt)
    vertices_z = list(zb) + list(zt)
    n = resolution
    i_idx = []
    j_idx = []
    k_idx = []
    
    # Side faces (two triangles per segment)
    for k in range(n):
        next_k = (k+1) % n
        # Triangle 1: bottom_k, bottom_next, top_k
        i_idx.extend([k])
        j_idx.extend([next_k])
        k_idx.extend([k+n])
        # Triangle 2: top_k, bottom_next, top_next
        i_idx.extend([k+n])
        j_idx.extend([next_k])
        k_idx.extend([next_k+n])
    
    # Top cap: add center vertex and create triangles fan
    center_top_index = len(vertices_x)
    vertices_x.append(x)
    vertices_y.append(y)
    vertices_z.append(z_top)
    for k in range(n):
        next_k = (k+1) % n
        i_idx.extend([center_top_index])
        j_idx.extend([k+n])
        k_idx.extend([next_k+n])
    
    # Bottom cap: add center vertex and create triangles fan (reverse order for proper orientation)
    center_bottom_index = len(vertices_x)
    vertices_x.append(x)
    vertices_y.append(y)
    vertices_z.append(z_bottom)
    for k in range(n):
        next_k = (k+1) % n
        i_idx.extend([center_bottom_index])
        j_idx.extend([next_k])
        k_idx.extend([k])
    
    return dict(x=vertices_x, y=vertices_y, z=vertices_z, i=i_idx, j=j_idx, k=k_idx)

def traverse(node, parent_coord=None, points=None, edges=None):
    if points is None:
        # list of tuples (phi, eta, pt, is_leaf, hover_text)
        points = []
    if edges is None:
        # list of tuples ((phi, eta, pt), (phi, eta, pt), hover_text)
        edges = []
    
    # Extract coordinates from the current node
    phi = node.get("phi")
    eta = node.get("eta")
    pt = node.get("pt")
    current_coord = (phi, eta, pt)
    
    # Determine if this node is a leaf (i.e. no children in d1 or d2 exist)
    has_child = node.get("d1") or node.get("d2")
    is_leaf = not bool(has_child)
    
    # Build hover text with information about index, pₜ, parent, and daughters if available
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
        edges.append((parent_coord, current_coord, f"Edge from {parent_index} to {node.get('index')}"))
    
    # Recursively process children (keys "d1" and "d2")
    for key in ["d1", "d2"]:
        child = node.get(key)
        if child:
            traverse(child, current_coord, points, edges)
    return points, edges

def main():
    # Load the JSON data
    with open("pythia_jets_lund.json", "r") as f:
        data = json.load(f)

    # Assuming the JSON has a single root with key "0"
    root = data["0"]

    # Traverse the JSON structure to get points and edges
    points, edges = traverse(root)

    # Prepare lists for plotting
    leaf_x, leaf_y, leaf_z, leaf_text = [], [], [], []
    nonleaf_x, nonleaf_y, nonleaf_z, nonleaf_text = [], [], [], []
    pts = []  # to determine axis range

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

    # Create the 3D figure
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
    for parent, child, hover_text in edges:
        fig.add_trace(go.Scatter3d(
            x=[parent[0], child[0]],
            y=[parent[1], child[1]],
            z=[parent[2], child[2]],
            mode='lines',
            line=dict(color='blue'),
            showlegend=False,
            text=hover_text,
            hoverinfo="text"
        ))

    # Determine pₜ range for the z-axis
    pt_root = root.get("pt")
    z_max = 1.1 * pt_root
    z_min = 0.0  # or min(pts) if your data requires

    # Add vertical cylinders for each leaf node
    for phi, eta, pt, is_leaf, hover in points:
        if is_leaf:
            # cyl = create_cylinder(phi, eta, z_min, pt, radius=pt/100/20., resolution=20)
            cyl = create_cylinder(phi, eta, z_min, pt, radius=0.005, resolution=20)
            fig.add_trace(go.Mesh3d(
                x=cyl['x'], y=cyl['y'], z=cyl['z'],
                i=cyl['i'], j=cyl['j'], k=cyl['k'],
                color='green',
                opacity=pt/pt_root,
                flatshading=True,
                showscale=True,
                name="particle",
                # hoverinfo="skip", # disable hover for cylinders
                # hovertemplate=None,
                text=hover,
                hoverinfo="text",
            ))

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

    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title="φ (rad)",
                showspikes=False  # disable spike projections for x
            ),
            yaxis=dict(
                title="η",
                showspikes=False  # disable spike projections for y
            ),
            zaxis=dict(
                title="pₜ (GeV/c)",
                range=[z_min, z_max],
                showspikes=False  # disable spike projections for z
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