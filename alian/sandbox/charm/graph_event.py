import networkx as nx
import matplotlib.pyplot as plt

def build_production_graph(idx, pythia, G=None):
    """
    Recursively builds a directed graph of the production tree for the particle at index `idx`.
    Each node is labeled with the particle's index, PDG id, and transverse momentum (pT).
    """
    if G is None:
        G = nx.DiGraph()
        
    particle = pythia.event[idx]
    node_label = f"{idx}\npid: {particle.id()}\npT: {particle.pT():.2f}"
    if idx not in G:
        G.add_node(idx, label=node_label)
        
    # Retrieve mothers (mother1 and mother2)
    mother1 = particle.mother1()
    mother2 = particle.mother2()
    
    for mother in set([mother1, mother2]):
        if mother != -1 and mother > 0:
            mother_particle = pythia.event[mother]
            m_label = f"{mother}\npid: {mother_particle.id()}\npT: {mother_particle.pT():.2f}"
            if mother not in G:
                G.add_node(mother, label=m_label)
            G.add_edge(mother, idx)
            build_production_graph(mother, pythia, G)
    return G

def show_production_graph(idx, pythia):
    """
    Builds and displays the production tree graph for the particle at index `idx`.
    """
    G = build_production_graph(idx, pythia)
    try:
        # Use Graphviz layout if available for a tree-like display.
        pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
    except Exception as e:
        print("Graphviz layout not available, using spring layout instead.", e)
        pos = nx.spring_layout(G)
        
    labels = nx.get_node_attributes(G, 'label')
    plt.figure(figsize=(10, 8))
    nx.draw(G, pos,
            with_labels=True,
            labels=labels,
            node_color='lightblue',
            node_size=2000,
            font_size=8,
            arrowstyle='->',
            arrowsize=10)
    plt.title("Pythia Production Tree")
    plt.show()

# Example usage:
# After generating an event and identifying a particular particle (e.g., a charm hadron),
# call show_production_graph with its event index, for example:
#
#    show_production_graph(charm_hadron_index, pythia)
#
# This will pop up a window with the graphical representation of the production tree.