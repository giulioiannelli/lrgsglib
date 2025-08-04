import numpy as np
import matplotlib.pyplot as plt; plt.ion()
import networkx as nx
import netgraph

# stochastic block model parameters
N1 = N2 = 10
sizes = [N1, N2]
probs = [[0.95, 0.05], [0.05, 0.95]]

G = nx.stochastic_block_model(sizes, probs, seed=1)

for edge in G.edges():
    G.add_edge(edge[0], edge[1], weight=1)
    sign = G.get_edge_data(edge[0], edge[1])

plt.figure(figsize=(8,3)) 
nx.draw(G, node_size=300, with_labels = True)

H = nx.quotient_graph(G, G.graph["partition"], relabel=False)

community1 = list(list(H.nodes())[0])
sub_community1 = community1[:len(community1)//2]
sub_community2 = community1[len(community1)//2:]
for vertex in sub_community1:
    for edge in G.edges(vertex):
        if edge[1] in sub_community2:
            G[edge[0]][edge[1]]['weight'] =-1

I = netgraph.InteractiveGraph(G, node_size=12)