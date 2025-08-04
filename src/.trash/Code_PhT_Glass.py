import os
import sys
import random
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
from numpy import linalg as LA


L=int(sys.argv[1])
p1=float(sys.argv[2])
navrg=1000
smax=np.zeros((navrg))
scnd=np.zeros((navrg))
trhd=np.zeros((navrg))
four=np.zeros((navrg))
weav=0.0
var=0.0
mu=1.0
file1='data/PhaseT_L_'+str(L)+'_sigma_'+str(p1)+'_'
cont=0
for avrg in range(navrg):
    H = nx.grid_graph(dim=(L, L, L),periodic=True)
    # Remove a fraction of edges
    fraction = 0.25  # Adjust this to change the fraction of edges to remove
    # Convert edge view to a list
    edges = list(H.edges())
    edges_to_remove = random.sample(edges, int(fraction * len(edges)))
    #Remove the selected edges
    H.remove_edges_from(edges_to_remove)
    # Relabel the nodes
    mapping = {old_label: new_label for new_label, old_label in enumerate(H.nodes())}
    G0 = nx.relabel_nodes(H, mapping)
    Gcc = sorted(nx.connected_components(G0), key=len, reverse=True)
    G = G0.subgraph(Gcc[0])
    weights = {edge: random.normalvariate(mu, p1) for edge in G.edges()}
    we=list(weights.values())
    we1=np.array(we)
    weav=weav+(len(we1[we1<0])/len(we1))
    nx.set_edge_attributes(G, values=weights, name='weight')
    N=len(G.nodes())
    # Calculate adjacency matrix
    adj = nx.adjacency_matrix(G, weight='weight',dtype=np.float64)
    # Compute the vector of degrees
    diabs = np.array(adj.sum(axis=1)).flatten()
    # Laplacian matrix
    slapl = diags(diabs, 0, format='csr', dtype=np.float64) - adj
    _, eigV1 = LA.eigh(slapl.todense())
    eigV=np.copy(eigV1.T[0])
    # Find the first eigenvector
    #_, eigV = eigsh(slapl, k=1, which='SM')
    # Count the number of negative and positive elements
    negative_count = np.sum(eigV < 0)
    positive_count = np.sum(eigV > 0)
    # Determine majority sign
    minority_sign = -1 if negative_count > positive_count else 1
    # Find positions of elements with majority sign
    if minority_sign == -1:
        min_pos = np.where(eigV < 0)[1]  # Get positions of negative values
    else:
        min_pos = np.where(eigV > 0)[1]  # Get positions of positive values

    # Extract subgraph with nodes at min_pos
    subG = G.subgraph(min_pos)
    # Find the connected components of the subgraph
    clust = np.array([len(component) for component in nx.connected_components(subG)])
    srtclust=sorted(clust)
    if len(clust) > 1:
        smax[cont] = np.max(clust) / (1.0 * N)
        scnd[cont] = srtclust[-2] / (1.0 * N)
        if len(clust) > 2:
            trhd[cont] = srtclust[-3] / (1.0 * N)
        if len(clust) > 3:
            four[cont] = srtclust[-4] / (1.0 * N)
        var += ((np.sum(np.multiply(clust,clust))- np.max(clust)*np.max(clust))/(np.sum(clust) - np.max(clust)))
    else:
        smax[cont] = 1.0*len(clust)
        var += 0.0
    if (avrg % 100) == 0:
        x = [np.mean(smax[0:cont]), np.std(smax[0:cont]), np.mean(scnd[0:cont]), np.std(scnd[0:cont]), np.mean(trhd[0:cont]), np.std(trhd[0:cont]),
         np.mean(four[0:cont]), np.std(four[0:cont]), var / (avrg + 1), p1,weav/(avrg+1.0), int(avrg + 1)]
        np.savetxt(file1, np.atleast_2d(x), delimiter=',')
    cont=cont+1
x = [np.mean(smax[0:cont]), np.std(smax[0:cont]), np.mean(scnd[0:cont]), np.std(scnd[0:cont]), np.mean(trhd[0:cont]), np.std(trhd[0:cont]),
         np.mean(four[0:cont]), np.std(four[0:cont]), var / (avrg + 1), p1, weav/(avrg+1.0), int(avrg + 1)]
np.savetxt(file1, np.atleast_2d(x), delimiter=',')
# Close the file
file1.close()
