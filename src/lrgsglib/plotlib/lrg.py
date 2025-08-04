import networkx as nx

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, fcluster

from ..nx_patches import get_giant_component
from ..utils.lrg import entropy, compute_laplacian_properties,\
    compute_normalized_linkage, compute_optimal_threshold
#
from .colorbars import imshow_colorbar_caxdivider
from .const_plotlib import *

__all__ = [
    'plot_graph_analysis',
]

def plot_graph_analysis(
    GG,
    ch_int_name_map,
    *,
    t1=-3,
    t2=5,
    steps=400,
    linkage_method='ward',
    entropy_func=entropy,
    laplacian_func=compute_laplacian_properties,
    linkage_func=compute_normalized_linkage,
):
    """
    Parameters
    ----------
    GG : networkx.Graph
        Weighted graph.
    ch_int_name_map : dict
        Node → label mapping used throughout plots.
    Returns
    -------
    fig : matplotlib.figure.Figure
    optimal_clusters : ndarray
        Cluster assignment for each node.
    """
    # adjacency matrix
    CC = nx.to_numpy_array(GG, weight='weight')
    GGtmp = GG.copy()
    GG= get_giant_component(GGtmp)
    # network entropy
    net_ent = entropy_func(GG, t1=t1, t2=t2, steps=steps)
    tau_scale = net_ent[-1]
    speC = net_ent[1] / net_ent[1].max()
    Sm1 = net_ent[0] / net_ent[0].max()

    # laplacian-based distances
    spectrum, L, rho, Trho, tau = laplacian_func(GG, tau=None)
    dists = squareform(Trho)
    linkage_matrix, label_list, _ = linkage_func(dists, GG, method=linkage_method)
    FlatClusteringTh, *_ = compute_optimal_threshold(linkage_matrix)
    optimal_clusters = fcluster(linkage_matrix, t=FlatClusteringTh, criterion='distance')

    # figure layout
    fig = plt.figure(constrained_layout=True, figsize=(9, 9))
    ax = fig.subplot_mosaic(
        """
        DDDDE
        DDDDE
        DDDDE
        BBCCE
        BBCCE
        """
    )

    # adjacency matrix
    im = ax['B'].imshow(CC, cmap='coolwarm_r', interpolation='none')
    imshow_colorbar_caxdivider(im, ax['B'], position='bottom', orientation='horizontal')
    ax['B'].axis('off')

    # entropy curves
    ax['C'].plot(tau_scale, Sm1, label=r'$1-S$')
    ax['C'].plot(tau_scale[:-1], speC, label=r'$C$')
    ax['C'].set_xscale('log')
    ax['C'].legend()
    ax['C'].set_xlabel(r'$\tau$')
    ax['C'].axvline(tau, color='r', linestyle='--')

    # dendrogram
    relabel_list = [ch_int_name_map[n] for n in label_list]
    dendro = dendrogram(
        linkage_matrix,
        ax=ax['E'],
        color_threshold=FlatClusteringTh,
        above_threshold_color='k',
        leaf_font_size=9,
        labels=relabel_list,
        orientation='right',
    )
    tmin = linkage_matrix[:, 2][0] * 0.8
    tmax = linkage_matrix[:, 2][-1] * 1.1
    ax['E'].axvline(FlatClusteringTh, color='b', linestyle='--', label=r'$\mathcal{D}_{\rm th}$')
    ax['E'].set_xscale('log')
    ax['E'].set_xlim(tmin, tmax)
    ax['E'].legend()
    ax['E'].set_xlabel(r'$\mathcal{D}/\mathcal{D}_{max}$')
    ax['E'].set_ylabel('Pin')

    # graph with cluster colors
    leaf_label_colors = {
        lbl: col for lbl, col in zip(dendro['ivl'], dendro['leaves_color_list'])
    }
    node_colors = [leaf_label_colors[ch_int_name_map[n]] if n in GG.nodes() else 'k' for n in GGtmp.nodes()]
    widths = [GGtmp[u][v].get('weight', 1.0) for u, v in GGtmp.edges()]
    edge_colors = [leaf_label_colors[ch_int_name_map[u]] if u in GG.nodes() else 'k' for u, v in GGtmp.edges()]

    pos = nx.spring_layout(GGtmp, seed=5)
    nx.draw(
        GGtmp,
        pos=pos,
        ax=ax['D'],
        node_size=200,
        font_size=10,
        width=widths,
        edge_color=edge_colors,
        node_color=node_colors,
        alpha=0.7,
        with_labels=True,
        labels=ch_int_name_map,
    )

    return fig, optimal_clusters