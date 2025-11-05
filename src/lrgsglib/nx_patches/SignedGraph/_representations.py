import networkx as nx
from typing import TYPE_CHECKING

from ...config.const import SG_REPR, SG_GRAPHINT_REPR

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def upd_graph_matrices(
    self: "SignedGraph",
    format: str = 'csr',
    on_g: str = SG_REPR
):
    # Store matrices in dictionaries keyed by graph representation
    self.adjacency_matrices[on_g] = self.get_adjacency_matrix(on_g=on_g)
    self.degree_matrices[on_g] = self.get_degree_matrix()
    self.signed_degree_matrices[on_g] = self.get_abs_degree_matrix()
    self.laplacian_matrices[on_g] = self.get_laplacian()
    self.signed_laplacian_matrices[on_g] = self.get_signed_laplacian()
    self.upd_Degree(on_g)


def upd_edge_sets(self: "SignedGraph", on_g: str = SG_REPR):
    self.eset[on_g] = set(self.gr[on_g].edges())
    self.fleset[on_g] = set([
        (u, v) 
        for u, v, data in self.gr[on_g].edges(data=True)
        if data.get('weight', 1) < 0
    ])
    self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])


def upd_Degree(self: "SignedGraph", on_g: str = SG_REPR):
    self.degrees = list(dict(self.gr[on_g].degree).values())


def zip_reprNodes(self: "SignedGraph", x, on_g: str = SG_REPR):
    return dict(zip(self.gr[on_g], self.gr[x]))


def zip_reprEdges(self: "SignedGraph", x, on_g: str = SG_REPR):
    return dict(zip(self.gr[on_g].edges(), self.gr[x].edges()))


def upd_NodeMap(self: "SignedGraph", on_g: str = SG_REPR):
    self.map_node[on_g] = {x: 
        {v: k for k, v in self.zip_reprNodes(x, on_g).items()} 
        for x in self.graph_reprs if x != on_g}


def upd_EdgeMap(self: "SignedGraph", on_g: str = SG_REPR):
    graph = self.gr[on_g]
    egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
    self.map_edge[on_g] = {x: {v: k for k, v in egraph_to(x).items()} 
                               for x in self.graph_reprs if x != on_g}


def upd_GraphRelabel(
        self: "SignedGraph",
        to_graph: str = SG_GRAPHINT_REPR, 
        on_g: str = SG_REPR
) -> None:
    node_map = self.map_node[to_graph][on_g]
    self.gr[to_graph] = nx.relabel_nodes(self.gr[on_g], node_map)
    self.eset[to_graph] = {x for e in self.eset[on_g] 
                        if (x := self.get_edge_mapping(e, to_graph, on_g))}
    self.fleset[to_graph] = {x for e in self.fleset[on_g]
                        if (x := self.get_edge_mapping(e, to_graph, on_g))}
    self.lfeset[to_graph] = {x for e in self.lfeset[on_g]
                        if (x := self.get_edge_mapping(e, to_graph, on_g))}
    setattr(self, to_graph, self.gr[to_graph])


def upd_ReprMaps(self: "SignedGraph", on_g: str = SG_REPR):
    graph = self.gr[on_g]
    ngraph_to = lambda x: dict(zip(graph, self.gr[x]))
    self.map_node[on_g] = {x: {v: k for k, v in ngraph_to(x).items()} 
                               for x in self.graph_reprs if x != on_g}
    egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
    self.map_edge[on_g] = {x: {v: k for k, v in egraph_to(x).items()} 
                               for x in self.graph_reprs if x != on_g}


def upd_GraphRepr_All(
    self: "SignedGraph",
    on_g: str = SG_REPR,
    also_itself: bool = True
):
    if also_itself:
        self.upd_ReprMaps(on_g=on_g)
    for i in self.graph_reprs:
        if i != on_g:
            self.upd_ReprMaps(on_g=i)
            self.upd_GraphRelabel(to_graph=i, on_g=on_g)
