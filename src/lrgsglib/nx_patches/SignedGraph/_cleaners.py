import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def remove_ising_clust_files(self: "SignedGraph"):
    for pname in self.clPname:
        os.remove(pname)

def remove_edgl_file(self: "SignedGraph"):
    os.remove(self.path_exp_edgl)

def remove_eigV_file(self: "SignedGraph"):
    os.remove(self.eigVPname)

def remove_exported_files(self: "SignedGraph"):
    if hasattr(self, "path_exp_edgl"):
        self.remove_edgl_file()
    if hasattr(self, "eigVPname"):
        self.remove_eigV_file()
    if hasattr(self, "clPname"):
        self.remove_ising_clust_files()

def clean_gclutil(self: "SignedGraph", k, val, on_g: str = 'graph'):
    self.graph_clustering_utility[k][val][on_g] = None

