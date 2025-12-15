import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def remove_ising_clust_files(self: "SignedGraph"):
    for pname in self.clPname:
        try:
            os.remove(pname)
        except FileNotFoundError as e:
            print(f"Warning: Could not remove Ising cluster file: {e}", file=sys.stderr)

def remove_edgl_file(self: "SignedGraph"):
    try:
        os.remove(self.path_exp_edgl)
    except FileNotFoundError as e:
        print(f"Warning: Could not remove edgelist file: {e}", file=sys.stderr)

def remove_eigV_file(self: "SignedGraph"):
    try:
        os.remove(self.eigVPname)
    except FileNotFoundError as e:
        print(f"Warning: Could not remove eigenvector file: {e}", file=sys.stderr)

def remove_adj_file(self: "SignedGraph"):
    try:
        os.remove(self.path_exp_adj)
    except FileNotFoundError as e:
        print(f"Warning: Could not remove adjacency file: {e}", file=sys.stderr)

def remove_exported_files(self: "SignedGraph"):
    if hasattr(self, "path_exp_edgl"):
        self.remove_edgl_file()
    if hasattr(self, "path_exp_adj"):
        self.remove_adj_file()
    if hasattr(self, "eigVPname"):
        self.remove_eigV_file()
    if hasattr(self, "clPname"):
        self.remove_ising_clust_files()

def clean_gclutil(self: "SignedGraph", k, val, on_g: str = 'graph'):
    self.graph_clustering_utility[k][val][on_g] = None

