import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX


def remove_ising_clust_files(self: "SignedGraphNX"):
    for pname in self.clPname:
        if os.path.exists(pname):
            try:
                os.remove(pname)
            except Exception as e:
                print(f"Warning: Could not remove Ising cluster file: {e}", file=sys.stderr)

def remove_edgl_file(self: "SignedGraphNX"):
    if os.path.exists(self.path_exp_edgl):
        try:
            os.remove(self.path_exp_edgl)
        except Exception as e:
            print(f"Warning: Could not remove edgelist file: {e}", file=sys.stderr)

def remove_eigV_file(self: "SignedGraphNX"):
    if os.path.exists(self.eigVPname):
        try:
            os.remove(self.eigVPname)
        except Exception as e:
            print(f"Warning: Could not remove eigenvector file: {e}", file=sys.stderr)

def remove_adj_file(self: "SignedGraphNX"):
    if os.path.exists(self.path_exp_adj):
        try:
            os.remove(self.path_exp_adj)
        except Exception as e:
            print(f"Warning: Could not remove adjacency file: {e}", file=sys.stderr)

def remove_exported_files(self: "SignedGraphNX"):
    if hasattr(self, "path_exp_edgl"):
        self.remove_edgl_file()
    if hasattr(self, "path_exp_adj"):
        self.remove_adj_file()
    if hasattr(self, "eigVPname"):
        self.remove_eigV_file()
    if hasattr(self, "clPname"):
        self.remove_ising_clust_files()

def clean_gclutil(self: "SignedGraphNX", k, val, on_g: str = 'graph'):
    self.graph_clustering_utility[k][val][on_g] = None

