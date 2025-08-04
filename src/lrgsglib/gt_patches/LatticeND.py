from .objects import *

class LatticeND(WeightedGraph):
    def __init__(self, shape: tuple, **kwargs):
        g = lattice(shape=shape, periodic=kwargs.pop('periodic', True))
        super().__init__(g=g, **kwargs)
