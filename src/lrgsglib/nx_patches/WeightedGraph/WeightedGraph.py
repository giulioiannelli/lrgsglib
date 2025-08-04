from ..common import *


class WeightedGraph(Graph):
    def __init__(self):
        super().__init__()
        self.weights = {}
