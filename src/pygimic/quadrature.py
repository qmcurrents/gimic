import numpy as np
from grid import Grid
from gauss_legendre import GaussLegendre
from gexecptions import NotImplemented

class GaussLegendre:
    def __init__(self, order):
        raise NotImplemented(GaussLegendre)

class GridQuadrature:
    def __init__(self, grid, func):
        self.grid = grid
        self.func = func
        raise NotImplemented(GridQuadrature)


# vim:et:ts=4:sw=4
