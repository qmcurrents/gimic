import numpy as np
from grid import Grid
from gauss_legendre import GaussLegendre

class GridQuadrature:
    def __init__(self, grid, func):
        self.grid = grid
        self.func = func

# vim:et:ts=4:sw=4
