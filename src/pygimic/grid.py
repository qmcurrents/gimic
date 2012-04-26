#
# Basic grid class
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np

class Grid:
    def __init__(self):
        self.basis = np.identity(3)
        self.l = np.zeros(3)
        self.origin = np.zeros(3)
        self.ortho = np.zeros(3)
        self.step = np.ones(3)
        self.distribution = 'even'
        self.radius = None
        self.grid = None

    def get_basis_vetors(self):
        return self.basis

    def get_i(self):
        return self.basis[:,0]

    def get_j(self):
        return self.basis[:,1]

    def get_k(self):
        return self.basis[:,2]

    def point(self):
        pass

    def realpoint(self):
        pass

    def normal(self):
        pass

    def normalize(self):
        pass

    def size(self):
        pass

    def lenght(self):
        pass

    def index(self):
        pass

    def get_center(self):
        pass

    def get_ortho(self):
        pass

    def get_range(self):
        pass

    def is3d(self):
        pass

    def check_orthogonal(self):
        pass

# vim:et:ts=4:sw=4
