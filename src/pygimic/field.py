#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

import numpy as np
from copy import deepcopy
from .grid import Grid, GridIterator
from .gimic_exceptions import NotReached

class Field(GridIterator):
    def __init__(self, grid):
        GridIterator.__init__(self)
        self.grid = grid
        self.npts = grid.size()

    def calc(self, func):
        raise NotReached('Baseclass method calc()')

    def __getitem__(self, idx):
        return self.getitem(idx)

    def __setitem__(self, idx, val):
        self.setitem(idx, val)

    def get_grid(self):
        return self.grid

    def size(self):
        return self.grid.size()

class ScalarField(Field):
    def __init__(self, grid, dtype=float):
        Field.__init__(self, grid)
        self.field = np.ndarray(self.grid.size(), dtype=dtype)

    def getitem(self, idx):
        return self.field[idx[0], idx[1], idx[2]]

    itervalue = getitem

    def setitem(self, idx, val):
        self.field[idx[0], idx[1], idx[2]] = val

    def calc(self, func):
        for i, j, k in self.grid.range():
            r = self.grid.gridpoint((i, j, k))
            w = self.grid.gridweight((i, j, k))
            if abs(w) > 0.0:
                self.field[i, j, k] = func(r)
            else:
                self.field[i, j, k] = 0.0


    def get_field(self, k=None):
        if k is not None:
            return self.field[:, :, k]
        else:
            return self.field

class VectorField(Field):
    dim = 3
    def __init__(self, grid, dtype=float):
        Field.__init__(self, grid)
        self.field = []
        for i in range(self.dim):
            self.field.append(np.ndarray(self.grid.size(), dtype=dtype))

    def getitem(self, idx):
        v = np.zeros(self.dim)
        for i in range(self.dim):
            v[i] = self.field[i][idx[0], idx[1], idx[2]]
        return v

    itervalue = getitem

    def setitem(self, idx, val):
        i, j, k = idx
        self.field[0][i, j, k] = val[0]
        self.field[1][i, j, k] = val[1]
        self.field[2][i, j, k] = val[2]

    def calc(self, func):
        for i, j, k in self.grid.range():
            r = self.grid.gridpoint((i, j, k))
            w = self.grid.gridweight((i, j, k))
            if abs(w) > 0.0:
                v = func(r)
            else:
                v = np.zeros(self.dim)
            self.field[0][i, j, k] = v[0]
            self.field[1][i, j, k] = v[1]
            self.field[2][i, j, k] = v[2]

    def get_field(self, k=None):
        v = []
        if k is not None:
            for i in range(self.dim):
                v.append(self.field[i][:, :, k])
            return v
        else:
            return self.field

class TensorField(VectorField):
    dim = 9
    def __init__(self, grid, dtype=float):
        VectorField.__init__(self, grid, dtype)

    def setitem(self, idx, val):
        raise NotImplemented()

if __name__ == '__main__':
    s=-1.0
    v=np.zeros(3)
    t=np.zeros(9)

    def stest(r):
        global s
        s += 1.0
        return s

    def vtest(r):
        global v
        v += 1.0
        return v

    def ttest(r):
        global t
        t += 1.0
        return t

    grid = Grid((2.0, 2.0, 2.0), npts=4)
    f = ScalarField(grid)
    f.calc(stest)
    for i in f:
        print(i)
    print(f.get_field())

#    f = VectorField(grid)
#    f.calc(vtest)
#    for i in f:
#        print i
#    print f.get_field()

# vim:et:ts=4:sw=4
