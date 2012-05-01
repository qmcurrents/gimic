import numpy as np
from copy import deepcopy
from grid import Grid, GridIterator
from gexceptions import NotReached

class Field(GridIterator):
    def __init__(self, grid):
        GridIterator.__init__(self)
        self.grid = grid
        self.npts = grid.size()

    def calc(self, func):
        raise NotReached('Baseclass method calc()')

    def __getitem__(self, idx):
        return self.getitem(idx)

class ScalarField(Field):
    def __init__(self, grid, dtype=float):
        Field.__init__(self, grid)
        self.field = np.ndarray(self.grid.size(), dtype=dtype)

    def getitem(self, idx):
        return self.field[idx[0], idx[1], idx[2]]

    itervalue = getitem

    def calc(self, func):
        for i, j, k in self.grid.range():
            r = self.grid.gridpoint((i, j, k))
            self.field[i, j, k] = func(r)

    def get(self, k=None):
        if k:
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

    def calc(self, func):
        for i, j, k in self.grid.range():
            r = self.grid.gridpoint((i, j, k))
            v = func(r)
            for n in range(self.dim):
                self.field[n][i, j, k] = v[n]

    def get(self, k=None):
        v = []
        if k:
            for i in range(self.dim):
                v[i] = self.field[i][:, :, k]
            return v
        else:
            return self.field

class TensorField(VectorField):
    dim = 9
    def __init__(self, grid, dtype=float):
        VectorField.__init__(self, grid, dtype)

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
        print i
    print f.get()

#    f = VectorField(grid)
#    f.calc(vtest)
#    for i in f:
#        print i
#    print f.get()

# vim:et:ts=4:sw=4
