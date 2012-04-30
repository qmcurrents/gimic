import numpy as np
from copy import deepcopy
from grid import Grid
from magnet import Magnet

class Field:
    def __init__(self, grid):
        self.grid = grid
        self.field = []
        self.it = 0
        self.iteridx = [0, 0, 0]

    def calc(self, func):
        for r in self.grid:
            self.field.append(func(r))

    def __getitem__(self, i, j, k=0):
        return self.get(i, j, k)

    def get(self, i, j, k=0):
        n = grid.size()
        return self.field[i + n[0] * j + n[0] * n[1] * k]

    def __iter__(self):
        return self

    def next(self):
        n = self.it
        if self.it < len(self.field):
            self.it += 1
        else:
            self.it = 0
            raise StopIteration()
        return self.field[n]

    def fieldgen(self):
        npts = self.grid.size()
        for i in self.field:
            cur = deepcopy(self.iteridx)
            self.iteridx[0] += 1
            if self.iteridx[0] == npts[0]:
                self.iteridx[0] = 0
                self.iteridx[1] += 1
                if self.iteridx[1] == npts[1]:
                    self.iteridx[0] = 0
                    self.iteridx[1] = 0
                    self.iteridx[2] += 1
                    if self.iteridx[2] == npts[2]:
                        self.iteridx = [0, 0, 0]
            yield self.grid.point(cur), self.get(cur[0], cur[1], cur[2])


if __name__ == '__main__':
    q=-1.0

    def ftest(r):
        global q
        q += 1.0 
        return q

    grid = Grid((2.0, 2.0, 2.0), npts=4)
    f = Field(grid)
    f.calc(ftest)
    for i in f.fieldgen():
        print i

# vim:et:ts=4:sw=4
