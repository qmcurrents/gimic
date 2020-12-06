#
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
import math
from .grid import Grid
from .field import VectorField, ScalarField

class CurrentField(VectorField):
    def __init__(self, grid, gimic):
        VectorField.__init__(self, grid)
        self.calc(gimic.jvector)

    def get_normal_flow(self, comp='total'):
        nf = ScalarField(self.grid)
        normal = self.grid.get_normal()
        for i, j, k in self.range():
            x = np.dot(self[i, j, k], normal)
            if comp == '+':
                if x > 0.0:
                    nf[i, j, k] = x
                else:
                    nf[i, j, k] = 0.0
            elif comp == '-':
                if x < 0.0:
                    nf[i, j, k] = x
                else:
                    nf[i, j, k] = 0.0
            else:
                nf[i, j, k] = x
        return nf

    def get_modulus(self, comp='total'):
        nf = ScalarField(self.grid)
        normal = self.grid.get_normal()
        for i, j, k in self.range():
            v = self[i, j, k]
            x = np.dot(v, normal)
            if comp == '+':
                if x > 0.0:
                    nf[i, j, k] = np.linalg.norm(v)
                else:
                    nf[i, j, k] = 0.0
            elif comp == '-':
                if x < 0.0:
                    nf[i, j, k] = np.linalg.norm(v)
                else:
                    nf[i, j, k] = 0.0
            else:
                nf[i, j, k] = np.linalg.norm(v)
        return nf

if __name__ == '__main__':
    class FakeGimic:
        def jvector(self, r):
            return np.array(r)

    grid = Grid((2.0, 2.0, 2.0), npts=4)
    foo = FakeGimic()
    f = CurrentField(grid, foo)
    nf = f.get_normal_flow('+')
    for i in nf:
        print(i)
    nf = f.get_modulus()
    for i in nf:
        print(i)


# vim:et:ts=4:sw=4
