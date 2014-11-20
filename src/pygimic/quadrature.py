#
# Simple quadrature on scalar fields
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
import math
from grid import Grid
from field import ScalarField

class FieldQuadrature:
    def integrate(self, field):
        npts = field.size()
        grid = field.get_grid()
        wgts_i = grid.get_weigths(0)
        wgts_j = grid.get_weigths(1)
        wgts_k = grid.get_weigths(2)

        x = np.zeros(npts[2])
        for k in range(npts[2]):
            v = field.get_field(k)
            w = np.dot(wgts_i, v)
            x[k] = np.dot(wgts_j, w)
        return np.dot(wgts_k, x)


if __name__ == '__main__':
    def gau(r):
        return math.exp(-1.0 * (r[0]**2 + r[1]**2 + r[2]**2))

    def fsin(r):
        return math.sin(r[0] + r[1] + r[2])

    q = FieldQuadrature()
    grid = Grid(l=(20.0, 20.0, 20.0), origin=(-10.0, -15.0, -15.0), npts=50, 
            distribution='gauss')
#    print grid.get_points(0)
    f = ScalarField(grid)
    f.calc(gau)
    x = q.integrate(f)
    print 'int(exp(-r^2) = ', x, math.sqrt(math.pi)**3

    grid = Grid(l=(2.0*math.pi, 0, 0), origin=(-math.pi, 0, 0), npts=40, 
            distribution='gauss')
#    print grid.get_points(0)
    f = ScalarField(grid)
    f.calc(fsin)
    x = q.integrate(f)
    print 'int(sin(r)) =', x 



# vim:et:ts=4:sw=4
