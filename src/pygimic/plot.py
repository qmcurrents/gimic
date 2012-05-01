#
# Plotting facilities
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
import math
from copy import deepcopy
from gimic_exceptions import NotImplemented
from grid import Grid
from field import ScalarField, VectorField

class BasePlot:
    def __init__(self, field):
        self.field = field

    def scalar_plot(self, filename, k=0):
        print filename
        raise NotImplemented

    def modulus_plot(self, filename, k=0):
        raise NotImplemented

    def vector_plot(self, filename, k=0):
        raise NotImplemented

    def stream_plot(self, filename, k=0):
        raise NotImplemented

    def contour_plot(self, filename, k=0):
        raise NotImplemented

class VtkPlot(BasePlot):
    def __init__(self, field):
        BasePlot.__init__(self, field)

class GnuPlot(BasePlot):
    def __init__(self, field):
        BasePlot.__init__(self, field)

class CubePlot(BasePlot):
    def __init__(self, field):
        BasePlot.__init__(self, field)

class MatPlot(BasePlot):
    def __init__(self, field):
        BasePlot.__init__(self, field)

class NglPlot(BasePlot):
    def __init__(self, field):
        BasePlot.__init__(self, field)

if __name__ == '__main__':
    s = 0.0
    def func(r):
        global s
        s += sum(r)
        return s

    grd = Grid(l=(6, 6, 0),
            origin=(-3.0, -3.0, 0.0),
            npts=30)
    f = ScalarField(grd)
    f.calc(func)
    plt = MatPlot(f)
    plt.scalar_plot('foo')


# vim:et:ts=4:sw=4
