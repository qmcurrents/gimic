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

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from streamplot import streamplot

class BasePlot:
    def __init__(self, field):
        self.field = field
        self.grid = field.get_grid()

    def scalar_plot(self, filename, k=0):
        raise NotImplemented

    def scalar_plot3d(self, fname=None, k=0):
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

    def scalar_plot(self, fname=None, k=0):
        if not isinstance(self.field, ScalarField):
            raise TypeError('Not a scalar field')
        fig = plt.figure()
#        ax = fig.gca(projection='2d')
        X = self.grid.get_axis()[0]
        Y = self.grid.get_axis()[1]
        X, Y = np.meshgrid(X, Y)
        Z = self.field.get(k)
        levels = np.arange(-1000, 0, 50.0)
        cset = plt.contourf(X, Y, Z, levels, cmap=cm.jet)

#        ax.zaxis.set_major_locator(LinearLocator(10))
#        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(cset, shrink=0.5, aspect=5)

        plt.show()

    def scalar_plot3d(self, fname=None, k=0):
        if not isinstance(self.field, ScalarField):
            raise TypeError('Not a scalar field')
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X = self.grid.get_axis()[0]
        Y = self.grid.get_axis()[1]
        X, Y = np.meshgrid(X, Y)
        Z = self.field.get(k)
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                        linewidth=0, antialiased=False)
        #ax.set_zlim(-1.01, 1.01)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()

    def vector_plot(self, fname=None, k=0):
        if not isinstance(self.field, VectorField):
            raise TypeError('Not a vector field')
        fig = plt.figure()
        x = self.grid.get_axis()[0]
        y = self.grid.get_axis()[1]
        x, y = np.meshgrid(x, y)
        u, v, w = self.field.get(k)
        q = plt.quiver(v, u)
#        qk = plt.quiverkey(q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
#                               fontproperties={'weight': 'bold'})
        l,r,b,t = plt.axis()
        dx, dy = r-l, t-b
        plt.axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])

        plt.title('Minimal vector plot')
        plt.show()

    def stream_plot(self, fname=None, k=0):
        if not isinstance(self.field, VectorField):
            raise TypeError('Not a vector field')
        fig = plt.figure()
        x = self.grid.get_axis()[0]
        y = self.grid.get_axis()[1]
#        x, y = np.meshgrid(x, y)
        u, v, w = self.field.get(k)
#        streamplot(x, y, u, v, density=1, INTEGRATOR='RK4', color='b')
        streamplot(x, y, v, u, density=1, INTEGRATOR='RK4', color=u)
        plt.show()


class NglPlot(BasePlot):
    def __init__(self, field):
        BasePlot.__init__(self, field)

if __name__ == '__main__':
    s=-1.0

    def sfunc(r):
        global s
        s += 1.0 
        return s

    def vfunc(r):
        return math.cos(r[0]), math.sin(r[0]), 0.0
      

    grd = Grid(l=(6, 6, 0),
            origin=(-3.0, -3.0, 0.0),
            npts=30)
#    f = ScalarField(grd)
#    f.calc(sfunc)
#    p = MatPlot(f)
#    p.scalar_plot('foo')

    f = VectorField(grd)
    f.calc(vfunc)
    p = MatPlot(f)
    p.vector_plot('foo')
    p.stream_plot('foo')

# vim:et:ts=4:sw=4
