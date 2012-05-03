#
# Basic grid class
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
import math
from copy import deepcopy
#from magnet import Magnet
from atom import Atom
from gimic_exceptions import NotImplemented

class GridIterator:
    def __init__(self):
        self.iteridx = [0, 0, 0]
        self.npts = [0, 0, 0]

    def __iter__(self):
        return self

    def next(self):
        'Iterator to return data points in normal order'
        if self.iteridx is None:
            self.iteridx = [0, 0, 0]
            raise StopIteration()
        cur = deepcopy(self.iteridx)
        self.iteridx[0] += 1
        if self.iteridx[0] == self.npts[0]:
            self.iteridx[0] = 0
            self.iteridx[1] += 1
            if self.iteridx[1] == self.npts[1]:
                self.iteridx[0] = 0
                self.iteridx[1] = 0
                self.iteridx[2] += 1
                if self.iteridx[2] == self.npts[2]:
                    self.iteridx = None
        return self.itervalue(cur)

    def reverse(self, k=None):
        'Iterator to return  data points in reverse order'
        a = range(self.npts[0])
        b = range(self.npts[1])
        if k is None:
            c = range(self.npts[2])
        else:
            c = (k,)
        for i in a:
            for j in b:
                for k in c:
                    yield self.itervalue((i, j, k))

    def range(self, k=None):
        'Iterator to return grid indeces in normal order'
        a = range(self.npts[0])
        b = range(self.npts[1])
        if k is None:
            c = range(self.npts[2])
        else:
            c = (k,)
        for k in c:
            for j in b:
                for i in a:
                    yield i, j, k

    def rrange(self, k=None):
        'Iterator to return grid indeces in reverse order'
        a = range(self.npts[0])
        b = range(self.npts[1])
        if k is None:
            c = range(self.npts[2])
        else:
            c = (k,)
        for i in a:
            for j in b:
                for k in c:
                    yield i, j, k

class Grid(GridIterator):
    def __init__(self, l, 
            npts=None,
            step=None,
            basis=None, 
            origin=(0.0, 0.0, 0.0), 
            distribution = 'even'):
        GridIterator.__init__(self)
        self.points = []
        if basis:
            self.basis = basis
        else:
            self.basis = np.identity(3)
        self.l = np.array(l)
        self.origin = np.array(origin)
        self.distribution = 'even'
        self.radius = None
        self.step = np.zeros(3)
        self.bond_ortho = np.zeros(3) # Needed for bond grids
        if step and npts:
            raise ValueError('Both step and number of points specified!')
        if not step and not npts:
            raise ValueError(
                    'Either step or number of points must be specified!')
        if npts:
            if isinstance(npts, int):
                npts = (npts, npts, npts)
            step = self._calc_step(npts)
        else:
            if isinstance(step, float) or isinstance(step, int):
                step = (step, step, step)
            npts = self._calc_npts(step)

        self.npts = npts
        self.step = step
        self._init_basis_vectors(True)
        self._init_grid(self.distribution)

    def _calc_step(self, npts):
        step = np.zeros(3)
        for i in range(3):
            if npts[i] < 1:
                npts[i] = 1
            if npts[i] == 1:
                s = self.l[i]
            else:
                s = float(self.l[i])/float(npts[i]-1)
            step[i] = s
        return step

    def _calc_npts(self, step):
        npts = (0, 0, 0)
        for i in range(3):
            if step[i] == 0 or step[i] < 0.0:
                step[i] = 1.0
            npts[i] = int(self.l[i]/step[i])
            if npts[i] == 0:
                npts[i] == 1
        return npts

    def __str__(self):
        s = 'origin  = {}\n'.format(str(self.origin))
        s += 'i       = {}\n'.format(str(self.basis[:, 0]))
        s += 'j       = {}\n'.format(str(self.basis[:, 1]))
        s += 'k       = {}\n'.format(str(self.basis[:, 2]))
        s += 'lenghts = {}\n'.format(str(self.l))
        s += 'ortho   = {}\n'.format(str(self.bond_ortho))
        return s

    def _init_basis_vectors(self, orthogonal = True):
        self.basis[:, 2] = np.cross(self.basis[:, 0], self.basis[:, 1])
        for i in np.arange(3):
            self.basis[:, i] = self._norm(self.basis[:, i])
        self.bond_ortho = self.basis[:, 2]
        if orthogonal:
            self.orthogonalize_basis()
    
    def _init_grid(self, distr = 'even'):
        self.points = []
        for i in range(3):
            self.points.append(np.zeros(self.npts[i]))
        o = self.origin
        if distr.lower() == 'even':
            for n in range(3):
                self.points.append(np.zeros(self.npts[n]))
                for i in np.arange(self.npts[n]):
                    self.points[n][i] = o[n] + float(i) * self.step[n]
        else:
            raise RuntimeError('Not implemented yet!')

    def orthogonalize_basis(self):
        '''Orthogonalize basis vecotrs. Returns True if the basis was
        orthogonal, and false otherwise.'''
        x = np.dot(self.basis[:, 0], self.basis[:, 1])
        if math.fabs(x) > 0.0:
            v = np.cross(self.basis[:, 0], self.basis[:, 2])
            self.basis[:, 1] = self._norm(v)
            return False
        return True

    def _norm(self, v):
        n = np.linalg.norm(v)
        return v / n

    def get_axis(self, n=None):
        if n:
            return self.points[n]
        return self.points

    def get_basis(self):
        return self.basis

    def get_i(self):
        return self.basis[:, 0]

    def get_j(self):
        return self.basis[:, 1]

    def get_k(self):
        return self.basis[:, 2]

    def set_radius(self, r):
        if r < 0.0:
            self.radius = None
        else:
            self.radius = r

    def rotate(self, angle):
        raise RuntimeError('Not implemented yet')

    def gridpoint(self, r):
        i, j, k = r
        return self.points[0][i] * self.basis[:, 0] + \
            self.points[1][j] * self.basis[:, 1] + \
            self.points[2][k] * self.basis[:, 2] 

    # Definition for the iterator base class
    itervalue = gridpoint

    def get_normal(self):
        return self.get_k()

    def size(self):
        return self.npts

    def get_steps(self):
        return self.step

    def get_origin(self):
        return self.origin

    def lenght(self):
        return self.l

    def get_center(self):
        v1 = self.gridpoint(self.npts[0] - 1, 0, 0)
        v2 = self.gridpoint(0, self.npts[1] - 1, 0)
        return (v1 + v2) * 0.5

    def get_ortho(self):
        return self.bond_ortho

    def is3d(self):
        if self.npts[2] > 1:
            return True
        return False
    

class BondGrid(Grid):
    def __init__(self, 
            bond, 
            height, 
            width, 
            npts,
            fixpoint=None,
            distance=None, 
            radius=None):
        self.l = np.zeros(3)
        self.basis = np.ndarray((3,3))
        self.distribution = 'even'
        self.radius = radius
        if len(bond) != 2:
            raise ValueError('Number of atoms in bond != 2')
        if fixpoint is None:
            fixpoint = np.zeros(3)
        self.l[0] = sum(height)
        self.l[1] = sum(width)
        self.l[2] = 0.0

        if distance is None:
            distance = bond[0].bond_distance(bond[1])*0.5

        # figure out the "orthogonal" axis wrt. the magnetic field
        v1 = bond[0].get_coord() - fixpoint
        v2 = bond[1].get_coord() - fixpoint
        self.bond_ortho = np.cross(v1, v2)
        if np.linalg.norm(self.bond_ortho) < 10.0e-6:
            raise RuntimeError('Basis vectors are linearly dependent!')
        self.bond_ortho = self._norm(self.bond_ortho)

        v3 = self._norm(v2-v1)
        v1 = -1.0 * self.bond_ortho 
        v2 = self._norm(np.cross(v3, v1))
        oo = bond[0].get_coord() + distance * v3
        self.origin = oo - width[1] * v2 - height[1] * v1
        self.center = oo

        self.basis[:, 0] = v1
        self.basis[:, 1] = v2
        self.basis[:, 2] = v3

        step = 0.5
        if npts:
            if isinstance(npts, int):
                npts = (npts, npts, npts)
            step = self._calc_step(npts)
        else:
            if isinstance(step, float) or isinstance(step, int):
                step = (step, step, step)
            npts = self._calc_npts(step)

        self.npts = npts
        self.step = step
        self._init_basis_vectors(True)
        self._init_grid(self.distribution)

    def __str__(self):
        s = Grid.__str__(self)
        s += 'center  = {}'.format(str(self.center))
        return s

if __name__ == '__main__':
    g = Grid(l=(4, 5, 6), npts=(2,3,2))
    print g
    a1 = Atom((-1.0, 0.0, 0.0))
    a2 = Atom((1.0, 0.0, 0.0))
    bg = BondGrid((a1, a2), height=(5,5), width=(2,2), npts=(2,3,2),
            fixpoint=(0,0,1))
    print bg

#    for i, j, k in g.range(0):
#        print i, j, k
#    print
#    for r in g.reverse():
#        print r
#    print
#    print g.gridpoint((0, 0, 0))

# vim:et:ts=4:sw=4
