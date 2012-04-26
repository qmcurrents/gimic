#
# Basic grid class
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
import math

class Grid:
    def __init__(self, l, 
            npts = None,
            step = None,
            basis=None, 
            origin=(0.0, 0.0, 0.0), 
            distribution = 'even'):
        self.points = []
        if basis:
            self.basis = basis
        else:
            self.basis = np.identity(3)
        self.l = np.array(l)
        self.origin = np.array(origin)
        self.ortho = np.zeros(3)
        self.distribution = 'even'
        self.radius = None
        self.step = np.zeros(3)
        self.npts = np.ones(3, dtype=np.int32)
        if step and npts:
            raise ValueError('Both step and number of points specified!')
        if not step and not npts:
            raise ValueError(
                    'Either step or number of points must be specified!')
        if npts:
            step = []
            if isinstance(npts, int):
                npts = (npts, npts, npts)
            for i in np.arange(3):
                if npts[i] == 0:
                    self.npts[i] = 1
                else:
                    self.npts[i] = npts[i]
                if self.npts[i] == 1:
                    s = self.l[i]
                else:
                    s = float(self.l[i])/float(npts[i]-1)
                step.append(s)

        if isinstance(step, float) or isinstance(step, int):
            step = (step, step, step)
        for i in np.arange(len(self.step)):
            self.step[i] = step[i]

        print self.step
        self._init_basis_vectors(True)
        self._init_grid(self.distribution)

    def __str__(self):
        return str(self.basis)

    def _init_basis_vectors(self, orthogonal = True):
        self.basis[:, 2] = np.cross(self.basis[:, 0], self.basis[:, 1])
        for i in np.arange(3):
            self.basis[:, i] = self._norm(self.basis[:, i])
        self.ortho = self.basis[:, 2]
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
        n = math.sqrt(np.dot(v, v.conj()))
        return v / n

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

    def point(self, i, j, k = 0):
        return self.points[0][i] * self.basis[:, 0] + \
            self.points[1][j] * self.basis[:, 1] + \
            self.points[2][k] * self.basis[:, 2] 

    def get_normal(self):
        return self.get_k()

    def size(self):
        return self.npts

    def lenght(self):
        return self.l

    def get_center(self):
        v1 = self.point(self.npts[0] - 1, 0, 0)
        v2 = self.point(0, self.npts[1] - 1, 0)
        return (v1 + v2) * 0.5

    def get_ortho(self):
        return self.ortho

    def is3d(self):
        if self.npts[2] > 1:
            return True
        return False

if __name__ == '__main__':
    g = Grid(l=(4, 5, 6), npts=(2,2,1))
    print g.point(0, 0, 0)
    print g.point(1, 1, 0)

# vim:et:ts=4:sw=4
