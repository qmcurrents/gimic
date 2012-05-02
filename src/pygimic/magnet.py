import numpy as np
from grid import Grid

class Magnet:
    def __init__(self, b, grid = None, ortho = False):
        self.grid = grid
        self.isortho = ortho
        if isinstance(b, str):
            if not self.grid:
                raise ValueError('Grid object not defined.')
            self.mag = self._get_direction(b)
        else:
            self.mag = np.array(b)
        self._check_field_direction()

    def __str__(self):
        return " ".join([str(i) for i in self.mag])

    def get_magnet(self):
        return self.mag

    def set_orthogonal(self, p = True):
        self.isortho = p
        self._check_field_direction()
    
    def is_orthogonal():
        return self.isortho

    def _get_direction(self, b):
        b = b.strip()
        if b[0] == '-':
            b = b[1:]
            sign = -1.0
        elif b[0] == '+':
            b = b[1:]
            sign = 1.0
        else:
            sign = 1.0

        self.mag = np.zeros(3)
        if b == 'T':
            b = 'ortho'
        b = b.lower()
        if b == 'i':
            self.mag = np.array(self.grid.get_i())
        elif b == 'j':
            self.mag = np.array(self.grid.get_j())
        elif b == 'k':
            self.mag = np.array(self.grid.get_k())
        elif b == 'x':
            self.mag[0] = 1.0
        elif b == 'y':
            self.mag[1] = 1.0
        elif b == 'z':
            self.mag[2] = 1.0
        elif b == 'ortho':
            self.mag = self.grid.get_ortho()
        else:
            raise ValueError('Invalid axis specification')
        self.mag *= sign
        return self.mag

    def _check_field_direction(self):
        if self.isortho or not self.grid:
            return
        kvec = self.grid.get_k()
        x = np.dot(kvec, self.mag)
        if x > 0.0:
            self.mag *= -1.0
            print 'Left handed coordinate system reversing magnetic field' 
        if abs(x) > 1.0e-12 and abs(x)-1.0 > 1.0e-12:
            print 'Warning! Magnetic field not orthogonal to the grid!' 


if __name__ == '__main__':
    g = Grid()
    b = Magnet([1.1, 2.2, 3.3])
    print b
    b = Magnet('i', g)

# vim:et:sw=4:ts=4
