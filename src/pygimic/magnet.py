import numpy as np
from grid import Grid

class Magnet:
    def __init__(self, b, grid = None, ortho = False):
        self.grid = grid
        self.ortho = ortho
        if isinstance(b, str):
            if not self.grid:
                raise ValueError('Grid object not defined.')
            self.mag = self._get_direction(b.lower())
        else:
            self.mag = np.array(b)
        self._check_field_direction()

    def __str__(self):
        return " ".join([str(i) for i in self.mag])

    def get_magnet(self):
        return self.mag

    def set_orthogonal(self, p = True):
        self.ortho = p
        self._check_field_direction()
    
    def is_orthogonal():
        return self.ortho

    def _get_direction(self, b):
        b = b.strip()
        if len(b) == 1:
            sign = 1.0
        elif len(b) == 2:
            if b[0] == '-':
                sign = -1.0
            elif b[0] == '+':
                sign = 1.0
            else:
                raise ValueError('Invalid axis specification')
            b = b[1]
        else:
            raise ValueError('Invalid axis specification')
            
        self.mag = np.zeros(3)
        if b == 'i':
            mag = self.grid.get_i()
        elif b == 'j':
            mag = self.grid.get_j()
        elif b == 'k':
            mag = self.grid.get_k()
        elif b == 'x':
            self.mag[0] = 1.0
        elif b == 'y':
            self.mag[1] = 1.0
        elif b == 'z':
            self.mag[2] = 1.0
        elif b == 't':
            self.mag = np.ones(3)
        else:
            raise ValueError('Invalid axis specification')
        self.mag *= sign
        return self.mag

    def _check_field_direction(self):
        if self.ortho or not self.grid:
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
