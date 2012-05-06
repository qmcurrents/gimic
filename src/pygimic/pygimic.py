#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
from grid import Grid
from field import VectorField, ScalarField
from currents import CurrentField
import gimic
from liblondon import london

class GimicDriver:
    def __init__(self, args, inkeys):
        self.args = args
        self.kw = inkeys
        if self.kw['backend'][0].arg[0] == 'london':
            self.gimic = london.London(self.kw['xdens'][0].arg[0])
        else:
            self.gimic = gimic.Gimic(self.kw['basis'][0].arg[0],
                self.kw['xdens'][0].arg[0])

    def run(self):
        print "This is PyGIMIC."

# vim:et:ts=4:sw=4
