#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
from grid import Grid
from field import VectorField, ScalarField
from currents import CurrentField

class GimicDriver:
    def __init__(self, args, inkeys):
        self.args = args
        self.kw = inkeys

    def run(self):
        print "This is PyGIMIC."

# vim:et:ts=4:sw=4
