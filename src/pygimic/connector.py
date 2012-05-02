# Generic GIMIC interface
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

from gimic_exceptions import NotAvailable
from connector cimport GimicConnector

cdef class GimicConnector:
    def jvector(self, r):
        raise NotAvailable('jvector()')

    def jtensor(self, r):
        raise NotAvailable('jtensor()')

    def divj(self, r):
        raise NotAvailable('divj()')

    def edens(self, r):
        raise NotAvailable('edens()')

    def set_property(self, prop, val):
        raise NotAvailable('set_property()')

if __name__ == '__main__':
    g = GimicConnector()
    g.jvector((0, 0, 0))

# vim:et:ts=4:
