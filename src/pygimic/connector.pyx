# Generic GIMIC interface
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

from gimic_exceptions import NotAvailable
from connector cimport GimicConnector

cdef class GimicConnector:
    cpdef jvector(self, r):
        raise NotAvailable('jvector()')

    cpdef jtensor(self, r):
        raise NotAvailable('jtensor()')

    cpdef divj(self, r):
        raise NotAvailable('divj()')

    cpdef edens(self, r):
        raise NotAvailable('edens()')

    cpdef set_property(self, prop, val):
        raise NotAvailable('set_property()')

if __name__ == '__main__':
    g = GimicConnector()
    g.jvector((0, 0, 0))

# vim:et:ts=4:
