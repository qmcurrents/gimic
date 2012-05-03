# Generic GIMIC interface
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

cdef class GimicConnector:
    cpdef jvector(self, r)
    cpdef jtensor(self, r)
    cpdef divj(self, r)
    cpdef edens(self, r)
    cpdef set_property(self, prop, val)

# vim:et:ts=4:
