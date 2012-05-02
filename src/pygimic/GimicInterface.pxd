# Generic GIMIC interface
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

cdef class GimicInterface:
    def jvector(self, r)
    def jtensor(self, r)
    def divj(self, r)
    def edens(self, r)
    def set_property(self, prop, val)

# vim:et:ts=4:
