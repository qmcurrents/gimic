# Generic GIMIC interface
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

from gexceptions import NotAvailable

class GimicInterface:
    def jvector(self, r):
        raise NotAvailable("jvector")

    def jtensor(self, r):
        raise NotAvailable("jtensor")

    def divj(self, r):
        raise NotAvailable("divj")

    def edens(self, r):
        raise NotAvailable("edens")

    def set_property(self, prop, val):
        raise NotAvailable("set_property")

# vim:et:ts=4:
