# Generic GIMIC interface
#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#

class NotAvailable(Exception):
    gimic_backend = 'GIMIC'
    def __init__(self, value):
        self.value = value
     def __str__(self):
         return "Requested method is not availbale for chosen backend:"
            "%s:\n  %s" % (gimic_backend, repr(self.value))

class GimicInterface:
    def jvector(self, r):
        raise NotAvailable

    def jtensor(self, r):
        raise NotAvailable

    def divj(self, r):
        raise NotAvailable

    def edens(self, r):
        raise NotAvailable

    def set_property(self, prop, val):
        raise NotAvailable

# vim:et:ts=4:
