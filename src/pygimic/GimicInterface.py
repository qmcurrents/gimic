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
    def calc_jvector(r, jv):
        raise NotAvailable

    def calc_jtensor(r, jt):
        raise NotAvailable

    def calc_divj(r, dj):
        raise NotAvailable

    def calc_edens(r, ed):
        raise NotAvailable

    def set_property(prop, val):
        raise NotAvailable

# vim:et:ts=4:
