#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import os

GIMIC_EXECUTABLE = 'gimic.bin'

class GimicDriver:
    def __init__(self, args, inkeys):
        self.args = args
        self.kw = inkeys

    def run(self):
        print "This is F-GIMIC."

        infile='GIMIC.in' + str(os.getpid())
        fd=open(infile,'w')
        print >>fd, self.kw

        fd.close()
        os.system(GIMIC_EXECUTABLE + ' < ' + infile)
        os.unlink(infile)

# vim:et:ts=4:sw=4
