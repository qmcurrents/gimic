#
# Jonas Juselius <jonas.juselius@uit.no> 2012
#
import numpy as np
from gimic_exceptions import NotImplemented
from grid import Grid, BondGrid
from magnet import Magnet
from atom import Atom
from molecule import Molecule
from field import VectorField, ScalarField
from currents import CurrentField
from plot import MatPlot
from quadrature import FieldQuadrature
from liblondon import london
import gimic

class GimicDriver:
    def __init__(self, args, inkeys):
        print "This is PyGIMIC."
        self.args = args
        self.kw = inkeys
        self.debug = self.kw.getkw('debug')
        self.openshell = False
        self.grid = None
        self.magnet = None
        title = self.kw.getkw('title')
        if len(title) > 0:
            print 'TITLE:', title[0]

        xdens = self.kw.getkw('xdens')[0]
        if self.kw.getkw('backend') == 'london':
            self.gimic = london.London(xdens)
        else:
            mol = self.kw.getkw('basis')[0]
            self.gimic = gimic.Gimic(mol, xdens)
            self.init_gimic_options()

        grid_type = self.kw.getsect('Grid').arg.arg[0] # ugly
        if grid_type == 'bond':
            self.init_bondgrid()
        elif grid_type == 'file':
            raise NotImplemented
        else:
            self.init_grid()
        print "Grid info:"
        print "=========="
        print self.grid
        
        self.init_magnet()
        b = tuple(self.magnet.get_magnet()) 
        self.gimic.set_property('magnet', b)

    def run(self):
        if self.args.dryrun:
            return
        calc = self.kw.getkw('calc')[0]

        if calc == 'cdens':
            j = CurrentField(self.grid, self.gimic)
            p = MatPlot(j)
            p.vector_plot()
            p.stream_plot()
        elif calc == 'rho':
            f = ScalarField(self.grid)
            f.calc(gimic.rho)
            p = MatPlot(f)
            p.scalar_plot()
            p.scalar_plot3d()
        elif calc == 'integral':
            j = CurrentField(self.grid, self.gimic)
            q = FieldQuadrature()
            jp = j.get_normal_flow('+')
            jn = j.get_normal_flow('-')
            Ip = q.integrate(jp)
            In = q.integrate(jn)
            c = self.au2nAt()
            print 'Current in a.u.: {0} (+{1} / {2})'.format(Ip + In, Ip, In)
            print 'Current in nA/T: {0} (+{1} / {2})'.format(
                    (Ip + In) * c, Ip * c, In * c)
        else:
            raise NotImplemented

    def init_bondgrid(self):
        sect = self.kw.getsect('Grid')
        distr = sect.getkw('type')[0]

        bond = map(int, sect.getkw('bond'))
        distance = float(sect.getkw('distance')[0])
        npts = map(int, sect.getkw('grid_points'))
        height = map(float, sect.getkw('height'))
        width = map(float, sect.getkw('width'))
        fixpoint = map(float, sect.getkw('fixpoint'))
        radius = map(float, sect.getkw('radius'))
        
        mol = Molecule('mol.xyz') # dirty
        atom1 = mol[bond[0]]
        atom2 = mol[bond[1]]
        
        self.grid = BondGrid(bond=(atom1, atom2), fixpoint=fixpoint, 
            npts=npts, distance=distance,
            height=height, width=width, 
            distribution=distr, radius=None)

    def init_grid(self):
        sect = self.kw.getsect('Grid')
        distr = sect.getkw('type')[0]
        origin = sect.getkw('origin')
        ivec = sect.getkw('ivec')
        jvec = sect.getkw('jvec')
        l = sect.getkw('lengths')
        npts = sect.getkw('grid_points')

        npts = map(int, npts)
        l = np.array(map(float, l))
        origin = np.array(map(float, origin))
        ivec = np.array(map(float, ivec))
        jvec = np.array(map(float, jvec))
        kvec = np.cross(ivec, jvec)
        basis = np.ndarray((3,3))
        basis[:, 0] = ivec
        basis[:, 1] = jvec
        basis[:, 2] = kvec

        self.grid = Grid(l=l, origin=origin, npts=npts, basis=basis, 
                distribution=distr)

    def init_magnet(self):
        if self.kw.findkw('magnet_axis').is_set():
            mag = self.kw.getkw('magnet_axis')[0]
        elif self.kw.findkw('magnet').is_set():
            mag = self.kw.getkw('magnet')[0]
        else:
            raise ValueError('No magnetic field specified')
        self.magnet = Magnet(mag, self.grid)

# Most options are still unimplemented. Easy to do, requires a bit time.
    def init_gimic_options(self):
        adv = self.kw.getsect('Advanced')

        if self.kw.getkw('openshell')[0] == 'True':
            self.openshell = True
        if self.openshell:
            self.gimic.set_property('uhf', 1)

        screening = adv.getkw('screening')[0]
        if screening:
            thrs = float(adv.getkw('screening_thrs')[0])
            self.gimic.set_property('screening', thrs)

## Future stuff
#        spherical = adv.getkw('spherical')[0]
#        giao = adv.getkw('GIAO')[0]
#        diamag = adv.getkw('diamag')[0]
#        paramag = adv.getkw('paramag')[0]

#        modens = self.kw.getkw('density')
#        mofile = self.kw.getkw('mofile')
#        molist = self.kw.getkw('mos')

    def init_london_options(self):
        pass

    def au2nAT(self, au=1.0):
        'Derivation of a.u. to nA/t'
        aulength=0.52917726e-10
        auspeedoflight=137.03599e0
        speedoflight=299792458.e0
        aucharge=1.60217733e-19
        hbar=1.05457267e-34

        autime=aulength*auspeedoflight/speedoflight
        autesla=hbar/aucharge/aulength/aulength
        audjdb=aucharge/autime/autesla

        return au*audjdb*1.e+09 # nA/T


# vim:et:ts=4:sw=4
