#
# Jonas Juselius <jonas.juselius@uit.no> 2004, 2012
#
# TODO: Read atom data from MOL files
#       Add classes for XYZMolecule, MOLMolecule, etc...
#
from atom import Atom
from elements import Element

class Molecule:
    def __init__(self, atoms):
        if isinstance(atoms, str):
            self.read_xyz(atoms)
        else:
            self.atoms = []
            for i in atoms:
                self.atoms.append(atoms[i])

    def __getitem__(self, i):
        return self.atoms[i]

    def __str__(self):
        s = ''
        for i in self.atoms:
            s += str(i) + '\n'
        return s

    def append(self, a):
        self.atoms.append(a)

    def get_atom(self, i):
        return self.atoms[i]

    def read_xyz(self, file):
        "Read XYZ file, return list of Atoms"
        atoms=[]
        with open(file, 'r') as f:
            natoms=int(f.readline())  # First line is numer of atoms
            f.readline()  # Empty/title line
            for i in f:
                data = i.split()
                sym = data[0]
                coord = map(float, data[1:])
                atoms.append(Atom(coord, sym))
        if len(atoms) != natoms:
            raise RuntimeError('Atom number mismatch in XYZ file.')
        self.atoms = atoms
        self.change_units('a2au')

    def write_xyz(self, file):
        "Write a XYZ file of atoms indexes"
        self.change_units('au2a')
        with open(file, 'w') as f:
            print >> f, len(self.atoms)
            print >> f
            for i in self.atoms:
                print >> f, i
        self.change_units('a2au')

    def change_units(self, conv):
        for i in self.atoms:
            i.change_units(conv)

if __name__ == '__main__':
    mol = Molecule('coord.xyz')
    print mol[0]
