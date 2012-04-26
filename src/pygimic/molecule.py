from atom import Atom
from elements import Element

class Molecule:
    def __init__(self, atoms):
        if isinstance(atoms, str):
            self.atoms = self.read_xyz(atoms)
        else:
            self.atoms = atoms

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
        return atoms

if __name__ == '__main__':
    mol = Molecule('coord.xyz')
