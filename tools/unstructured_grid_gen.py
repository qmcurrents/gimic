#!/usr/bin/python3
#!@PYTHON_EXECUTABLE@
#

import sys
import string
import re
import math
sys.path.append('@PYTHON_INSTDIR@')
from QCTools import QCTools
from QCTools.Elements import PeriodicTable as PTab
from argparse import ArgumentParser

def parseCommandline():
    parser = ArgumentParser(prog="unstructured_grid_gen.py", description="first create a body centred cubic lattice, then prune using atomic coordinates", add_help=True)
    parser.add_argument("--coord-file",   help="",      nargs=1, dest="coord_file", action="store", type=str, required=False, default=['coord'])
    parser.add_argument("--cutoff",       help="[a_0]", nargs=1, dest="cutoff",   action="store", type=float,  required=False, default=5.0)
    # parser.add_argument("--cutoff-perpB", help="[a_0]", nargs=1, dest="cutoff_pB",  action="store", type=float,  required=False, default=8.0 )
    # parser.add_argument("--B",            help="",      nargs=3, dest="B",          action="store", type=float,  required=False, default=[0, 0, 1] )
    parser.add_argument("--grid-spacing", help="[a_0]", nargs=1, dest="spacing",    action="store", type=float,  required=False, default=0.5 )
    argparse = parser.parse_args()
    argparse.spacing /= math.sqrt(3)/2.0
    print(argparse.coord_file)
    return argparse

def logical_xor(str1, str2):
    return bool(str1) ^ bool(str2)

def main():
    args = parseCommandline()
    print(args)  

    atoms = QCTools.readCoord(args.coord_file[0])
    atoms = atoms[1:]
    print(atoms)

    # find min/max x/y/z
    min_x = sys.maxsize
    min_y = sys.maxsize
    min_z = sys.maxsize
    max_x = -sys.maxsize
    max_y = -sys.maxsize
    max_z = -sys.maxsize
    print(min_x, max_x, min_y, max_y, min_z, max_z)

    for atom in atoms:
        print(atom)
        min_x = min(atom.coord[0], min_x)
        min_y = min(atom.coord[1], min_y)
        min_z = min(atom.coord[2], min_z)
        max_x = max(atom.coord[0], max_x)
        max_y = max(atom.coord[1], max_y)
        max_z = max(atom.coord[2], max_z)
    print(min_x, max_x, min_y, max_y, min_z, max_z)
    
    # number of layers in each direction
    x_layers = 2 * (max_x - min_x + 2*args.cutoff) / float(args.spacing)
    y_layers = 2 * (max_y - min_y + 2*args.cutoff) / float(args.spacing)
    z_layers = 2 * (max_z - min_z + 2*args.cutoff) / float(args.spacing)
    print(x_layers, y_layers, z_layers)
    
    grid = [] 
    for z in range(0, math.ceil(z_layers)):
      z_coord = z * float(args.spacing) / 2.0 + (min_z - args.cutoff)
      z_even = (z%2 == 0)
      for y in range(0, math.ceil(y_layers)):
        y_coord = y * float(args.spacing) + (0.0 if z_even else float(args.spacing) / 2.0) + (min_y - args.cutoff)
        y_even = (y%2 == 0)
        for x in range(0, math.ceil(x_layers)):
          x_coord = x * float(args.spacing) + (0.0 if z_even else float(args.spacing) / 2.0) + (min_x - args.cutoff)
          for atom in atoms:
            dist = math.sqrt( math.pow((atom.coord[0] - x_coord),2) +
                              math.pow((atom.coord[1] - y_coord),2) +
                              math.pow((atom.coord[2] - z_coord),2) )
            # print (dist, args.cutoff)
            if( dist < args.cutoff):
              grid.append([x_coord, y_coord, z_coord])
              print('X', x_coord, y_coord, z_coord)
              break
          
    write_xyz(grid)

    # for i in range(1,len(atoms)):
    #     a=atoms[i]
    #     sym=str.upper(a.element.symbol)
    #     if sym == 'Q':
    #         sys.stderr.write("\n*** ERROR ***\nAtoms without basis functions should be removed from the control file \n\n")
    #         sys.exit(1)
    #     print(("%2s%2i%20.12f%20.12f%20.12f" % \
    #         (sym, 1, a.coord[0], a.coord[1], a.coord[2])))


if __name__ == '__main__':
    main()

