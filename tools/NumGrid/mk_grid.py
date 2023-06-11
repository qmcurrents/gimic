#!/usr/bin/env python

# function for avoiding CFOUR trouble 
def add_blank(filename, old, new):
    with open(filename) as f:
        s = f.read()

    with open(filename, 'w') as f:
        s = s.replace(old, new)
        f.write(s)

# replace "#" with " #" because of CFOUR problem
ll = add_blank("MOL", "#", " #")


# generate grid for molecule using numgrid external library

import numgrid
import numpy as np
import math
from itertools import islice

# some fix variables 
radial_precision = 1.0e-12
min_num_angular_points = 86
max_num_angular_points = 302
hardness = 3

# read in MOL file and extract relevant information
# number of centres
# number of proton charges
# coordinates in bohr
# from the basis set for each atom we need:
# maximum exponent alpha_max only s 
# maximum l quantum number
# minimum exponent alpha_min, min s, min p, min d.... 

center_coordinates_bohr = []
proton_charges = []
max_l_quantum_numbers = []
block = []
element = []
alpha = []
alpha_max = []
alpha_min = []
alpha_min_tmp = []
label = []
coordinates = []
weights = []

idxb = 0

fin = "MOL"
# replace "#" with " #" because of CFOUR problem
# add_blank(fin, "#", " #")

with open(fin) as f:
    # skip first three lines
    for _ in range(3): 
        next(f)
    l = f.readline().strip().split() 
    num_centers = int(l[0])
    # skip one line
    for _ in range(1):
        next(f)
    for n in range(num_centers):
        l = f.readline().strip().split() 
        tmp = float(l[0])
        proton_charges.append(int(tmp))
        lmax_plus_one = int(l[2])
        max_l_quantum_numbers.append(lmax_plus_one -1)
        for i in range(lmax_plus_one):
            block.append(int(l[3+i]))
            idxb = idxb + 1

        l = f.readline().strip().split() 
        element.append(l[0])
#       must be a list of tuples ()
        xyz  = (float(l[2]), float(l[3]), float(l[4]))
        center_coordinates_bohr.append(xyz)
    
        idxk = 0
        for j in range(idxb):
            for i in range(block[j]):
                l = f.readline().strip().split()
                nfunc = int(l[0]) 
                for k in range(nfunc):
                    l = f.readline().strip().split()
                    alpha.append(float(l[0]))
                    idxk = idxk + 1
            alpha_min_tmp.append(alpha[idxk-1])
            label.append(j)
#       build requested dictionary for alpha_min
        tmp = dict(zip(label,alpha_min_tmp)) 
        alpha_min.append(tmp)
        alpha_max.append(alpha[0])
        # clean up for new element
        idxb = 0
        block[:] = []
        alpha[:] = []
        alpha_min_tmp[:] = []

# save coordinates in au
with open('coord.au', 'w') as f1:
    for n in range(num_centers):
       # map tuple to string for printing
       string = " ".join(map(str, center_coordinates_bohr[n]))
       f1.write(string+"\n")

# now calculate the grid after all input has been extracted from MOL

fout_grid = "gridfile.grd"
fout_weights = "grid_w.grd" 
fout_nelpts = "nelpts.info" 
f1 = open(fout_grid, "w")
f2 = open(fout_weights, "w")
f3 = open(fout_nelpts, "w")

for n in range(num_centers):
    print("n", n)
    print("proton_charges[n]", proton_charges[n])
    print("alpha_max[n]", alpha_max[n])
    print("max_l_quantum_numbers[n]", max_l_quantum_numbers[n])
    print("alpha_min[n]", alpha_min[n])
#   get atom grid using explict basis set parameters
    coordinates, weights = numgrid.atom_grid(
        alpha_min[n],
        alpha_max[n],
        radial_precision,
        min_num_angular_points,
        max_num_angular_points,
        proton_charges,
        n,
        center_coordinates_bohr,
        hardness,
    )

    num_points = len(coordinates)
    print("center", n)
    print("num_points", num_points)
    # collect center index and related number of points
    f3.write(str(n+1) + " " + str(num_points) + "\n")
    # print grid coordinates and weights on file
    for k in range(num_points): 
        # map tuple to string for printing
        string = " ".join(map(str, coordinates[k]))
        f1.write(string+"\n")
        f2.write(str(weights[k]) + "\n")

f1.close()
f2.close()
f3.close()
