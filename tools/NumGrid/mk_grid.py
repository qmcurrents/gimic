#!/usr/bin/env python

# generate grid for molecule using numgrid external library

import numgrid
import numpy as np
import math
from itertools import islice

# some fix variables 
radial_precision = 1.0e-12
min_num_angular_points = 86
max_num_angular_points = 302

# read in MOL file and extract relevant information
# number of centres
# number of proton charges
# coordinates in bohr
# from the basis set for each atom we need:
# maximum exponent alpha_max only s 
# maximum l quantum number
# minimum exponent alpha_min, min s, min p, min d.... 

proton_charges = []
max_l_quantum_numbers = []
block = []
element = []
x_coordinates_bohr = []
y_coordinates_bohr = []
z_coordinates_bohr = []
alpha = []
alpha_max = []
alpha_min = []
alpha_min_tmp = []
idxb = 0

fin = "MOL"
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
        x_coordinates_bohr.append(float(l[2]))
        y_coordinates_bohr.append(float(l[3]))
        z_coordinates_bohr.append(float(l[4]))
    
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

        alpha_max.append(alpha[0])
        alpha_min.append(alpha_min_tmp[:])
        # clean up for new element
        idxb = 0
        block[:] = []
        alpha[:] = []
        alpha_min_tmp[:] = []

# save coordinates in au
with open('coord.au', 'w') as f1:
    for n in range(num_centers):
        f1.write(str(x_coordinates_bohr[n])+" "+str(y_coordinates_bohr[n])+" "+str(z_coordinates_bohr[n])+"\n")

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

    context = numgrid.new_atom_grid(radial_precision, 
            min_num_angular_points, max_num_angular_points, 
            proton_charges[n], alpha_max[n], max_l_quantum_numbers[n], 
            alpha_min[n] )

    num_points = numgrid.get_num_grid_points(context)
    print("number of points", num_points)
    # collect center index and related number of points
    f3.write(str(n+1) + " " + str(num_points) + "\n")
    # generate an atomic grid in the molecular environment
    x, y, z, w = numgrid.get_grid(context, num_centers, n, 
            x_coordinates_bohr, y_coordinates_bohr, 
            z_coordinates_bohr, proton_charges)

    for k in range(num_points):
        f1.write(str(x[k]) + " " + str(y[k]) + " " + str(z[k]) + "\n")
        f2.write(str(w[k]) + "\n")

    numgrid.free_atom_grid(context)

f1.close()
f2.close()
f3.close()



