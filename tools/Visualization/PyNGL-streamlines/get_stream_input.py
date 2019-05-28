#! /usr/bin/env python

# purpose: generate input file for streamline 2D script of Radovan Bast
#          
# usage: best in combination with plot script
# input: jvec.txt and coord.xyz (Angstroem) 
# output: jvec_2d.txt and coord.xy
#

import numpy as np
import sys

filename = "jvec.txt"
fout = "jvec_2d.txt"

infile   = open(filename, 'r')

x = []
y = []
z = []

jx = []
jy = []
jz = []

for line in infile:
    # skip empty lines in input file
    if not line.strip():
        continue
    else:
        a = float(line.split()[0]) 
        b = float(line.split()[1]) 
        c = float(line.split()[2]) 
        u = float(line.split()[3]) 
        v = float(line.split()[4]) 
        w = float(line.split()[5])
    
        x.append(a)
        y.append(b)
        z.append(c)
        jx.append(u)
        jy.append(v)
        jz.append(w)
    
        icount = len(jz)
# close file
infile.close()
xx = np.array(x)
yy = np.array(y)
jxx = np.array(jx)
jyy = np.array(jy)

npts = icount
# print npts

f1 = open(fout, 'w')
for i in range(npts):
    print >> f1, ("%14.6e  %14.6e  %14.6e  %14.6e" % ( yy[i],xx[i],jyy[i],jxx[i]))
f1.close()

filename = "coord.xyz"
fout = "coord.xy"

label = []
x = []
y = []
z = []

# now the coord file
with open(filename, 'r') as infile:
    #skip first line
    next(infile)
    for line in infile:
        # skip empty lines in input file
        if not line.strip():
            continue
        else:
            a = line.split()[0] 
            b = float(line.split()[1]) 
            c = float(line.split()[2]) 
            u = float(line.split()[3]) 
        
            label.append(a)
            x.append(b)
            y.append(c)
            z.append(u)
        
            icount = len(z)
# close file
infile.close()

xx = np.array(x)
yy = np.array(y)
ll = np.array(label)

npts = icount
# print npts

f1 = open(fout, 'w')
for i in range(npts):
    print >> f1, ("%s  %14.6f  %14.6f " % ( ll[i], yy[i], xx[i] ))

f1.close()


