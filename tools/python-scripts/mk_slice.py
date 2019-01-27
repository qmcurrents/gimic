#!/usr/bin/python 
import sys, os
import argparse
from itertools import islice
from subprocess import call
import numpy as np
import math
import time

# usage: 
# python file.py --atoms at1,at2 --down val --up val --ins val --out val
# needs coord.xyz, XDENS and MOL in the same directory
#
start_time = time.time()

# read input from command line
parser = argparse.ArgumentParser()
parser.add_argument('--atoms', type=str)
parser.add_argument('--fix', type=str)
parser.add_argument('--down', type=float)
parser.add_argument('--up', type=float)
parser.add_argument('--ins', type=float)
parser.add_argument('--out', type=float)
args = parser.parse_args() 
var = args.atoms.split(',')
idx = [int(x) for x in var]
if args.fix != "F":
    fixatom = args.fix
    idxfix = int(args.fix)
    print "fixatom", fixatom
#
height = [args.down,args.up]
width = [args.ins,args.out]
print height
print width

# assume height is kept fixed and vary only width
step = 0.1
length = np.abs(width[0]) + width[1]
print "length =", length
nstep = length/0.1
print "number of steps", nstep
imax = int(nstep)
print imax
start=width[0]
end=width[1]

for i in range(0,imax):
    workdir = "sc"+"_%04d"%(i)
    os.mkdir(workdir) 
    os.chdir(workdir)
    os.system(" cp ../coord.xyz . ")
    os.system(" cp ../get_input.py . ")
    os.system(" ln -s ../XDENS ")
    os.system(" ln -s ../MOL ")
    l0 =" python get_input.py"
    l1 =" --atoms "+str(idx[0])+","+str(idx[1])
    l1f =" --fix "+str(idxfix)
    l2 =" --down "+str(height[0])+" --up "+str(height[1])
    l3 =" --ins "+str(start)+" --out "+str(start+step)
    line = l0+l1+l1f+l2+l3
    print line
    start = start + step
    os.system(line) 
    runcalc = "gimic > gimic_h"+"_%04d"%(i)
    os.system(runcalc)
    os.system("cp gimic_h* ../")
    os.system("rm coord.xyz")
    os.system("rm show_fixpoint.xyz")
    os.chdir("../")
    print("--- %s seconds ---" % (time.time() - start_time))

os.system("rm -r sc* ")


