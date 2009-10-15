#!/usr/bin/env python

import sys, os, re

nsteps = 20
delta = 0.25
start = 1.0
end = -start + delta

for i in range(nsteps):
	xstart = start - delta*i
	xend = end + delta*i
	sed="sed 's/@start@/%f/; s/@end@/%f/' gimic.int >gimic.%d.inp" % \
	(xstart, xend, i)
	os.system(sed)
	

