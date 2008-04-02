#!/usr/bin/env python

import string, math, sys, re
from getopt import getopt

vmin=0.0
vmax=1.e10
scale=1.e0
fd=None
ofd=sys.stdout

def setup():
	global scale, vmin, vmax, fd, ofd
	jvecf="JVEC.txt"
	outfile=None
	try:
		(opts,args)=getopt(sys.argv[1:], "o:s:m:M:h")
	except Exception, inst:
		print "Error: ", inst
		sys.exit(1)

	for i in opts:
		if i[0] == '-o':
			outfile=i[1]
			continue 
		if i[0] == '-s':
			scale=float(i[1])
			continue 
		if i[0] == '-M':
			vmin=float(i[1])
			continue 
		if i[0] == '-m':
			vmax=float(i[1])
			continue 
		if i[0] == '-h':
			print "usage: %s [-s scale] [-M min] [-m max] [file]" % (
					sys.argv[0])
			sys.exit(1)

	if len(args) != 0:
		jvecf=args[0]

	fd=open(jvecf,'r')
	if outfile is not None:
		ofd=open(outfile,'w')

def main():
	global scale, vmin, vmax, fd, ofd
	setup()
	amin=0.e10
	amax=0.e0

	tmp=fd.readlines()
	lines=[] 
	for i in tmp:
		if re.search(r'^\s*#', i):
			continue
		lines.append(i)

	for i in lines:
		(x,y,z,a,b,c)=map(float, string.split(i))
		xlen=math.sqrt(a**2+b**2+c**2)
		if xlen > amax:
			amax=xlen
		if xlen < amin:
			amin=xlen
		xlen=xlen*scale
		if (xlen < vmax and xlen > vmin):
			print >> ofd, x, y, z, a*scale, b*scale, c*scale

	print >> sys.stderr, "max: %f (%f)" % (amax, amax*scale)

	fd.close()
	ofd.close()

main()
