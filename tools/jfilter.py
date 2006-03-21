#!/usr/bin/env python

import string, math, sys, re

tol=3.0
mint=1.e-2
amin=0.e10
amax=0.e0
scale=1.e0
foo=[3.0]

defvals="""tol=3.0
mint=1.e-2
amin=0.e10
amax=0.e0
scale=1.e0
foo=[3.0]"""

try:
	execfile('jfilter.inp')
except:
	fd=open('jfilter.inp','w')
	fd.write(defvals)
	fd.close()
	pass

def main():
    	jvecf=sys.argv[1]
	fd=open(jvecf,'r')
	while 1: 
		t=fd.tell()
	    	l=fd.readline()
		if not re.search(r'^\s*#', l):
		    fd.seek(t)
		    break
	l=fd.readlines()
	global amax, amin

	print amax, amin, scale

	fds=[]

	for i in range(len(foo)):
		fds.append(open(jvecf+'.'+str(i+1),'w'))

	for i in l:
		(x,y,z,a,b,c)=map(float, string.split(i))
		xlen=math.sqrt(a**2+b**2+c**2)
		if xlen > amax:
			amax=xlen
		if xlen < amin:
			amin=xlen
		if (xlen > tol or xlen < mint):
			a=0.e0; b=0.e0; c=0.e0 
			print >> fds[0], x, y, z, a, b, c
			continue
		ff=1
		for j in foo: 
			if xlen <= j:
				#n=ff/scale
				#n=j/scale
				n=1/scale
				print >> fds[ff-1], x, y, z, a/n, b/n, c/n
				break
			ff+=1

	print >> sys.stderr, "%f (%f)" % (amax, amax*scale)

	for i in fds:
		i.close()

main()
