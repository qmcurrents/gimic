import re
import string
import sys
import Atom

def readAtoms(file, fmt=None):
	fmts={'coord' : readCoord, 'xyz' : readXYZ}
	try:
		func=fmts[fmt]
	except KeyError, x:
		print "Unknown format", x
		sys.exit()
	
	return func(file)

def readCoord(file):
	"Read TURBOMOLE coord file, return list of Atoms"
	
	atoms=[Atom.Atom((0,0,0), 'None')]
	
	try:
		f=open(file, 'r')
	except IOError, x:
		print '%s: %s' % (x[1], infile)
		sys.exit(1)
	
	blankline=re.compile('[ \t]*$')
	dollar=re.compile('\$')
	
	while 1:
		line=f.readline()
		
		if not line: break # EOF
		
		if blankline.match(line):
			continue
		elif dollar.match(line):  			# starts with $
			if re.match("\$coord", line): 	# $coord -> continue
				continue
			else: break # end of coords
		else:
			try:
				atoms.append(Atom.TurboAtom(line))
			except:
				print "Error in file", file
				sys.exit()

	f.close()
	return atoms

def readXYZ(file):
	"Read XYZ file, return list of Atoms"
	
	atoms=[]
	atoms=[Atom.Atom((0,0,0), 'None')]
	
	try:
		f=open(file, 'r')
	except IOError, x:
		print '%s: %s' % (x[1], infile)
		sys.exit(1)
	
	blankline=re.compile('[ \t]*$')

	# Parse beginning of file
	while 1:
		line=f.readline()
		
		if not line: break # EOF
		
		if blankline.match(line):
			continue
		else:
	# First line is numer of atoms
			try:
				n=string.atoi(line)
			except:
				print "Error in file", file
				sys.exit()
			f.readline()
			break
		
	for i in range(1,n):
		line=f.readline()
		
		if not line: break # EOF
		
		if blankline.match(line):
			i=i-1
			continue
		else:
			try:
				atoms.append(Atom.XYZAtom(line))
			except:
				print "Error in file", file
				sys.exit()
	f.close()	
	return atoms

def readData(file, cols=(1, 2)):
	"Read data file, return list of tuples"
	
	data=[]
	
	
	try:
		f=open(file, 'r')
	except IOError, x:
		print '%s: %s' % (x[1], infile)
		sys.exit(1)
	
	blankline=re.compile('[ \t]*$')
	comment=re.compile('[ \t]*#.*$')
	float=re.compile('\.')
	
	while 1:
		line=f.readline()
		
		if not line: break # EOF
		
		if blankline.match(line):
			continue
		elif comment.match(line):  			# starts with #
			continue
		else:
			tmpdata=[]
			line=string.split(line)
			for i in cols:
				if float.search(line[i-1]):	
					try:
						tmpdata.append(string.atof(line[i-1]))
					except:
						print "Error in data!"
				else:
					try:
						tmpdata.append(string.atoi(line[i-1]))
					except:
						print "Error in data!"

			data.append(tmpdata)
	f.close()
	return data

def basename(name):
	"Split name at first '.'"
	name=string.strip(name)
	name=string.split(name, '.')
	return name[0]
