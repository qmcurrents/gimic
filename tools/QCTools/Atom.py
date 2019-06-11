import string
import math
from . import Elements
from . import Const

_cfact={'au2m' : Const.au2m, 'au2nm' : Const.au2nm, 
	'au2a' : Const.au2a, 'au2pm' : Const.au2pm, 
	'None' : 1, 'none' : 1, 'a2au' : Const.a2au, 'nm2au' : Const.nm2au,
	'pm2au' : Const.pm2au}

class Atom:
	coord=(0,0,0)
	conv=_cfact['None']
	
	def __init__(self, coord=(0,0,0), sym='None', cf='None'):
		self.coord=coord
		self.conv=_cfact[cf]
		self._setelement(sym)
	
	def __add__(self, other):
		if hasattr(other, 'coord'): # is an Atom
			c1=self.coord[0]+other.coord[0]
			c2=self.coord[1]+other.coord[1]
			c3=self.coord[2]+other.coord[2]
			if self.element.symbol == other.element.symbol:
				return Atom((c1, c2, c3), self.element.symbol)
			else:
				return Atom((c1, c2, c3), 'None')
		else:
			c1=self.coord[0]+other[0]
			c2=self.coord[1]+other[1]
			c3=self.coord[2]+other[2]
			return (c1, c2, c3)
	
	def __sub__(self, other):
		if hasattr(other, 'coord'): # is an Atom
			c1=self.coord[0]-other.coord[0]
			c2=self.coord[1]-other.coord[1]
			c3=self.coord[2]-other.coord[2]
			if self.element.symbol == other.element.symbol:
				return Atom((c1, c2, c3), self.symbol)
			else:
				return Atom((c1, c2, c3), 'None')
		else:
			c1=self.coord[0]-other[0]
			c2=self.coord[1]-other[1]
			c3=self.coord[2]-other[2]
			return (c1, c2, c3)
	
	def __mod__(self, other):
		"Returns bond length between two Atoms."
		if hasattr(other, 'coord'): # is an Atom
			c1=self.coord[0]-other.coord[0]
			c2=self.coord[1]-other.coord[1]
			c3=self.coord[2]-other.coord[2]
		else: # is a tuple
			c1=self.coord[0]-other[0]
			c2=self.coord[1]-other[1]
			c3=self.coord[2]-other[2]
		return self.conv*math.sqrt(c1**2+c2**2+c3**2)
			
	def __repr__(self):
		s='{Atom: ' + self.element.name + ' at ' + str(self.coord) + '}'
		return s
		
	def _setelement(self, s):
		s=str.lower(s)
		#s=string.capitalize(s)
		try:
			self.element=Elements.PeriodicTable[s]
		except KeyError as x:
			print(("Invalid element:", x));

	def setconv(self, conv='None'):
		self.conv=_cfact[conv]
		
class TurboAtom(Atom):
	def __init__(self, str="0 0 0 q"):
		"Initialze an Atom from a line from a TM coord file"
		str=str.split()
		#atof=float(str)
		#print(atof[0])
		self.coord=(float(str[0]), float(str[1]), float(str[2]))
		self._setelement(str[3])

class XYZAtom(Atom):
	def __init__(self, str="x 0 0 0"):
		"Initialze an Atom from a line from a xyz file"
		str=string.split(str)
		atof=string.atof
		self.coord=(atof(str[1]), atof(str[2]), atof(str[3]))
		self._setelement(str[0])

