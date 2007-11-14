#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim:syntax=python
#
# getkw -- a simple input parser for Fortran 95
#
# Written by Jonas Juselius <jonas.juselius@chem.uit.no> 
# University of Tromsø, 2006
#
# TODO: 
#       o general cleanup
#       o python interface
#
# Known bugs: names with '-' mess things up, either in py or f90...
#

import pdb
import sys,os,inspect
import re, string, copy
from copy import deepcopy
from pyparsing import \
	Literal, Word, ZeroOrMore, Group, Dict, Optional, removeQuotes, \
	printables, ParseException, restOfLine, alphas, alphanums, nums, \
	pythonStyleComment, oneOf, quotedString, SkipTo, Forward, \
	commaSeparatedList, OneOrMore, Combine, srange, delimitedList, \
	downcaseTokens, line, lineno

verbose=True
strict=True

ival=re.compile(r'[-]?\d+$')
dval=re.compile(r'-?\d+\.\d*([dDeE][+-]?\d+)?',re.I)
lval=re.compile(r'(0|1|yes|no|true|false|on|off)$',re.I)


class Section:
	def __init__(self,name,arg=None,req=False,multi=False, callback=None):
		self.name=name
		self.sect={}
		self.kw={}
		self.arg=None
		self.req=req
		self.multi=multi
		self.set=False
		self.callback=callback
		if arg is not None:
			self.set_kwarg(arg)

	def validate(self,path=None):
		dlm=''
		if path is None:
			path=''
		else:
			path=path+dlm+self.name
			dlm='.'
		if self.req and not self.set:
			print '>>> Required section not set: %s \n' % (path)
			sys.exit(0)
		if self.arg is not None:
			self.arg.validate(path)
		for i in self.kw:
			if len(self.kw[i]) > 1 and not self.kw[i][0].multdef():
				print '>>> Mulitply defined key: %s ' % (path+'.'+i)
				if strict:
					sys.exit(1)
			for j in self.kw[i]:
				j.validate(path)
		for i in self.sect:
			if len(self.sect[i]) > 1 and not self.sect[i].multdef():
				print '>>> Mulitply defined section: %s ' % (path+'.'+i)
				if strict:
					sys.exit(1)
			for j in self.sect[i]:
				j.validate(path)

	def __cmp__(self, other):
		return cmp(self.name,other.name)

	def __getitem__(self, key):
		if self.sect.has_key(key):
			foo=self.sect
		elif self.kw.has_key(key):
			foo=self.kw
		else:
			return None
		return foo[key]

	def is_set(self, key=None):
		if key is None:
			return self.set
		if self.kw.has_key(key):
			return self.kw[key][0].is_set()
		if self.sect.has_key(key):
			return self.sect[key][0].is_set()

	def required(self):
		return self.req

	def multdef(self):
		return self.multi

	def arg_required(self):
		if arg is None:
			return False
		return self.arg.required()

	def add_sect(self, sect, set=False):
		if not self.sect.has_key(sect.name):
			self.sect[sect.name]=[]
		sect.set=set
		self.sect[sect.name].append(sect)

	def add_kwkw(self, kw, set=False):
		if not self.kw.has_key(kw.name):
			self.kw[kw.name]=[]
		kw.set=set
		self.kw[kw.name].append(kw)

	def add_kw(self, name, typ, arg, req=False, multi=False, 
			set=False, callback=None):
		if not self.kw.has_key(name):
			self.kw[name]=[]
		kw=Keyword(name,typ,arg,req,multi,callback)
		kw.set=set
		self.kw[name].append(kw)

	def set_kwarg(self, kw, set=False):
		if isinstance(kw,Keyword):
			kw.set=set
			self.arg=kw
		else:
			raise TypeError

	def set_arg(self, typ, arg, req=False, set=False):
		kw=Keyword(self.name,typ,arg,req)
		kw.set=set
		self.arg=kw

	def findkw(self, name):
		if self.kw.has_key(name):
			return self.kw[name][0]
		return None

	def getkw(self, name):
		if self.kw.has_key(name):
			return self.kw[name][0].arg
		return None

	def setkw(self, name, arg):
		if self.kw.has_key(name):
			self.kw[name][0].setkw(arg)
		else:
			print 'Error: invalid kw: ', name

	def findsect(self, name):
		if self.sect.has_key(name):
			return self.sect[name][0]
		return None

	def getsect(self, name):
		if self.sect.has_key(name):
			return self.sect[name][0]
		return None

	def get_keys(self):
		return self.kw

	def get_sects(self):
		return self.sect

	def get_arg(self):
		return self.arg

	def status(self):
		return self.set
	
	def set_status(self, set):
		if set:
			self.set=True
		else:
			self.set=False

	def sanitize(self, templ):
		self.equalize(templ)
		self.xvalidate(templ)

	# add missing keys
	def equalize(self, templ):
		for i in templ.kw:
			if not self.kw.has_key(i):
				self.kw[i]=templ.kw[i]
		for i in templ.sect:
			if not self.sect.has_key(i):
				self.sect[i]=templ.sect[i]
			for j in self.sect[i]:
				j.equalize(templ.sect[i][0])

	def run_callbacks(self, templ):
		if templ.callback is not None:
			templ.callback(self)
		for i in templ.kw:
			cb=templ.kw[i][0]
			if cb.callback is not None:
				cb.callback(self.kw[i])
		for i in templ.sect:
			for j in self.sect[i]:
				j.run_callbacks(templ.sect[i][0])

	#cross-validate against a template
	def xvalidate(self,templ,path=None):
		dlm=''
		if path is None:
			path=''
		else:
			path=path+dlm+self.name
			dlm='.'
		if templ.req and not self.set:
			print '>>> Required section not set: %s \n' % path
			sys.exit(1)
		self.xvalidate_arg(templ.arg,path)
		for i in self.kw:
			j=templ.findkw(i) 
			if j is None:
				print '>>> Invalid keyword: %s ' % (path+dlm+i)
				sys.exit(1)
			if len(self.kw[i]) > 1 and not j.multdef():
				print '>>> Mulitply defined key: %s ' % (path+dlm+i)
				if strict:
					sys.exit(1)
			for k in self.kw[i]:
				k.xvalidate(j,path)
		for i in self.sect:
			j=templ.findsect(i) 
			if j is None:
				print '>>> Invalid section: %s ' % (path+dlm+i)
				sys.exit(1)
			if len(self.sect[i]) > 1 and not j.multdef():
				print '>> Multiply defined section: %s ' % (path+dlm+i)
				if strict:
					sys.exit(1)
			for k in self.sect[i]:
				k.xvalidate(j,path)
				
		
	def xvalidate_arg(self,templ,path):
		kw=self.arg
		if templ is not None:
			if kw is not None:
				kw.xvalidate(templ)
			elif templ.req:
				print '>>> Required arg not set: %s \n' % (path+
						'.'+kw.name)
				sys.exit(1) # return state...
		elif kw is not None:
			print '>>> Invalid argument: %s ' % (path+'.'+kw.name)
			sys.exit(1)

	def __str__(self):
		nsect=0
		for i in self.sect:
			nsect=nsect+len(self.sect[i])
		nkw=0
		for i in self.kw:
			nkw=nkw+len(self.kw[i])
		s="SECT %s %d %s\n" % (self.name, nsect, self.set)
		if self.arg is not None:
			s=s+"ARG T KW %d\n" % (nkw)
			s=s+str(self.arg)
		else:
			s=s+"ARG F KW %d\n" % (nkw)

		for i in self.kw.values():
			for j in i:
				s=s+str(j)
		for i in self.sect.values():
			for j in i:
				s=s+str(j)
		return s
		
class Keyword:
	yes=re.compile(r'(1|yes|true|on)$',re.I)
	no=re.compile(r'(0|no|false|off)$',re.I)

	def __init__(self, name, typ, arg, req=False, multi=False, callback=None):
		self.name=name
		self.type=typ
		self.req=req
		self.multi=multi
		self.nargs=None
		self.arg=[]
		self.callback=callback
		
		if arg is None: # unlimited arg length
			self.nargs=-1
			arg=None
		elif isinstance(arg,int): # number of elements in self.arg == arg
			self.nargs=arg  
			arg=None
		if isinstance(arg,tuple) or isinstance(arg,list):
			self.nargs=len(arg)
		self.setkw(arg)
		self.set=False # reset the self.set flag

	def __cmp__(self, other):
		return cmp(self.name,other.name)
	
	def setkw(self, arg):
		if isinstance(arg,tuple):
			pass
		elif isinstance(arg,list):
			arg=tuple(arg)
		else:
			arg=(arg,)

		if arg[0] is None: # init stage
			self.arg=[]
			return

		if self.nargs > 0: 
			if len(arg) != self.nargs:
				print "keyword lenght mismatsh %s(%i): %i" % (self.name,
						self.nargs, len(arg))
				sys.exit(1)
		# store everyting as strings internally
		self.arg=[]
		for i in arg:
			self.arg.append(str(i))
		try:
			self.sanity()
		except TypeError:
			print 'Invalid argument:', arg
			sys.exit(1)
		self.set=True
	
	def is_set(self):
		return self.set

	def validate(self,path):
		if path is None:
			path=self.name
		else:
			path=path+'.'+self.name
		if self.req  and not self.set:
			print '>>> Required key not set: %s \n' % (path)
			if strict:
				sys.exit(1)

	def xvalidate(self,templ,path=None):
		if path is None or path == '':
			path=self.name
		else:
			path=path+'.'+self.name
		if templ.req and not self.set:
			print '>>> Required key not set: %s \n' % (path)
			sys.exit(1)
		if templ.type != self.type:
			print '>>> Invalid data type in: %s \n' % (path)
			sys.exit(1)
		if self.nargs < 0:  #  < 0 == unlimited arg length
#            print self.nargs, self.arg
			if len(templ.arg) != len(self.arg):
				print '>>> Invalid data length in: %s \n' % (path)
				sys.exit(1)
		return True
				
	def sanity(self):
		if self.arg[0] == 'None':
			return True
		if self.type == None:
			return True
		if (self.type == 'INT' or self.type == 'INT_ARRAY'):
			for i in self.arg:
				if not ival.match(i):
					print 'getkw: Not an integer: ', self.name, '=',i
					raise TypeError
		elif (self.type == 'DBL' or self.type == 'DBL_ARRAY'):
			for i in self.arg:
				if not dval.match(i):
					print 'getkw: Not a real: ', self.name, '=',i
					raise TypeError
		elif (self.type == 'BOOL' or self.type == 'BOOL_ARRAY'):
			for i in range(len(self.arg)):
				if not lval.match(self.arg[i]):
					print 'getkw: Not a bool: ', self.name, '=',self.arg[i]
					raise TypeError
				if self.yes.match(self.arg[i]):
					self.arg[i]='True'
				elif self.no.match(self.arg[i]):
					self.arg[i]='False'
		elif self.type == 'STR':
			tmp=self.arg
			self.arg=[]
			for i in tmp:
				if len(i.strip()) != 0:
					self.arg.append(i)
			return True
		else:
			print 'getkw: Unknown type: ', self.name, '=', self.type
			raise TypeError
		return True

	def multdef(self):
		return self.multi

	def required(self):
		return self.req

	def status(self):
		return self.set
	
	def set_status(self, set):
		if set:
			self.set=True
		else:
			self.set=False

	def __str__(self):
		if self.type == 'STR' and self.arg == []: # empty string
			nargs=-1 # flags as empty for Fortran code
		else:
			nargs=len(self.arg)
		s="%s %s %d %s\n" % (self.type, self.name, nargs, self.set)
		if self.arg != []:
			for i in self.arg:
				s=s+str(i)+'\n'
		return s

class GetkwParser:
	bnf=None
	caseless=False

	def __init__(self,templ=None):
		self.top=Section('top')
		self.stack=[self.top]
		self.cur=self.stack[0]
		self.templ=templ
		self.strg=None
		self.loc=None
		if templ is not None:
			self.path=[self.templ]
		else:
			self.path=None
		if GetkwParser.bnf == None:
			GetkwParser.bnf=self.getkw_bnf()
		self.parseString=self.bnf.parseString
	
	def set_caseless(self, arg):
		if arg is True:
			self.caseless=True
		else:
			self.caseless=False
	
	def getkw(self, path):
		pass

	def parseFile(self,fil):
		self.bnf.parseFile(fil)
		return self.top

	def add_sect(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		arg=None
		if len(q) > 1:
			if len(q[1]) > 0:
				arg=q[1][0]
		if self.caseless:
			name=name.lower()
		k=Section(name)
		self.cur.add_sect(k, set=True)  
		self.push_sect(k)

		if arg is not None:
			if self.templ is None:
				argt=self.guess_type(arg) 
			else:
				argt=self.check_type(arg, self.path[-1].arg.type)
			kw=Keyword(name, argt, arg)
		else:
			kw=None
		if kw is not None:
			k.set_kwarg(kw, set=True)

	def add_vecsect(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		arg=None
		if len(q) > 1:
			if len(q[1]) > 0:
				arg=q[1][0]
		if self.caseless:
			name=name.lower()
		k=Section(name)
		self.cur.add_sect(k, set=True)  
		self.push_sect(k)
		if arg is not None:
			if self.templ is None:
				argt=self.guess_type(arg) 
			else:
				argt=self.check_vectype(arg, self.path[-1].arg.type)
			kw=Keyword(name, argt, arg)
		else:
			kw=None
		if kw is not None:
			k.set_kwarg(kw, set=True)

	def push_sect(self,k):
		self.stack.append(k)
		self.cur=self.stack[-1]
		if self.templ is not None:
			x=self.path[-1].findsect(k.name)
			if x is None:
				print "Invalid section on line %d: \n%s" % (
						lineno(self.loc,self.strg), line(self.loc,self.strg))
				sys.exit(1)
			self.path.append(x)

	def pop_sect(self,s,l,t):
		if self.templ is not None:  
#            if self.path[-1].callback is not none:
#                self.path[-1].callback(self.stack[-1])
			del self.path[-1]
		del self.stack[-1]
		self.cur=self.stack[-1]
	
	def store_key(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		arg=q[1]
		if self.caseless:
			name=name.lower()
		if self.templ is None:
			argt=self.guess_type(arg)
		else:
			k=self.path[-1].findkw(name)
			if k is None:
				print "Unknown keyword '%s' line: %d" % (name, 
						lineno(self.loc,self.strg))
				if strict:
					sys.exit(1)
				argt=None
			else:
				argt=self.check_type(arg,k.type)
		k=Keyword(name,argt,(arg,))
		self.cur.add_kwkw(k,set=True)

	def store_vector(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		arg=q[1:]
		if self.caseless:
			name=name.lower()
		if self.templ is None:
			argt=self.guess_vectype(arg)
		else:
			k=self.path[-1].findkw(name)
			if k is None:
				print "Unknown keyword '%s' line: %d" % (name, 
						lineno(self.loc,self.strg))
				if strict:
					sys.exit(1)
				argt=None
			else:
				argt=self.check_vectype(arg,k.type)
		k=Keyword(name,argt,arg)
		self.cur.add_kwkw(k,set=True)

	def store_data(self,s,l,t):
		name=t[0]
		self.strg=s
		self.loc=l
		if self.caseless:
			name=name.lower()
		dat=t[1].split('\n')
		arg=[]
		for i in dat:
			arg.append(i.strip())
		k=Keyword(name,'STR',arg)
		self.cur.add_kwkw(k,set=True)

	def check_vectype(self, arg, argt):
		if argt == 'INT_ARRAY':
			for i in arg:
				if not ival.match(i):
					print 'Invalid type on line %d: Not an int: \n -> %s' % (
							lineno(self.loc,self.strg), line(self.loc,
								self.strg).strip())
					sys.exit(1)
		elif argt == 'DBL_ARRAY':
			for i in arg:
				if not dval.match(i):
					print 'Invalid type on line %d: Not a float: \n -> %s' % (
							lineno(self.loc,self.strg), line(self.loc,
								self.strg).strip())
					sys.exit(1)
		elif argt == 'BOOL_ARRAY':
			for i in arg:
				if not lval.match(i):
					print 'Invalid type on line %d: Not a bool: \n -> %s' % (
							lineno(self.loc,self.strg), line(self.loc,
								self.strg.strip()))
					sys.exit(1)
		elif argt != 'STR':
			print 'Invalid type on line %d: Not a %s: \n -> %s' % (
					lineno(self.loc,self.strg), argt, line(self.loc,
						self.strg).strip())
			sys.exit(1)
		return argt

	def check_type(self, arg, argt):
		if argt == 'INT':
			if not ival.match(arg):
				print 'Invalid type on line %d: Not an int: \n -> %s' % (
						lineno(self.loc,self.strg), line(self.loc,
							self.strg).strip())
				sys.exit(1)
		elif argt == 'DBL':
			if not dval.match(arg):
				print 'Invalid type on line %d: Not a float: \n -> %s' % (
						lineno(self.loc,self.strg), line(self.loc,
							self.strg).strip())
				sys.exit(1)
		elif argt == 'BOOL':
			if not lval.match(arg):
				print 'Invalid type on line %d: Not a bool: \n -> %s' % (
						lineno(self.loc,self.strg), line(self.loc,
							self.strg).strip())
				sys.exit(1)
		elif argt != 'STR':
			print 'Invalid type on line %d: Not a %s: \n -> %s' % (
					lineno(self.loc,self.strg), argt, line(self.loc,
						self.strg).strip())
			sys.exit(1)
		return argt

	def guess_vectype(self,arg):
		if ival.match(arg[0]):
			type='INT_ARRAY'
			for i in arg:
				if not ival.match(i):
					print 'invalid ', type
					sys.exit(1)
		elif dval.match(arg[0]):
			type='DBL_ARRAY'
			for i in arg:
				if not dval.match(i):
					print 'invalid ', type
					sys.exit(1)
		elif lval.match(arg[0]):
			type='BOOL_ARRAY'
			for i in arg:
				if not lval.match(i):
					print 'invalid ', type
					sys.exit(1)
		else:
			type='STR'
		return type
	
	def guess_type(self,arg):
		if ival.match(arg):
			return 'INT'
		if dval.match(arg):
			return 'DBL'
		if lval.match(arg):
			return 'BOOL'
		return 'STR'

	def getkw_bnf(self):
		lcb  = Literal("{").suppress()
		rcb  = Literal("}").suppress()
		lsb  = Literal("[").suppress()
		rsb  = Literal("]").suppress()
		lps  = Literal("(").suppress()
		rps  = Literal(")").suppress()
		eql  = Literal("=").suppress()
		dmark=Literal('$').suppress()
		end_sect=rcb
		end_data=Literal('$end').suppress()
		prtable = srange("[0-9a-zA-Z]")+r'!$%&*+-./<>?@^_|~'

		kstr=Word(prtable) ^ quotedString.setParseAction(removeQuotes)

		name = Word(alphas+"_",alphanums+"_")
		
		vec=lsb+delimitedList(Word(prtable) ^ Literal("\n").suppress() ^\
				quotedString.setParseAction(removeQuotes))+rsb
		key=kstr ^ vec
		keyword = name + eql + kstr
		vector = name + eql + vec
		data=Combine(dmark+name)+SkipTo(end_data)+end_data
		data.setParseAction(self.store_data)
		sect=name+lcb
		sect.setParseAction(self.add_sect)
		key_sect=name+Group(lps+kstr+rps)+lcb
		key_sect.setParseAction(self.add_sect)
		vec_sect=name+Group(lps+vec+rps)+lcb
		vec_sect.setParseAction(self.add_vecsect)
		end_sect.setParseAction(self.pop_sect)

		keyword.setParseAction(self.store_key)
		vector.setParseAction(self.store_vector)

		section=Forward()
		input=section ^ data ^ keyword ^ vector

		sectdef=sect ^ key_sect ^ vec_sect
		section << sectdef+ZeroOrMore(input)+rcb
		
		bnf=ZeroOrMore(input)
		bnf.ignore( pythonStyleComment )

		return bnf


def test( strng ):
	bnf = GetkwParser()
	try:
		tokens=bnf.parseString(strng)
	except ParseException, err:
		print err.line
		print " "*(err.column-1) + "^"
		print err
	print bnf.top
	return tokens

if __name__ == '__main__':
	teststr="""
title = foo
title2="fooo bar"

perturbation (apa) { 
foo=[1, 2, 3,
4,5, 6,7,8,9, 
10] 
bar=22.0
}

dalron {
	foo=1 
	uar=22 
}

verbose=yes #(yes|true|on|1)

DAlron2 {
	foo=1
	bar=1

	
	vafan(3) {
		foo=0 
	}

	$COORD
	  o 0.0 0.0 0.0
	  h 1.0 1.0 0.0
	  h -1.0 1.0 0.0
	$end
}

test(off) {
	zip=" asfasf adsf asf asfd "
	$as
		"asdf sadwa:
		asf  ddddddddddddddddddddddd"

		sdfff
		sdff
	$end
}

"""
	ini = test(teststr)
	print ini

