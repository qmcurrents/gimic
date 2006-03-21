#!/usr/bin/env python
#
# $Id: mkdep90.py 1.1 03/10/31 15:14:43-00:00 jonas@ $
#

import sys, string, re

use=re.compile('\s*use\s+(\w+)[^!,]*', re.I)
mod=re.compile('\s*module(?!\s+procedure)\s+(\w+)[^!,]*', re.I)
include=re.compile('\s*#?include\s+["<\']?([\w\._]+)[">\']?[^!,]*', re.I)

class depfile:
	def __init__(self, name):
		self.name=name
		self.oname=""
		self.uses={}
		self.includes={}
		ri=string.rindex(name, '.')
		self.oname=name[:ri]+'.o'
	
	def adduses(self, u):
		self.uses[u]=''

	def addincludes(self, u):
		self.includes[u]=''

def main():
	allmods={}
	dfiles=[]
	for ff in sys.argv[1:]:
		fd=open(ff, 'r')
		buf=fd.readlines()
		fd.close()
		dfiles.append(depfile(ff))

		for ln in buf:
			m=mod.match(ln)
			if m is not None:
				allmods[string.lower(m.group(1))]=ff
				continue

			m=use.match(ln)
			if m is not None:
				dfiles[-1].adduses(string.lower(m.group(1)))
				continue

			m=include.match(ln)
			if m is not None:
				dfiles[-1].addincludes(string.lower(m.group(1)))
				continue

	for df in dfiles:
		deps=[]
		for dd in df.uses.keys():
			try:
				if (allmods[dd] != df.name):
					ri=string.rindex(allmods[dd], '.')
					omod=allmods[dd][:ri]+'.o'
					deps.append(omod)
			except:
				print >> sys.stderr, 'Missing dependency for', dd, 'in',\
				df.name

#        for dd in df.includes.keys():
#                deps.append(dd)

		if deps:
			dstr=df.oname+': '
			for i in deps:
				dstr=dstr+i+' '
			print dstr


if __name__ == '__main__':
	main()
