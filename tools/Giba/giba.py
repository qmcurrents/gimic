#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Giba: Gimic Bond Analysis
#       
#       Copyright 2011 Raul Mera <rmera@chem.helsinki.fi>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#       
#       
#A little thing to manage all the calibrations needed for bond strenght analysis with GIMIC.
#

import sys, os, shutil, heapq
from optparse import OptionParser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def setgimicfiles(namedir,command):
	if command.startswith("qsub") or command.startswith("sbatch"):
		shutil.copy("./"+command.split()[1],namedir)
#	os.system("ln -s ./XDENS "+namedir) #not nice but 
#	os.system("ln -s ./MOL "+namedir)
	shutil.copy("./XDENS",namedir)
	shutil.copy("./MOL",namedir)
	

#makes sample slices for with or 
def slices(woh,step,job,fake):
	llimit=-5.0
	rlimit=llimit+step
	cont=0
	while(rlimit<=5.0):
		founddistance=False
		namedir=woh+"_%04d_%3.1fT%3.1f"%(cont,llimit,rlimit)
		os.mkdir(namedir)
		setgimicfiles(namedir,command)
		os.chdir(namedir)
		fin=open("../gimic.inp","r")
		fout=open("gimic.inp","w")
		for i in fin:
			if i.lstrip(" ").lstrip("\n").startswith("width") and woh=="w":
				fout.write("      width=[%3.1f,%3.1f]       \n"%(llimit,rlimit))
				founddistance=True
			elif i.lstrip(" ").lstrip("\n").startswith("height") and woh=="h":
				fout.write("      height=[%3.1f,%3.1f]       \n"%(llimit,rlimit))
				founddistance=True		
			else:
				fout.write(i)
		if not founddistance:
			print  "no distance found in file!!"
			sys.exit(1)
		fout.close()
		fin.close()
		if not fake:
			os.system(job) 
		#os.system("gimic > gimic.out")
		llimit+=step
		rlimit=llimit+step
		cont+=1
		os.chdir("../")	
	return True 	
			
#leng: lenght of the bond, nsteps: number of steps in which the bond will be divided.
#This creates inputs for gimic jobs at different steps along a bond.
def distances(leng,step,job,fake):
	distance=0.0
	print ""
	print "DISTANCES   Total: ", distance, "Total steps:", nsteps, "Step lenght:", step
	print""	
	cont=0
	while(distance<=leng):
		founddistance=False
		namedir="d_%04.3f"%distance
		os.mkdir(namedir)
		setgimicfiles(namedir,command)
		os.chdir(namedir)
		fin=open("../gimic.inp","r")
		fout=open("gimic.inp","w")
		for i in fin:
			if i.lstrip(" ").lstrip("\n").startswith("distance"):
				fout.write("        distance=%7.6f        # place grid 'distance' between atoms\n"%distance)
				founddistance=True
			else:
				fout.write(i)
		if not founddistance:
			print  "no distance found in file!!"
			sys.exit(1)
		fout.close()
		fin.close()
		if not fake:
			os.system(job)
		#os.system("gimic > gimic.out")
		distance+=step
		cont+=cont
		os.chdir("../")
	return True

#gets currents or mod currents in a.u. from a file called "gimic.out" which must exist in the current directory
#if positive or negative contributions differ by more than the given threshold, raises ValueError exception.
#returns both positive and negative values for the currents.
def get_current(curr=False,threshold=0.001):
	ReadingCurrents=False
	positive=0.0
	negative=0.0
	f=open("gimic.out","r")
	for i in f:
		if not curr:
			if "Induced mod current (au)" in i:
				ReadingCurrents=True
		else:
			if "Induced current (au)" in i:
				ReadingCurrents=True
		if ReadingCurrents and (" Positive contribution:" in i):
			positive=float(i.split()[2])
		if ReadingCurrents and (" Negative contribution:" in i):
			negative=float(i.split()[2])
			break
	if positive+negative>=threshold:
		passed=positive
	else:
		passed=0
	return positive,negative,passed
	
#plots positive and negative currents against something. labelx is the label for that something			
def plot_currents(x,ypos,yneg,yerr,labelx):
	plt.figure(1)
	plt.subplot(211)
	plt.plot(x,ypos,"r")
	plt.plot(x,yerr,"ko")
	plt.title("Currents Vs "+labelx)
	plt.xlabel(labelx+" (a.u.)")
	plt.ylabel("Positive currents (a.u.)")
	plt.subplot(212)
	plt.plot(x,yneg,"b")
	yerr2=-1*yerr
	plt.plot(x,yerr2,"ko")
	plt.xlabel(labelx+" (a.u.)")
	plt.ylabel("Negative currents (a.u).")
	labelx.rstrip("(radians)") #remove this part from the image name for the B-field scan
	name="".join(labelx.split())
	plt.savefig(name)
	return True

#Naive way of finding the local miniun of a list that is closest to the begining
def local_minimun(l):
	if l[0]<=l[1]:
		return 0
	for i in range(len(l))[1:-1]:
		if l[i]<l[i-1] and l[i]<l[i+1]:
			return i



#calculates a vector between paralel to the integration plane and to the plane of the molecule.
#(angle=0) and paralel to the inegration plane but perpendicular to the molecular plane (angle=pi/2)
def orient_magnetic_field(a1,a2,a3,angle):
	coord=open("coord","r")
	cont=0
	for i in coord:
		if cont==a1:
			temp=i.split()
			first=np.array([float(temp[0]),float(temp[1]),float(temp[2])])
		elif cont==a2:
			temp=i.split()
			second=np.array([float(temp[0]),float(temp[1]),float(temp[2])])
		elif cont==a3:
			temp=i.split()
			third=np.array([float(temp[0]),float(temp[1]),float(temp[2])])
		#note that a1, etc are atoms numbers starting from 1 and in the coord file there is one
		#unimportant line. Incrementing cont here compensates both facts.
		cont+=1 
	if angle<0:
		stderr.write("Angle must be positive, setting to positive\n")
		angle=abs(angle)
	while angle>(np.pi/2.0):
		sys.stderr.write("Angle must be between 0 and pi/2, setting to pi/2\n") # pi-angle")
		angle=np.pi/2.0
		#angle=np.pi-angle
	coeff=((angle*2.0)/np.pi)
	#print first, second, third
	coord.close()
	bond=second-first
	Tvector=np.cross(bond,(third-first))
	Pvector=np.cross(bond,Tvector)
	Finalvec=coeff*Tvector + (1-coeff)*Pvector  
	Finalvec=Finalvec/np.linalg.norm(Finalvec) #unitarize the vector
	vectorline="        magnet=[%4.2f, %4.2f, %4.2f]  \n"%(Finalvec[0],Finalvec[1],Finalvec[2])
	fin=open("gimic.inp","r")
	r=fin.readlines()
	fin.close()
	fout=open("gimic.inp","w")
	for i in r:
		if "magnet=[" in i:
			fout.write(vectorline)
		else:
			fout.write(i)
	fout.close()
	return Finalvec

	
				
#handles plotting of distances and slices obtained with this same program. In addition, it returns the distance
#at which (one of the) the minimun value for the positive contribution occurs (which in principle should be used 
#for calculations). In case of slices analyses, it returns 2 local minima: starting from the front and from the
#back of the distance sequence. It stop execution if the difference of positive and negative contributions 
#is larger than the given threshold.
#argument key indicates wether bond distances, width slices or height slices are being measured. Its value must
#be of "_d", "_w" or "_h" in each respective case.
def analize_distances(key,mod=True,threshold=0.001):
	positives=[] #positive currents
	negatives=[] 
	errors=[]  #zeroes when positive+negative<threshold, ==positive otherwise
	dists=[]
	founddir=False
	dirs=os.listdir("./")
	dirs.sort() #for plotting and for estimating the minima we want this in order, they will be from most negative to most positive
	for i in dirs:
		if i.startswith(key):
			founddir=True
			os.chdir(i)
			currents=get_current(mod,threshold) #currents=[positive,negative,correct_difference?]
			if currents[2]: #if error != 0
				sys.stderr.write("Difference between negative and positive contributions to current is larger than ")				  
				sys.stderr.write("current threshold of "+str(threshold)+" for file in "+i+"\n")
				print currents[0], currents[1], currents[0]+currents[1], currents[2]###############################
			positives.append(currents[0])
			negatives.append(currents[1])
			errors.append(currents[2])
			if key=="d_":
				dists.append(float(i.split("_")[1]))
			else:
				dists.append(i.split("_")[-1].split("T")) # i will be h_x.y-w.z, so the resulting is ["x.y","w.z"] 
				dists[-1]=[float(dists[-1][0]),float(dists[-1][1])] #they were strings
			os.chdir("../")
	print "Currents", positives, "\n", negatives, "\n", errors
	if not founddir:
		return 9999999
	parray=np.array(positives) #will be used for plotting
	narray=np.array(negatives)
	earray=np.array(errors)
	#masks for the errors array
	m = earray==0
	maskearray=np.ma.masked_array(earray,mask=m) #masked array, only the faulty values are unmasked
	if key == "d_":
		smallest=heapq.nsmallest(1,positives) #Look for the smallest value in the positive parts
		#print positives, smallest######################################################################################
		smallest=positives.index(smallest[0]) #
		dist_smallest=dists[smallest] #distance of smallest positive current
		darray=np.array(dists) #for plotting
#in the following case the curve should have the shape of a W so we need 2 minima, one from the 
#front and one from the back.
	else:
		dist_smallest=[]
		dist_smallest.append(dists[local_minimun(positives)]) #local min starting from the front		
		positives.reverse() #now we start from the back
		mini=local_minimun(positives) 
		#the following is required to transform the index in the reverse positive  
		#sequence into an index in the distance sequence
		mini=-1*(mini+1)
		dist_smallest.append(dists[mini])
		positives.reverse() #return the positives sequence to normal order 		
		darray=np.array(dists)
		darray=(darray[:,:1]+darray[:,1:])*0.5 #now darray contains the average of each pair of dists which will be plotted
		darray=np.reshape(darray,darray.shape[0]) #make one continuous array with all  the elements
	if key=="d_":
		xlabel="Bond dists"
	elif key=="w_":
		xlabel = "Width slices"
	else:
		xlabel = "Height slices"
	plot_currents(darray,parray,narray,maskearray,xlabel)
	return dist_smallest




#handles plotting of bfields obtained with this same program. In addition, it returns the angle
#at which (one of the) the minimun value for the positive contribution occurs (which in principle should be used 
#for calculations). 
def analize_field(mod=True,threshold=0.01):
	positives=[] #positive currents
	negatives=[] 
	errors=[]  #zeroes when positive+negative<threshold, ==positive otherwise
	dists=[]
	founddir=False
	dirs=os.listdir("./")
	dirs.sort() #for plotting and for estimating the minima we want this in order, they will be from most negative to most positive
	for i in dirs:
		if i.startswith("a_"):
			founddir=True
			os.chdir(i)
			currents=get_current(mod,threshold) #currents=[positive,negative,correct_difference?]
			if currents[2]: #this should work
				sys.stderr.write("Difference between negative and positive contributions to current is larger than ")				  
				sys.stderr.write("current threshold of "+str(threshold)+" for file in "+i+"\n")	
				print currents[0], currents[1], currents[0]+currents[1], currents[2]###############################
			positives.append(currents[0])
			negatives.append(currents[1])
			errors.append(currents[2])
			dists.append(float(i.split("_")[1]))
			os.chdir("../")
	if not founddir:
		return 9999999
	parray=np.array(positives) #will be used for plotting
	narray=np.array(negatives)
	earray=np.array(errors)
	#masks for the errors array
	m = earray==0
	maskearray=np.ma.masked_array(earray,mask=m) #masked array, only the faulty values are unmasked
	smallest=heapq.nsmallest(1,positives) #Look for the smallest value in the positive parts
	print positives, smallest######################################################################################
	smallest=positives.index(smallest[0]) #
	dist_smallest=dists[smallest] #distance of smallest positive current
	darray=np.array(dists) #for plotting
	darray=darray
	xlabel="B-field angle (radians)"
	plot_currents(darray,parray,narray,maskearray,xlabel)
	return dist_smallest

#handles plotting of bfields obtained with this same program. In addition, it returns the angle
#at which (one of the) the minimun value for the positive contribution occurs (which in principle should be used 
#for calculations).  Is a bit more readable than the old analyze_all
def analyze_all_pretty(mod, threshold, toanali,refatom):
	for j in toanali:
		if j=="a_":  #this one is handled rather differently
			smallest=analize_field(mod,threshold)
			if len(refatom)<3:
				continue
			if smallest<9999999:
				print "Smallest value for positive current found at "+str(smallest)+" radians "
				print "This value will be written to the gimic.inp file"
				print "refatom", len(refatom)
				orient_magnetic_field(refatom[0],refatom[1],refatom[2],smallest)
			continue
		if j=="w_":
			turn="Widthslices"
			search="width=["
		elif j=="h_":
			turn="Heightslices"
			search="height=["
		else:
			turn="Distances"
			search="distance="
		smallest=analize_distances(j,mod,threshold)
		if smallest==9999999: #most likely, the results were not there
			continue
		print "Analizing", turn
		if j!="d_": #analysis for plane width and height
			if abs(smallest[0][0])<abs(smallest[1][1]):
				smallest[1][1]=abs(smallest[0][0])
			else:
				smallest[0][0]=-1*(smallest[1][1])
			print "corner local minima for positive current found at", smallest[0][0], "and"
			print smallest[1][1],
			print "This will be written to gimic.inp"
			restline="%3.1f,%3.1f]    \n"%(smallest[0][0],smallest[1][1]) #line to be written to gimic.inp
		else:  ##analysis for bond distances
		 	print "minima for positive current found at", smallest,
			print "This will be written to gimic.inp"
		 	restline="%5.4f    \n"%smallest #line written to gimic.inp
		#common part, this just writes the obatined results to gimic.inp
		fin=open("gimic.inp","r")
		r=fin.readlines()
		fin.close()
		fout=open("gimic.inp","w")			
		for i in r:
			if search in i:
				fout.write("         "+search+restline)
			else:
				fout.write(i)
		fout.close()
	return True




def add_xyz(newlines,inxyz,outxyz):
	xyzin=open(inxyz,"r")
	xyzout=open(outxyz,"w")
	number=int(xyzin.readline())
	xyzin.readline() #the white line
	xyzout.write(str(number+len(newlines))+"\n\n") #number of atoms and the white line
	for i in xyzin:
		xyzout.write(i)
	for i in newlines:
		xyzout.write(i)
	return True


def scan_field(step, command, a1, a2, a3, fake):
	angle=0.0
	angleatoms=[]
	while(angle<=np.pi/2.0):
		namedir="a_%5.3f"%angle
		os.mkdir(namedir)
		shutil.copy("./gimic.inp", namedir)
		setgimicfiles(namedir,command)
		shutil.copy("coord",namedir) #additional file required
		os.chdir(namedir)
		vector=orient_magnetic_field(a1,a2,a3,angle)
		if not fake:
			os.system(command)
		if "integral.xyz" in os.listdir("../"): #help analyze the results
			fileout="integral%5.3f.xyz"%(angle)
			angleatoms.append("Ca     %4.2f   %4.2f   %4.2f\n"%(vector[0],vector[1],vector[2]))
			add_xyz([angleatoms[-1]],"../integral.xyz", fileout)		
		os.chdir("../")	
		angle=angle+step
	fileout="integral_allangles.xyz"
	add_xyz(angleatoms,"./integral.xyz", fileout)
		
		
	

#############################################################################################################
#############################################################################################################
#All the options and help stuff
#############################################################################################################
#############################################################################################################

usage="Orient B-field: "+sys.argv[0]+" -o atom1 atom2 atom3\n Other uses: "+sys.argv[0]+" [options]" 

parser = OptionParser(usage=usage)

parser.add_option("--fake", action="store_true", dest="fake", default=False,
				help="Set, but don't run the gimic calculations, not valid for analyses runs")

parser.add_option("-L", "--lazy", action="store_true",
				dest="lazy", default=False,
				help="asks only for three atoms and does everything automatically")

#This one is not really an option since is the default procedure if no other options are present.
parser.add_option("-A", "--analyze", action="store_true",
				dest="analyze", default=True,
				help="Analyze results from gimic calibration trials ran with this program. Default behavior")

parser.add_option("-F", "--field-angle", action="store_true",
				dest="fieldscan", default=False,
				help="Scans different angles for B-field")

parser.add_option("-D", "--distances", action="store_true",
				dest="rundistances", default=False,
				help="Set and execute calibration gimic runs setting integration plane in different points along a bond at constant plane size")

parser.add_option("-H", "--height", action="store_true",
				dest="runheights", default=False,
				help="Set and execute calibration runs in slices covering plane height from -5.0 to 5.0 a.u. at constant position and width")

parser.add_option("-W", "--width", action="store_true",
				dest="runwidth", default=False,
				help="Set and execute calibration runs in slices covering plane width from -5.0 to 5.0 a.u. at constant position and height")

parser.add_option("-O", "--orient-field", action="store_true",
				dest="orient_field", default=False,
				help="Set the B field in gimic input so is paralel to a molecular plane defined by 3 atoms")

parser.add_option("-c", "--currents", action="store_true",
				dest="module", default=False,
				help="use currents, not modules, in analyses (used with option -A). Default False")

parser.add_option("-t", "--threshold", action="store", type="float",
				dest="threshold", default=0.01,
				help="Maximun difference allowed between positive and negative contributions to the current in a.u. Used with option -A. Default 0.001")

parser.add_option("-l", "--bond-lenght", action="store", type="float",
				dest="bondleng", default=4.0,
				help="Lenght of the bond to scan, in a.u. Used with option -D. Default 4")

#parser.add_option("-n", "--number-of-steps", action="store", type="int",
#				dest="nsteps", default=30,
#				help="Number of steps in the bond scan if using -D or -F. Default 20")                 

parser.add_option("-s", "--step-lengh", action="store", type="float",
				dest="steplen", default=0.1,
				help="Lenght of the step for slice, distance andangle scans in a.u. or radians Used with options -H and -W. Default 0.1.")

parser.add_option("-a", "--angle", action="store", type="float",
				dest="angle", default=0.0,
				help="Angle of the final b field (radians).  Used with option -O. Default 0, paralell to the molecule plane.")

parser.add_option("-j", "--job-script", action="store",
				dest="jobscript", default="job.cmd",
				help="Jobscript for runing gimic calculations in the qsub system. Used with options -D, -H and -W. Default job.cmd")

parser.add_option("--job-system", action="store",
				dest="jobsystem", default="gimic",
				help="Type of queuing system. Supported: sbatch, qsub, gimic (i.e none). Default gimic")
				
parser.add_option("--group", action="store",
				dest="group", default="d_",
				help="To be used with option -A. Group to be analyzed: 'd_', 'a_', 'h_' or 'w_' for distance, angle, height or width, respectively")


(options, args) = parser.parse_args()

print 
print "For all the procedures in this program, a gimic.inp file",
print "must exist in the base directory. Each specific procedure may have other requirements."
print "use", sys.argv[0], "--help for help with the options and procedures" 


#
#End options and help
#

if options.jobsystem=="gimic":
	command="nohup gimic > gimic.out &"  #UNIX only
elif options.jobsystem=="sbatch":
	command="sbatch "+options.jobscript
	print command
elif options.jobsystem=="qsub":
	command="qsub "+options.jobscript
else: 
	print "Queing system not supported, will not run the executables" #better not to run than to run without queuing
	command=""
	options.fake=True

if options.rundistances:
	os.system("rm -R d_*")
	distances(options.bondleng, options.steplen, command,options.fake)

elif options.runheights:
	os.system("rm -R h_*")
	slices("h",options.steplen,command, options.fake)

elif options.runwidth:
	os.system("rm -R w_*")
	slices("w",options.steplen,command, options.fake)


elif options.fieldscan:
	os.system("rm -R a_*")
	scan_field(options.steplen,command,int(args[0]),int(args[1]),int(args[2]),options.fake)
	
elif options.orient_field:
	status=orient_magnetic_field(int(args[0]),int(args[1]),int(args[2]), options.angle)
	if status:
		print "The magnetic field has been set so it is perpendicular to the bond in study",
		print "(first two parameters) but paralell to the plane of the molecule as defined ",
		print "by the three atoms given. "



#This option doesnt work very well, especially, it doesn't work with 
#sbatch job system.
elif options.lazy:
	os.system("rm -R a_*") #delete previous results
	os.system("rm -R d_*")
	os.system("rm -R h_*")
	os.system("rm -R w_*")
	distances(options.bondleng, options.nsteps, options.jobscript,options.fake)
	analyze_all_pretty(options.module,options.threshold,("d_"),[0])
	slices("h",options.steplen,options.jobscript, options.fake)
	analyze_all_pretty(options.module,options.threshold,("h_"),[0])
	slices("w",options.steplen,options.jobscript, options.fake)
	analyze_all_pretty(options.module,options.threshold,("w_"),[0])
	scan_field((np.pi/(2*options.nsteps)),command,int(args[0]),int(args[1]),int(args[2]),options.fake)
	analyze_all_pretty(options.module,options.threshold,("a_"),(int(args[0]),int(args[1]),int(args[2])))
	print "The bond distance, plane size and field angle have been optimized. For extra safety", 
	print "run gibosa with option --lazy again."

elif options.analyze:
	if len(args)<3:
		refatoms=[0]
	else:
		print "Angles"
		refatoms=(int(args[0]),int(args[1]),int(args[2]))
	print options.group
	analyze_all_pretty(options.module,options.threshold,[options.group],refatoms)
	

	

