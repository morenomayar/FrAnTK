
from __future__ import division
import math
import gzip
import sys
import os
import subprocess
import copy

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

transitionalleles=["CT", "GA", "TC", "AG"]

#convert freqs to treemix
#this will create an ALL file and a NT file
#required arguments: freqpref, tmpref, popsofint

allowedpars=["freqpref", "tmpref", "popsofint"]

USAGE="""
Convert freqs file to treemix input.
Arguments should be passed as: arg=value.
freqpref\tString. Prefix from prefix_freqs.gz file created by BuildFreqs.py. Required. 
tmpref\tString. Prefix for output file names. Should be different from plinkpref. Required. 
popsofint\tString. List of (at least 3) populations of interest for treemix file. This is one population per line. Required. 
Sample call: Freqs2Treemix.py freqpref=prefix tmpref=newprefix popsofint=poi
Run without arguments to get this message. 
"""

argvect=sys.argv[1:]
if len(argvect)==0:
	print(USAGE)
	sys.exit(1)

arghash={}
for i in argvect:
	if "=" not in i:
		print(i+" is not a valid arg sntx\n")
		print(USAGE)
		sys.exit(1)
	if len(i.split("="))!=2 or "" in i.split("="):
		print(i+" is not a valid arg, use arg=val")
		print(USAGE)
		sys.exit(1)
	argname=i.split("=")[0].split()[0]
	if argname not in allowedpars:
		print(i+" is not a valid arg\n")
		print(USAGE)
		sys.exit(1)
	arghash[argname]=i.split("=")[1].split()[0]

for i in allowedpars:
	if i not in arghash:
		print(i+" is missing\n")
		print(USAGE)
		sys.exit(1)

freqpref=arghash["freqpref"]
tmpref=arghash["tmpref"]
popsofint=arghash["popsofint"]

if freqpref == tmpref:
	print("freqpref and tmpref should be different")
	print(USAGE)
	sys.exit(1)

#read list of pops of interest

l=open(popsofint, "r").readlines()
poi=[]
for i in l:
	poi.append(i.split()[0])

#check that there are at least 3 pops

if len(poi) < 3:
	print("There should be at least 3 pops")
	sys.exit(1)

#read _pop file
l=open(freqpref+"_pop", "r").readlines()
nl=[]
for i in l:
	nl.append(i.split()[0])

#check that pops of interest are in _pop

for i in poi:
	if i not in nl:
		print(i+" not found in "+freqpref+"_pop")
		sys.exit(1)


#identify positions in freqs file based on _pop file
#NOTE that now everything will be with respect to this ordering
#pops in pop (freqs) file are alphabetical due to strat, but poi will change the order

p=[]
for j in poi:
	for i in xrange(len(nl)):
		if nl[i] == j:
			p.append(i)
			break

#open ALL and NT pos outfiles for writing

Oap=open(tmpref+"_ALL_tmpos", 'w')
Onp=open(tmpref+"_NT_tmpos", 'w')

#open ALL and NT treemix outfiles for writing (gzipped)

Oa=subprocess.Popen("gzip -c > "+tmpref+"_ALL_tm.gz", shell=True, stdin=subprocess.PIPE)
On=subprocess.Popen("gzip -c > "+tmpref+"_NT_tm.gz", shell=True, stdin=subprocess.PIPE)

#print(header for tm outs (note that it is space-delimited)

if pyver==2:
	Oa.stdin.write(" ".join(poi)+"\n")
	On.stdin.write(" ".join(poi)+"\n")
else:
	Oa.stdin.write((" ".join(poi)+"\n").encode())
	On.stdin.write((" ".join(poi)+"\n").encode())

##not necessary any more
##open ncop file
#fnc=subprocess.Popen(["zcat",freqpref+"_ncop.gz"],stdout=subprocess.PIPE).stdout

#for each line:
with subprocess.Popen(["zcat",freqpref+"_freqs.gz"],stdout=subprocess.PIPE).stdout as f:
	for line in f:
		#line=f.readline()
		##not necessary any more
		##read a line from ncop file for each freq line
		#linenc=fnc.readline()
		#turn nt flag off
		ntflag=0
		#get first member of triplet for all pops of interest
		if pyver==2:
			l=line.split()
		else:
			l=line.decode().split()
		#freqs=[]
		A1cops=[]
		A2cops=[]
		for i in p:
			#freqs.append(l[5+i*3])
			A1cops.append(l[5+i*3+1])
			A2cops.append(l[5+i*3+2])
		
		#if it is not a fully covered site, skip
		#if 'N' in freqs:
		if 'N' in A1cops:
			continue
		
		#if site is monomorphic, skip
		#if set(['0.0']) == set(freqs) or set(['1.0']) == set(freqs):
		if set(['0']) == set(A1cops) or set(['0']) == set(A2cops):
			continue
		
		#if it is fully covered, write site info to ALL pos outfile
		Oap.write("\t".join(l[0:5])+"\n")
		
		#check if it is ti or tv, turn flag on
		if "".join(l[3:5]) not in transitionalleles:
			ntflag=1
		#print("\t".join(l[0:5])+"\n")
		#print(A1cops)
		#print(A2cops)
		##Not necessary any more
		##for each pop of interest, get treemix majcounts/mincounts format
		##first, get number of total copies from ncop file
		#lnc=linenc.split()
		#ncops=[]
		#for i in p:
		#	ncops.append(lnc[i])
		#
		##then get number of "first allele" by doing freq*ncops in a vector
		##this SHOULD NOT fail because we are already sure it is 
		##a site where all populations have coverage
		#cnts=[]
		#for i in xrange(len(freqs)):
			#cnts.append(round(float(freqs[i])*float(ncops[i])))
		#
		##finally, get number of "second allele" by doing copies-first vectors
		#cnts2=[]
		#for i in xrange(len(freqs)):
		#	cnts2.append(round(float(ncops[i])-cnts[i]))
		#
		#format line to write
		tmcols=[]
		#for i in xrange(len(freqs)):
		for i in xrange(len(A1cops)):
			tmcols.append(str(int(A1cops[i]))+","+str(int(A2cops[i])))
		
		fmtline=" ".join(tmcols)+"\n"
		
		#write to ALL treemix outfile
		if pyver==2:
			Oa.stdin.write(fmtline)
		else:
			Oa.stdin.write((fmtline).encode())
		
		#if nt flag is on, write pos to NT pos outfile
		#if nt flag is on, write site to NT outfile
		if ntflag==1:
			Onp.write("\t".join(l[0:5])+"\n")
			if pyver==2:
				On.stdin.write(fmtline+"")
			else:
				On.stdin.write((fmtline+"").encode())

#close everything

Oap.close()
Onp.close()
Oa.communicate()
On.communicate()



