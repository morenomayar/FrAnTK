
from __future__ import division
import gzip
import sys
import os
import subprocess

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["plinkpref", "prefout"]

USAGE="""
Get dummy allele frequencies from a bim file. Helpful when we only know the sites
where we want to produce pseudo-haploids. 
Arguments should be passed as: arg=value.
plinkpref\tString. Prefix of a bim plink file. Required. 
prefout\tString. Prefix for output file names. Should be different from plinkpref. Required. 
NOTE: You should have permissions to write in the directory where the plink file is located. 
Sample call: BuildDummyFreqs.py plinkpref=plinkfilename prefout=prefix
Run without arguments to get this message. 
"""

argvect=sys.argv[1:]
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

plinkpref=arghash["plinkpref"]
prefout=arghash["prefout"]

if plinkpref == prefout:
	print("plinkpref and prefout should be different")
	print(USAGE)
	sys.exit(1)

#write empty pop file

print("Parsing pop info... ")
Oh=open(prefout+"_pop", 'w')
Oh.close()

#Ob=open(plinkpref+".bim", "r")

print("Parsing bim... ")
O=subprocess.Popen("gzip -c > "+prefout+"_freqs.gz", shell=True, stdin=subprocess.PIPE)
Or=open(prefout+"_regions", "w")

cont=0
site=[]
nsites=0
chrs4chrs=[]
pos4chrs=[]
starts4chrs=[]
idxlines=[]
lastchr=0
lastpos=10000000000
#cops=[]
A1copsums=0
A2copsums=0
currminpos=1

with open(plinkpref+".bim", "r") as Ob:
	for line in Ob:
		bimline=line.split()
		nsites+=1
		#print nsites
		#write
		if pyver==2:
			O.stdin.write(bimline[0]+"\t"+bimline[3]+"\t"+bimline[1]+"\t"+bimline[4]+"\t"+bimline[5]+"\t0\t0\n")
		else:
			O.stdin.write((bimline[0]+"\t"+bimline[3]+"\t"+bimline[1]+"\t"+bimline[4]+"\t"+bimline[5]+"\t0\t0\n").encode())
		try:
			regchr=str(int(bimline[0]))
		except:
			regchr=str(bimline[0])
		regposs=str(int(bimline[3])-1)
		regpose=str(int(bimline[3]))
		regmaj=str(bimline[4])
		regmin=str(bimline[5])
		Or.write(regchr+"\t"+regposs+"\t"+regpose+"\t"+regmaj+"\t"+regmin+"\n")
		regposen=int(regpose)
		if regchr != lastchr:
			chrs4chrs.append(lastchr)
			pos4chrs.append(lastpos+1)
			idxlines.append(nsites)
			currminpos=regposen
			starts4chrs.append(currminpos)
			#print(str(lastchr)+"\t"+str(currminpos)+"\t"+str(lastpos+1)+"\t"+str(idxlines[-1])+"\n")
			lastchr=regchr
		lastpos=regposen

try:
	regchr=str(int(bimline[0]))
except:
	regchr=str(bimline[0])

regposs=str(int(bimline[3])-1)
regpose=str(int(bimline[3]))
regmaj=str(bimline[4])
regmin=str(bimline[5])
regposen=int(regpose)
lastpos=regposen
#Or.write(regchr+"\t"+regposs+"\t"+regpose+"\t"+regmaj+"\t"+regmin+"\n")
chrs4chrs.append(lastchr)
pos4chrs.append(lastpos+1)

#O.close()
O.communicate()
#Oc.communicate()
Or.close()

print("Parsing chr info... ")
#print chrs file
del chrs4chrs[0]
del pos4chrs[0]

idxlinesend=[]
for i in idxlines[1:]:
	idxlinesend.append(i-1)
#IMPORTANT DIFFERENCE TO NORMAL BuildFreqs.py: Here, we do nsites
#instead of nsites+1 since we are looping with the bimfile, and
#not with the freqs (and reading bim lines each time)
idxlinesend.append(nsites)

cchrs=map(str, chrs4chrs)
cfirstpos=map(str, starts4chrs)
cpos=map(str, pos4chrs)
cstarts=map(str, idxlines)
cends=map(str, idxlinesend)

Oc=open(prefout+"_chrs", 'w')
for chr, firstpos, pos, istart, iend in zip(cchrs, cfirstpos, cpos, cstarts, cends):
	Oc.write(chr+"\t"+firstpos+"\t"+pos+"\t"+istart+"\t"+iend+"\n")

Oc.close()

