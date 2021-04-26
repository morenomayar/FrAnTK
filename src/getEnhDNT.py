
from __future__ import division
import math
import gzip
import sys
import os
import subprocess

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["h1", "h2", "h3", "h4", "bsize", "freqpref", "resdir", "category"]

USAGE="""
Compute enhanced D statistics of the form D(h1, h2; h3, h4) on sites where one allele is fixed in h4, using (allele frequency) definition and normalization in Patterson et al, 2012 (Genetics). Works better with many sites and not too many individuals in h4. 
Arguments should be passed as: arg=value.
h1\tString. Name of population to set in h1 (should be present in prefix_pop). Required. 
h2\tString. Name of population to set in h2 (should be present in prefix_pop). Required. 
h3\tString. Name of population to set in h3 (should be present in prefix_pop). Required. 
h4\tString. Comma-separated (no spaces) list of populations to set in h4 (should be present in prefix_pop). Note that the number of sites will decrease with the number of individuals/pops in h4. Required. 
bsize\tInt. Length of the blocks (in bp) for block-jackknife procedure. This value should correspond approximately to the size of the linkage blocks in the reference genome. Default: 5000000.
freqpref\tString. Prefix from prefix_freqs.gz file created by BuildFreqs.py.
resdir\tString. Name of directory where results will be placed. Default: res.
category\tString. Used for automatic plotting. This will be appended to the result. Default: Other.
Sample call: getEnhD.py h1=popnameh1 h2=popnameh2 h3=popnameh3 h4=popnameh4_1,popnameh4_2,popnameh4_3 bsize=3500000 freqpref=prefix resdir=dirname category=somecategory
Run without arguments to get this message. 
"""


argvect=sys.argv[1:]
if len(argvect)==0:
	print(USAGE)
	sys.exit(1)

arghash={}
arghash["bsize"]=5000000
arghash["resdir"]="res"
arghash["category"]="Other"
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

h1=arghash["h1"]
h2=arghash["h2"]
h3=arghash["h3"]
h4=arghash["h4"]
bsize=arghash["bsize"]
freqpref=arghash["freqpref"]
resdir=arghash["resdir"]
Cat=arghash["category"]

try:
	bsize=int(bsize)
except:
	print("bsize should be an int\n")

if os.path.exists(resdir) and os.path.isdir(resdir):
	r=0
elif os.path.exists(resdir) and not os.path.isdir(resdir):
	print(resdir+" is a file, not overwriting.\n")
	print(USAGE)
	sys.exit(1)
elif not os.path.exists(resdir):
	r=os.system("mkdir "+resdir)

chrsfile=freqpref+"_chrs"
l=open(chrsfile, "r").readlines()
chr=[]
firstsnps=[]
lens=[]
for i in l:
	chr.append(i.split()[0])
	firstsnps.append(int(i.split()[1]))
	lens.append(int(i.split()[2]))

bchrs=[]
bstarts=[]
bends=[]
for i in xrange(len(chr)):
	#starts=range(1,lens[i], bsize)
	starts=range(firstsnps[i],lens[i], bsize)
	ends=[]
	for j in starts[1:]:
		ends.append(j-1)
	ends.append(lens[i])
	chrs=[]
	for j in starts:
		chrs.append(chr[i])
	bchrs.append(chrs)
	bstarts.append(starts)
	bends.append(ends)


Bchrs=[]
for i in bchrs:
	for j in i:
		Bchrs.append(j)

Bstarts=[]
for i in bstarts:
	for j in i:
		Bstarts.append(j)

Bends=[]
for i in bends:
	for j in i:
		Bends.append(j)

l=open(freqpref+"_pop", "r").readlines()
nl=[]
for i in l:
	nl.append([i.split()[0]])

h1p=-1
for i in xrange(len(nl)):
	if nl[i] == [h1]:
		h1p=i
		break

if h1p == -1:
	print("population 1 not found in pop names")
	print(USAGE)
	sys.exit(1)

h2p=-1
for i in xrange(len(nl)):
	if nl[i] == [h2]:
		h2p=i
		break

if h2p == -1:
	print("population 2 not found in pop names")
	print(USAGE)
	sys.exit(1)

h3p=-1
for i in xrange(len(nl)):
	if nl[i] == [h3]:
		h3p=i
		break

if h3p == -1:
	print("population 3 not found in pop names")
	print(USAGE)
	sys.exit(1)

h4name=h4
h4=h4.split(",")
noutgs=len(h4)
h4p=[-1]
for outg in h4:
	for i in xrange(len(nl)):
		if nl[i] == [outg]:
			h4p.append(i)
			break

h4p=h4p[1:]

if h4p == [-1] or len(h4p) != noutgs:
	print("population 4 not found in pop names")
	print(USAGE)
	sys.exit(1)

freqprefbase=freqpref.split("/")[len(freqpref.split("/"))-1]

btops=[]
bbots=[]
currb=0
currbtops=[]
currbbots=[]
abba=0
baba=0
invarsites=0
transitionalleles=["CT", "GA", "TC", "AG"]
#with gzip.open(freqpref+"_freqs.gz") as f:
with subprocess.Popen(["zcat",freqpref+"_freqs.gz"],stdout=subprocess.PIPE).stdout as f:
	for line in f:
		if pyver==2:
			l=line.split()
		else:
			l=line.decode().split()
		chr=l[0]
		pos=int(l[1])
		alleles=l[3]+l[4]
		if alleles in transitionalleles:
			continue
		h1f=l[5+h1p*3]
		h2f=l[5+h2p*3]
		h3f=l[5+h3p*3]
		h4fvec=[]
		for outg in h4p:
			h4fvec.append(l[5+outg*3])
		try:
			h4fvec=list(map(float, h4fvec))
		except:
			continue
		h4f=sum(h4fvec)/noutgs
		try:
			h1f=float(h1f)
			h2f=float(h2f)
			h3f=float(h3f)
			h4f=float(h4f)
			#if not ((h4f==0 and h3f==1) or (h4f==1 and h3f==0)):
			if not ((h4f==0 and h3f!=0) or (h4f==1 and h3f!=1)):
				continue
			invarsites+=1
			top=(h1f-h2f)*(h3f-h4f)
			bot=(h1f+h2f-2*h1f*h2f)*(h3f+h4f-2*h3f*h4f)
			abba+=((1-h1f)*h2f*h3f*(1-h4f))+(h1f*(1-h2f)*(1-h3f)*h4f)
			baba+=(h1f*(1-h2f)*h3f*(1-h4f))+((1-h1f)*h2f*(1-h3f)*h4f)
		except:
			continue
		if top == 0 and bot == 0:
			continue
		if chr == Bchrs[currb] and pos >= Bstarts[currb] and pos <= Bends[currb]:
			currbtops.append(top)
			currbbots.append(bot)
		else:
			btops.append(currbtops)
			bbots.append(currbbots)
			currbtops=[]
			currbbots=[]
			currb+=1
			while not(chr == Bchrs[currb] and pos >= Bstarts[currb] and pos <= Bends[currb]):
				btops.append(currbtops)				
				bbots.append(currbbots)
				currbtops=[]
				currbbots=[]
				currb+=1
			currbtops.append(top)
			currbbots.append(bot)

btops.append(currbtops)
bbots.append(currbbots)

ToRm=[]
for i in xrange(len(btops)):
	if len(btops[i]) == 0 or len(bbots[i]) == 0:
		ToRm.append(i)

for i in sorted(ToRm, reverse=True):
	del btops[i]
	del bbots[i]

btopsums=[]
bbotsums=[]
Props=[]
for i in xrange(len(btops)):
	Props.append(len(btops[i]))
	btopsums.append(sum(btops[i]))
	bbotsums.append(sum(bbots[i]))

try:
	Dstat=sum(btopsums)/sum(bbotsums)
	jEsts=[]
	for i in xrange(len(btopsums)):
		jEsts.append(sum(btopsums[:i] + btopsums[(i + 1):])/sum(bbotsums[:i] + bbotsums[(i + 1):]))
	nBlocks=len(jEsts)
	nSNPs=sum(Props)
	nProps=[]
	for i in Props:
		nProps.append(i/nSNPs)
	Sumin=0
	for i in xrange(len(nProps)):
		Sumin+= (1-nProps[i])*jEsts[i]
	Sumout=0
	for i in xrange(len(nProps)):
		Sumout+= 1/(1/nProps[i]-1) * (1/nProps[i]*Dstat-(1/nProps[i]-1)*jEsts[i] - nBlocks*Dstat+ Sumin)**2
	jackSE=math.sqrt(1/nBlocks * Sumout)
	Z=Dstat/jackSE
	Z
	O=open(resdir+"/D_"+"_".join([freqprefbase, h1, h2, h3, h4name])+".txt", 'w')
	#print("\t".join(map(str, [h1, h2, h3, h4name, Dstat, jackSE, Z, nSNPs, nBlocks, "TRUE", Cat]))+"\n")
	O.write("\t".join(map(str, [h1, h2, h3, h4name, Dstat, jackSE, Z, nSNPs, nBlocks, "TRUE", Cat, abba, baba]))+"\n")
	O.close()
except:
	O=open(resdir+"/D_"+"_".join([freqprefbase, h1, h2, h3, h4name])+".txt", 'w')
	#print("\t".join(map(str, [h1, h2, h3, h4name, "NA", "NA", "NA", 0, 0, "FALSE", Cat]))+"\n")
	O.write("\t".join(map(str, [h1, h2, h3, h4name, "NA", "NA", "NA", 0, 0, "FALSE", Cat, abba, baba]))+"\n")
	O.close()

