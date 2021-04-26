
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

allowedpars=["h1", "h2", "h3", "h4", "anc", "minder", "maxder", "mincops", "bsize", "freqpref", "resdir", "category"]

USAGE="""
Compute D statistics of the form D(h1, h2; h3, h4), using (allele frequency) definition and normalization in Patterson et al, 2012 (Genetics). We only consider sites where there are (minder, maxder] derived alleles in the complete panel. 
Arguments should be passed as: arg=value.
h1\tString. Name of population to set in h1 (should be present in prefix_pop). Required. 
h2\tString. Name of population to set in h2. This is the pop where we ascertain. (should be present in prefix_pop). Required. 
h3\tString. Name of population to set in h3 (should be present in prefix_pop). Required. 
h4\tString. Name of population to set in h4 (should be present in prefix_pop). Required. 
anc\tString. Name of population with the ancestral allele (should be a single hap ind and should be present in prefix_pop). Required. 
minder\tInt. Minimum number of derived alleles in h1+h2+h3+h4 together. This is the left endpoint of an interval and should be open, e.g., set to 1 for at least 2 derived alleles in h1+h2+h3+h4 together. Required. 
maxder\tInt. Maximum number of derived alleles in h1+h2+h3+h4 together. This is the right endpoint of an interval and should be closed, e.g., set to 3 for at most 3 derived alleles in h1+h2+h3+h4 together. Required. 
mincops\tString. Comma-separated (no spaces) list of minimum number of total copies in h1,h2,h3,h4, e.g., 1,1,1,1 for no filter. Required. 
bsize\tInt. Length of the blocks (in bp) for block-jackknife procedure. This value should correspond approximately to the size of the linkage blocks in the reference genome. Default: 5000000.
freqpref\tString. Prefix from prefix_freqs.gz file created by BuildFreqs.py.
resdir\tString. Name of directory where results will be placed. Default: res.
category\tString. Used for automatic plotting. This will be appended to the result. Default: Other. 
IMPORTANT: Make sure that pseudo-haploid populations were specified in the _clust file before running BuildFreqs.py. Otherwise, derived alleles in pseudo-haploid pops will be counted twice. 
Sample call: getDstrat.py h1=popnameh1 h2=popnameh2 h3=popnameh3 h4=popnameh4 anc=popnameanc bsize=3500000 minder=2 maxder=3 mincops=20,20,20,20 freqpref=prefix resdir=dirname category=somecategory
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
anc=arghash["anc"]
minder=arghash["minder"]
maxder=arghash["maxder"]
mincops=arghash["mincops"]
bsize=arghash["bsize"]
freqpref=arghash["freqpref"]
resdir=arghash["resdir"]
Cat=arghash["category"]

try:
	bsize=int(bsize)
except:
	print("bsize should be an int\n")

try:
	minder=int(minder)
	maxder=int(maxder)
except:
	print("minder and maxder should be an int\n")

if maxder <= minder:
	print("maxder should be greater than minder")
	print(USAGE)
	sys.exit(1)

if minder < 1 or maxder < 2:
	print("maxder and minder should be >=2 and >=1, respectively")
	print(USAGE)
	sys.exit(1)

mincopsvec=mincops.split(",")
mincopslen=len(mincopsvec)
if mincopslen != 4:
	print("mincops should be a comma-separated (no spaces) list of four ints, e.g., 12,13,14,13")
	print(USAGE)
	sys.exit(1)

for i in xrange(mincopslen):
	try:
		mincopsvec[i]=int(mincopsvec[i])
	except:
		print("elements of mincops should be int\n")

for i in xrange(mincopslen):
	if mincopsvec[i] <= 0:
		print("elements of mincops should be positive")
		print(USAGE)
		sys.exit(1)

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

h4p=-1
for i in xrange(len(nl)):
	if nl[i] == [h4]:
		h4p=i
		break

if h4p == -1:
	print("population 4 not found in pop names")
	print(USAGE)
	sys.exit(1)

ancp=-1
for i in xrange(len(nl)):
	if nl[i] == [anc]:
		ancp=i
		break

if ancp == -1:
	print("population anc not found in pop names")
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
#with gzip.open(freqpref+"_freqs.gz") as f:
with subprocess.Popen(["zcat",freqpref+"_freqs.gz"],stdout=subprocess.PIPE).stdout as f:
	for line in f:
		if pyver==2:
			l=line.split()
		else:
			l=line.decode().split()
		chr=l[0]
		pos=int(l[1])
		h1f=l[5+h1p*3]
		h2f=l[5+h2p*3]
		h3f=l[5+h3p*3]
		h4f=l[5+h4p*3]
		ancf=l[5+ancp*3]
		try:
			h1f=float(h1f)
			h2f=float(h2f)
			h3f=float(h3f)
			h4f=float(h4f)
			ancf=float(ancf)
			top=(h1f-h2f)*(h3f-h4f)
			bot=(h1f+h2f-2*h1f*h2f)*(h3f+h4f-2*h3f*h4f)
		except:
			continue
		if top == 0 and bot == 0:
			continue
		if ancf == 1:
			#arrange like this because anc has a1, so der=a2:
			#SECOND of the last two elements of l
			#this is the total number of copies of a2
			nofder=int(l[-1])
			#FIRST of the last two elements of l
			#total number of copies of a1
			nofanc=int(l[-2])
			nofcop=[]
			for i in (h1p, h2p, h3p, h4p):
				nofcop.append(int(l[2+5+i*3])+int(l[1+5+i*3]))
		elif ancf == 0:
			#arrange like this because anc has a2, so der=a1:
			#FIRST of the last two elements of l
			#this is the total number of copies of a1
			nofder=int(l[-2])
			#SECOND of the last two elements of l
			#total number of copies of a2
			nofanc=int(l[-1])
			nofcop=[]
			for i in (h1p, h2p, h3p, h4p):
				nofcop.append(int(l[1+5+i*3])+int(l[2+5+i*3]))
		else:
			print("freqs in anc should be 0 or 1")
			print(str(chr))
			print(str(pos))
			sys.exit(1)
		ncopsok=1
		for i in xrange(mincopslen):
			if nofcop[i] < mincopsvec[i]:
				ncopsok=0
		if nofder <= minder or nofder > maxder or ncopsok == 0:
			continue
		abba+=((1-h1f)*h2f*h3f*(1-h4f))+(h1f*(1-h2f)*(1-h3f)*h4f)
		baba+=(h1f*(1-h2f)*h3f*(1-h4f))+((1-h1f)*h2f*(1-h3f)*h4f)
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
	O=open(resdir+"/Dstrat2_"+"_".join([freqprefbase, h1, h2, h3, h4, anc, str(minder), str(maxder), str(mincops)])+".txt", 'w')
	#print("\t".join(map(str, [h1, h2, h3, h4, Dstat, jackSE, Z, nSNPs, nBlocks])))
	O.write("\t".join(map(str, [h1, h2, h3, h4, Dstat, jackSE, Z, nSNPs, nBlocks, "TRUE", Cat, abba, baba]))+"\n")
	O.close()
except:
	O=open(resdir+"/Dstrat2_"+"_".join([freqprefbase, h1, h2, h3, h4, anc, str(minder), str(maxder), str(mincops)])+".txt", 'w')
	#print("\t".join(map(str, [h1, h2, h3, h4, Dstat, jackSE, Z, nSNPs, nBlocks])))
	O.write("\t".join(map(str, [h1, h2, h3, h4, "NA", "NA", "NA", 0, 0, "FALSE", Cat, abba, baba]))+"\n")
	O.close()

