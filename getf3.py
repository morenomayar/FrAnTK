
from __future__ import division
import math
import gzip
import sys
import os
import subprocess

#deal with py2 vs py3
pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["h1", "h2", "target", "bsize", "freqpref", "resdir", "category", "singlehapout"]

USAGE="""
Compute f3 statistics of the form f3(h1, h2; target), using definition and normalization in Patterson et al, 2012 (Genetics).
Arguments should be passed as: arg=value.
h1\tString. Name of population to set in h1 (should be present in prefix_pop). Required. 
h2\tString. Name of population to set in h2 (should be present in prefix_pop). Required. 
target\tString. Name of population to set in targe (should be present in prefix_pop). Required. 
bsize\tInt. Length of the blocks (in bp) for block-jackknife procedure. This value should correspond approximately to the size of the linkage blocks in the reference genome. Default: 5000000.
freqpref\tString. Prefix from prefix_freqs.gz file created by BuildFreqs.py.
resdir\tString. Name of directory where results will be placed. Default: res.
category\tString. Used for automatic plotting. This will be appended to the result. Default: Other.
singlehapout\t(0 | 1). Set to 1 when target is a single pseudo-haploid outgroup. This will disable het. normalisation from Patterson 2012 and set f3 denominator to 1. Default: 0. 
Sample call: getf3.py h1=popnameh1 h2=popnameh2 target=popnametarget bsize=3500000 freqpref=prefix resdir=dirname category=somecategory
Run without arguments to get this message. 
"""


argvect=sys.argv[1:]
if len(argvect)==0:
	print(USAGE)
	sys.exit(1)

arghash={}
arghash["bsize"]=5000000
arghash["singlehapout"]=0
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
h3=arghash["target"]
bsize=arghash["bsize"]
freqpref=arghash["freqpref"]
resdir=arghash["resdir"]
Cat=arghash["category"]
singlehapout=arghash["singlehapout"]

try:
	bsize=int(bsize)
except:
	print("bsize should be an int\n")

try:
	singlehapout=int(singlehapout)
except:
	print("singlehapout should be 0 or 1\n")

if singlehapout!=1 and singlehapout!=0:
	print("singlehapout should be 0 or 1")
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

#get contig names, starts and ends from _chrs file
chrsfile=freqpref+"_chrs"
l=open(chrsfile, "r").readlines()
chr=[]
firstsnps=[]
lens=[]
for i in l:
	chr.append(i.split()[0])
	firstsnps.append(int(i.split()[1]))
	lens.append(int(i.split()[2]))

#define blocks. Each block has a chrname, a start and an end
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


#read population names from panel
l=open(freqpref+"_pop", "r").readlines()
nl=[]
for i in l:
	nl.append([i.split()[0]])

#find index of h1
h1p=-1
for i in xrange(len(nl)):
	if nl[i] == [h1]:
		h1p=i
		break

if h1p == -1:
	print("population 1 not found in pop names")
	print(USAGE)
	sys.exit(1)

#find index of h2
h2p=-1
for i in xrange(len(nl)):
	if nl[i] == [h2]:
		h2p=i
		break

if h2p == -1:
	print("population 2 not found in pop names")
	print(USAGE)
	sys.exit(1)

#find index of h3
h3p=-1
for i in xrange(len(nl)):
	if nl[i] == [h3]:
		h3p=i
		break

if h3p == -1:
	print("target population not found in pop names")
	print(USAGE)
	sys.exit(1)

#get base name of freqs file (no path, and no extensions)
freqprefbase=freqpref.split("/")[len(freqpref.split("/"))-1]

#btops and bbots will be lists of lists. Each element has all the tops and
#bots from each block
btops=[]
bbots=[]
#currb keeps track of which block we are filling
currb=0
currbtops=[]
currbbots=[]
#with gzip.open(freqpref+"_freqs.gz") as f:
with subprocess.Popen(["zcat",freqpref+"_freqs.gz"],stdout=subprocess.PIPE).stdout as f:
	for line in f:
		#read freqs line and parse chr, pos and freqs (using indices from above)
		if pyver==2:
			l=line.split()
		else:
			l=line.decode().split()
		chr=l[0]
		pos=int(l[1])
		h1f=l[5+h1p*3]
		h2f=l[5+h2p*3]
		h3f=l[5+h3p*3]
		try:
			if singlehapout == 0:
				#this is hhat=(min*maj)/(tot*(tot-1))
				h3h=(float(l[5+h3p*3+1])*float(l[5+h3p*3+2]))/((float(l[5+h3p*3+1])+float(l[5+h3p*3+2]))*((float(l[5+h3p*3+1])+float(l[5+h3p*3+2]))-1))
				#this is the extraterm=hhat/tot
				h3e=h3h/(float(l[5+h3p*3+1])+float(l[5+h3p*3+2]))
				top=((float(h3f)-float(h1f))*(float(h3f)-float(h2f)))-float(h3e)
				bot=2*float(h3h)
			else:
				tot=float(l[5+h3p*3+1])+float(l[5+h3p*3+2])
				if tot != 1:
					h3h=(float(l[5+h3p*3+1])*float(l[5+h3p*3+2]))/((float(l[5+h3p*3+1])+float(l[5+h3p*3+2]))*((float(l[5+h3p*3+1])+float(l[5+h3p*3+2]))-1))
				else:   
					h3h=(float(l[5+h3p*3+1])*float(l[5+h3p*3+2]))
				h3e=h3h/(float(l[5+h3p*3+1])+float(l[5+h3p*3+2]))
				top=((float(h3f)-float(h1f))*(float(h3f)-float(h2f)))-float(h3e)
				bot=1
		except:
			#the try block could fail because float casts did not work (missing data). If so, we skip the site
			continue
		if top == 0 and bot == 0:
			#if numerator and denom are 0, we skip the site
			continue
		if chr == Bchrs[currb] and pos >= Bstarts[currb] and pos <= Bends[currb]:
			#if the current site is within the bounds of the block we are 
			#standing on (currb), we append to its tops and bots
			currbtops.append(top)
			currbbots.append(bot)
		else:
			#if the current site is NOT within the bounds of the block
			#we are standing on (currb), we take the current block
			#(tops and bots) and keep them as an element (a list) of 
			#the btops and bbots lists. 
			#then we get ready for a new block, and stand on the next
			#block (currb+=1). Up to this point, we have not added the
			#current site to its block (because we are not standing on the
			#right block yet). 
			btops.append(currbtops)
			bbots.append(currbbots)
			currbtops=[]
			currbbots=[]
			currb+=1
			#we have to make sure that the increase in currb left us on
			#the right block for the current site, otherwise we have to
			#keep increasing currb until we are on the right block
			while not(chr == Bchrs[currb] and pos >= Bstarts[currb] and pos <= Bends[currb]):
				#if we are not standing on the right block, we append
				#empty lists in the btops and bbots lists. These 
				#correspond to blocks for which we have no data. Then
				#we increase currb to try the next block. We dont
				#go into this part once we are in the right block
				btops.append(currbtops)				
				bbots.append(currbbots)
				currbtops=[]
				currbbots=[]
				currb+=1
			#once we are on the right block (we got out of the while)
			#we add the current site to the current block
			currbtops.append(top)
			currbbots.append(bot)

#Above, we only append blocks to the btops and bbots lists when we change
#blocks. Since we dont change after the last block, we have to append the last 
#block manually. 
btops.append(currbtops)
bbots.append(currbbots)

#Find indices of empty blocks for removing them
ToRm=[]
for i in xrange(len(btops)):
	if len(btops[i]) == 0 or len(bbots[i]) == 0:
		ToRm.append(i)

#Remove empty blocks from btops and bbots. We do this in reverse
#order to preserve indices
for i in sorted(ToRm, reverse=True):
	del btops[i]
	del bbots[i]

#Jackknife
#For each block we are going sum tops and bots over all snps in the block
#props will have the number of snps per block for now
btopsums=[]
bbotsums=[]
Props=[]
for i in xrange(len(btops)):
	Props.append(len(btops[i]))
	btopsums.append(sum(btops[i]))
	bbotsums.append(sum(bbots[i]))

try:
	#point estimate is the sum of the sums of tops and bots
	f3=sum(btopsums)/sum(bbotsums)
	#get leave-1-block-out jackknife estimates for each block
	#recall we are only using blocks with data (empty ones removed above)
	jEsts=[]
	for i in xrange(len(btopsums)):
		jEsts.append(sum(btopsums[:i] + btopsums[(i + 1):])/sum(bbotsums[:i] + bbotsums[(i + 1):]))
	#get number of blocks and number of snps
	nBlocks=len(jEsts)
	nSNPs=sum(Props)
	#nprops is the number of snps per block divided by the total number of snps
	#this is how much each block weighs
	nProps=[]
	for i in Props:
		nProps.append(i/nSNPs)
	#compute jackknife SE. Took this one from the angsd dstat. 
	Sumin=0
	for i in xrange(len(nProps)):
		Sumin+= (1-nProps[i])*jEsts[i]
	Sumout=0
	for i in xrange(len(nProps)):
		Sumout+= 1/(1/nProps[i]-1) * (1/nProps[i]*f3-(1/nProps[i]-1)*jEsts[i] - nBlocks*f3+ Sumin)**2
	jackSE=math.sqrt(1/nBlocks * Sumout)
	Z=f3/jackSE
	Z
	#print results
	O=open(resdir+"/f3_"+"_".join([freqprefbase, h1, h2, h3])+".txt", 'w')
	#print "\t".join(map(str, [h1, h2, h3, f3, jackSE, Z, nSNPs, nBlocks]))
	O.write("\t".join(map(str, [h1, h2, h3, f3, jackSE, Z, nSNPs, nBlocks, "TRUE", Cat]))+"\n")
	O.close()
except:
	#If we couldnt compute the above, means we had no data so we print NA res
	O=open(resdir+"/f3_"+"_".join([freqprefbase, h1, h2, h3])+".txt", 'w')
	O.write("\t".join(map(str, [h1, h2, h3, "NA", "NA", "NA", "0", "0", "FALSE", Cat]))+"\n")
	O.close()
