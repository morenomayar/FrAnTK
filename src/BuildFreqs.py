
from __future__ import division
import gzip
import sys
import os
import subprocess

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["plinkpref", "clustfile", "npops", "prefout"]

USAGE="""
Get allele frequencies from a plink file for f-statistics computation. 
Arguments should be passed as: arg=value.
plinkpref\tString. Prefix of a bed-bim-fam plink file. Required. 
clustfile\tString. Tab-separated list of individuals and populations to which they belong, passed to --within in plink. This file may also include information about the ploidy of each individual, in a potential fourth column. If the fourth column is set to 1, then we assume that individual only contributes one copy of the allele (it is haploid). Any other value will be set to the default (diploid). Note that for rows without a fourth column, the individual will be set to diploid. IMPORTANT: The ploidy of all individuals in a population should be the same. If different values are specified, the ploidy of the pop is set to the most frequent value (diploid for ties). Format: FamName IndName PopName 0|1. Required. 
npops\tInt. Total number of populations in the plink file. Required. 
prefout\tString. Prefix for output file names. Should be different from plinkpref. Required. 
NOTE: You should have permissions to write in the directory where the plink file is located. 
IMPORTANT: Turning real diploids to haploids is quite dangerous as that may lead to 0 counts. For example, if a single individual pop is het at a site (1,1), that would give rise to (0,0) when dividing by 2. In general there will be issues when dividing odd numbers, since it is an int division, int(3/2)=1. This is not a problem with real haploids as numbers will be even all the time. 
Sample call: BuildFreqs.py plinkpref=plinkfilename clustfile=popfilename npops=120 prefout=prefix
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
clustfile=arghash["clustfile"]
npops=arghash["npops"]
prefout=arghash["prefout"]

if plinkpref == prefout:
	print("plinkpref and prefout should be different")
	print(USAGE)
	sys.exit(1)

try:
	npops=int(npops)
except:
	print("npops should be an int")

print("Running plink cmd... ")
r=os.system("plink --allow-extra-chr --silent --allow-no-sex --bfile "+ plinkpref+" --freq --within "+clustfile+" --out "+plinkpref)

#initialize popdict with ancient/modern info
popdict={}

print("Parsing pop info... ")
Oh=open(prefout+"_pop", 'w')
Ob=open(plinkpref+".bim", "r")
with open(plinkpref+".frq.strat") as f:
	noline=f.readline()
	for i in xrange(npops):
		line=f.readline().split() 
		Oh.write(line[2]+"\n")
		popdict[line[2]]=[]

Oh.close()

#parse clust file and fill the above dictionary with keys:populations 
#and values:ancient 0|1. Since all individuals in a population SHOULD be 
#the same, the most frequent value will be chosen. In case of a tie, the 
#default will be 0 (diploid). If no value is found, the default will 
#be 0 (diploid). If something is NOT 1, it will be turned to 0. 

print("Parsing hap/dip info... ")
with open(clustfile) as f:
	for line in f:
		l=line.split()
		if len(l)==4:
			if l[3]=="1":
				popdict[l[2]].append(1)
			elif l[3]=="0":
				popdict[l[2]].append(0)
			else:
				popdict[l[2]].append(0)
				print(" ".join(l)+" set to diploid by default")
		else:
			print(" ".join(l)+" set to diploid by default (no 4th column)")
			popdict[l[2]].append(0)

#print len(popdict)
for i in popdict:
	#print i
	#print popdict[i]
	vals=list(set(popdict[i]))
	#print vals
	if len(vals)==1:
		popdict[i]=vals[0]
	else:
		#this takes advantage of the fact that the dict can
		#only be filled with 0 or 1 (see above)
		cnts=[0,0]
		for j in popdict[i]:
			cnts[j]+=1
		#print cnts
		if cnts[0]>=cnts[1]:
			print("Setting "+i+" to diploid")
			print(popdict[i])
			popdict[i]=0
		else:
			print("Setting "+i+" to haploid")
			print(popdict[i])
			popdict[i]=1
	#print popdict[i]

print("Writing pop info again (with hap/dip info)... ")
Oh=open(prefout+"_pop", 'w')
with open(plinkpref+".frq.strat") as f:
	noline=f.readline()
	for i in xrange(npops):
		line=f.readline().split() 
		Oh.write(line[2]+"\t"+str(popdict[line[2]])+"\n")

Oh.close()


print("Parsing strat file (this might take a while)... ")
#O=gzip.open(prefout+"_freqs.gz", 'w')
O=subprocess.Popen("gzip -c > "+prefout+"_freqs.gz", shell=True, stdin=subprocess.PIPE)
#Oc=subprocess.Popen("gzip -c > "+prefout+"_ncop.gz", shell=True, stdin=subprocess.PIPE)
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
with open(plinkpref+".frq.strat") as f:
	noline=f.readline()
	for line in f:
		cont+=1
		if cont<npops+1:
			sl=line.split()
			#the columns in this new version of the
			#freqs file will look like this:
			#frequency\tncopsA1\tncopsA2
			#it is very important to note that the
			#freq output by the plink freq cmd is the freq
			#of the allele in column A1 (which is the minor allele)
			#MAC is the count that corresponds to the freqs and A1
			maf=float(sl[5])
			A1cops=int(sl[6])
			A2cops=int(sl[7])-A1cops
			#print str(A1cops)+"\t"+str(A2cops)+"\n"
			#check if this population is non-missing for this site
			if int(sl[7])!=0:
				#check if population is haploid, 
				#if it is sampled then divide cops by two
				#population in strat file SHOULD be in popdict
				#print sl[2]
				#print popdict[sl[2]]
				if popdict[sl[2]]==1:
					#print "entre para"+sl[2]+"\n"
					A1cops=int(A1cops/2)
					A2cops=int(A2cops/2)
				#create line
				l=[maf, A1cops, A2cops]
				#print "A1copsums was "+str(A1copsums)+" "+" adding "+str(A1cops)+"\n"
				A1copsums+=A1cops
				#print "A1copsums is now "+str(A1copsums)+"\n"
				A2copsums+=A2cops
			else:
				l=["N", "N", "N"]
				#ncop="0"
			site.append(l)
			#cops.append(ncop)
		else:
			bimline=Ob.readline().split()
			nsites+=1
			#print nsites
			#write
			fmtsite=[]
			for i in site:
				fmtsite.append(map(str, i))
			#O.write(bimline[0]+"\t"+bimline[3]+"\t"+sl[1]+"\t"+sl[3]+"\t"+sl[4]+"\t")
			#O.write("\t".join(map("\t".join, fmtsite))+"\n")
			if pyver==2:
				O.stdin.write(bimline[0]+"\t"+bimline[3]+"\t"+sl[1]+"\t"+sl[3]+"\t"+sl[4]+"\t")
				O.stdin.write("\t".join(map("\t".join, fmtsite))+"\t")
				#print "jjjjjjjjjjjjjjjjjjjjjjjjj\n"
				fmtcops=str(int(A1copsums))+"\t"+str(int(A2copsums))
				O.stdin.write(fmtcops+"\n")
			else:
				O.stdin.write((bimline[0]+"\t"+bimline[3]+"\t"+sl[1]+"\t"+sl[3]+"\t"+sl[4]+"\t").encode())
				O.stdin.write(("\t".join(map("\t".join, fmtsite))+"\t").encode())
				#print "jjjjjjjjjjjjjjjjjjjjjjjjj\n"
				fmtcops=str(int(A1copsums))+"\t"+str(int(A2copsums))
				O.stdin.write((fmtcops+"\n").encode())
			try:
				regchr=str(int(bimline[0]))
			except:
				regchr=str(bimline[0])
			regposs=str(int(bimline[3])-1)
			regpose=str(int(bimline[3]))
			regmaj=str(sl[3])
			regmin=str(sl[4])
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
			cont=1
			site=[]
			A1copsums=0
			A2copsums=0
			#cops=[]
			sl=line.split()
			maf=float(sl[5])
			A1cops=int(sl[6])
			A2cops=int(sl[7])-A1cops
			#print str(A1cops)+"\t"+str(A2cops)+"\n"
			#check if this population is non-missing for this site
			if int(sl[7])!=0:
				#check if population is haploid, 
				#if it is sampled then divide cops by two
				#population in strat file SHOULD be in popdict
				#print sl[2]
				#print popdict[sl[2]]
				if popdict[sl[2]]==1:
					#print "entre para"+sl[2]+"\n"
					A1cops=int(A1cops/2)
					A2cops=int(A2cops/2)
				#create line
				l=[maf, A1cops, A2cops]
				#print "A1copsums was "+str(A1copsums)+" "+" adding "+str(A1cops)+"\n"
				A1copsums+=A1cops
				#print "A1copsums is now "+str(A1copsums)+"\n"
				A2copsums+=A2cops

			else:
				l=["N", "N", "N"]
				#ncop="0"
			site.append(l)
			#cops.append(ncop)

bimline=Ob.readline().split()
fmtsite=[]
for i in site:
	fmtsite.append(map(str, i))

#O.write(bimline[0]+"\t"+bimline[3]+"\t"+sl[1]+"\t"+sl[3]+"\t"+sl[4]+"\t")
#O.write("\t".join(map("\t".join, fmtsite))+"\n")
if pyver==2:
	O.stdin.write(bimline[0]+"\t"+bimline[3]+"\t"+sl[1]+"\t"+sl[3]+"\t"+sl[4]+"\t")
	O.stdin.write("\t".join(map("\t".join, fmtsite))+"\t")
	fmtcops=str(int(A1copsums))+"\t"+str(int(A2copsums))
	O.stdin.write(fmtcops+"\n")
else:
	O.stdin.write((bimline[0]+"\t"+bimline[3]+"\t"+sl[1]+"\t"+sl[3]+"\t"+sl[4]+"\t").encode())
	O.stdin.write(("\t".join(map("\t".join, fmtsite))+"\t").encode())
	fmtcops=str(int(A1copsums))+"\t"+str(int(A2copsums))
	O.stdin.write((fmtcops+"\n").encode())
try:
	regchr=str(int(bimline[0]))
except:
	regchr=str(bimline[0])
regposs=str(int(bimline[3])-1)
regpose=str(int(bimline[3]))
regmaj=str(sl[3])
regmin=str(sl[4])
regposen=int(regpose)
lastpos=regposen
Or.write(regchr+"\t"+regposs+"\t"+regpose+"\t"+regmaj+"\t"+regmin+"\n")
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
idxlinesend.append(nsites+1)

cchrs=map(str, chrs4chrs)
cfirstpos=map(str, starts4chrs)
cpos=map(str, pos4chrs)
cstarts=map(str, idxlines)
cends=map(str, idxlinesend)

Oc=open(prefout+"_chrs", 'w')
for chr, firstpos, pos, istart, iend in zip(cchrs, cfirstpos, cpos, cstarts, cends):
    Oc.write(chr+"\t"+firstpos+"\t"+pos+"\t"+istart+"\t"+iend+"\n")

Oc.close()

#cleanup
cmd="rm "+plinkpref+".frq.strat"
r=os.system(cmd)


















