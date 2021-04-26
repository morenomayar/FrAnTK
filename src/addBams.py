
from __future__ import division
import sys
import os
import subprocess
import string
import random

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["listname", "freqpref", "newpref", "nthr", "regionsfile"]

USAGE="""
Add sequencing data in one or more bam files to a precomputed frequency file (output from BuildFreqs.py). For each bam file, one random allele is sampled at each site present in the frequency file.
Arguments should be passed as: arg=value.
listname\tString. Required. Tab-separated file where each row contains:
\t\t1. Path to bam file. 
\t\t2. Individual name. 
\t\t3. Population name. 
\t\t4. Minimum mapping quality
\t\t5. Minimum base quality
\t\t6. Number of bases to trim from both ends of the reads.
freqpref\tString. Prefix from prefix_freqs.gz file created by BuildFreqs.py.
newpref\tString. Prefix for output file names. Required.
nthr\tInt. Number of threads to use. Each bam file is sent to one cpu. 
#regionsfile\tString. A _regions (.bed) file, similar to the one created by BuildFreqs, except contigs may have different names. This file should have the same number of rows as the original, as well as the same alleles (also same order). The first column of the file corresponds to contig names in bam. This argument is used if the bam and plink files have different contig/chromosome names. IMPORTANT: Different contig names in freqs and bam is BETA, so use at your own risk (same for bam2plink.py). Optional.
Sample call: addBams.py listname=list.txt freqpref= prefix newpref=prefixforout nthr=2
Run without arguments to get this message. 
"""

argvect=sys.argv[1:]
arghash={}
arghash["nthr"]=1
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

try:
	regionsfile=arghash["regionsfile"]
except:
	arghash["regionsfile"]=None
	regionsfile=None

for i in allowedpars:
	if i not in arghash:
		print(i+" is missing\n")
		print(USAGE)
		sys.exit(1)

listname=arghash["listname"]
freqpref=arghash["freqpref"]
newpref=arghash["newpref"]
nthr=arghash["nthr"]

if freqpref == newpref:
	print("freqpfref and newpref should be different")
	print(USAGE)
	sys.exit(1)

try:
	nthr=int(nthr)
except:
	print("nthr should be an int")
	print(USAGE)

def all2num(x):
	nall=[]
	for i in x:
		try:
			nall.append(int(i))
		except:
			nall.append(None)
	return nall

def numall(x):
	nall=0
	for i in x:
		if i != None:
			nall+=1
	return nall

def mincop(x):
	nmin=0
	for i in x:
		if i == 1:
			nmin+=1
	return nmin

freqprefbase=freqpref.split("/")[len(freqpref.split("/"))-1]

bamlist=open(listname, "r").readlines()

RandExt=''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in range(6))
O=open(freqprefbase+"_"+newpref+"_"+RandExt+"_samplecmds", "w")

nline=0
filenames={}
#note that filenames is a hash whose keys corresponds to population names
#and values correspond to arrays that contain the files in each population
for line in bamlist:
	nline+=1
	try:
		l=line.split()
		#O.write("bam2freqs.py bamfile="+l[0]+" freqpref="+freqpref+" trim="+l[5]+" MinMQ="+l[3]+" MinBQ="+l[4]+" indname="+l[1]+" popname="+l[2]+"\n")
		if regionsfile is None:
			O.write("bamfile="+l[0]+" freqpref="+freqpref+" trim="+l[5]+" MinMQ="+l[3]+" MinBQ="+l[4]+" indname="+l[1]+" popname="+l[2]+"\n")
		else:
			O.write("bamfile="+l[0]+" freqpref="+freqpref+" trim="+l[5]+" MinMQ="+l[3]+" MinBQ="+l[4]+" indname="+l[1]+" popname="+l[2]+" regionsfile="+regionsfile +"\n")
		try:
			filenames[l[2]].append(l[1]+"_"+l[2]+"_"+freqprefbase+"_spu.gz")
		except:
			filenames[l[2]]=[]
			filenames[l[2]].append(l[1]+"_"+l[2]+"_"+freqprefbase+"_spu.gz")
	except:
		print("not doing a command for line "+str(nline)+"\n")

O.close()

#xargcmd="cat "+freqpref+"_"+newpref+"_"+RandExt+"_samplecmds"+" | xargs -n 8 -P "+str(nthr)+" python"
if regionsfile is None:
	xargcmd="cat "+freqprefbase+"_"+newpref+"_"+RandExt+"_samplecmds"+" | xargs -n 7 -P "+str(nthr)+" bam2freqs.py"
else:
	xargcmd="cat "+freqprefbase+"_"+newpref+"_"+RandExt+"_samplecmds"+" | xargs -n 8 -P "+str(nthr)+" bam2freqs.py"

r=os.system(xargcmd)

freqfilenames=[]
#freqfilenamesnc=[]
poporder=[]
for i in filenames:
	#these two are the names of the freqs and ncop files, which will be
	#used later for adding them in order to the original freqs and ncop files
	freqfilenames.append(i+"_"+newpref+".gz")
	#freqfilenamesnc.append(i+"_"+newpref+"_ncop.gz")
	poporder.append(i)
	O=subprocess.Popen("gzip -c > "+i+"_"+newpref+".gz", shell=True, stdin=subprocess.PIPE)
	#On=subprocess.Popen("gzip -c > "+i+"_"+newpref+"_ncop.gz", shell=True, stdin=subprocess.PIPE)
	handles=[]
	for j in filenames[i]:
		handles.append(subprocess.Popen(["zcat",j],stdout=subprocess.PIPE).stdout)
	
	if len(handles)>1:
		for line in handles[0]:
			alleles=[]
			#oriall=[]
			if pyver==2:
				newall=line.split()[0]
			else:
				newall=line.decode().split()[0]
			alleles.append(newall)
			#oriall.append(newall)
			#print oriall
			for h in handles[1:]:
				if pyver==2:
					newall=h.readline().split()[0]
				else:
					newall=h.readline().decode().split()[0]
				alleles.append(newall)
				#oriall.append(newall)
				#print oriall
			
			alleles=all2num(alleles)
			nalleles=numall(alleles)
			#print str(nalleles)+"\n"
			if nalleles>1:
				#print "sisi\n"
				nminall=mincop(alleles)
				frqline=[str(nminall/nalleles)]
				#hhat=(((2*nminall)*((2*nalleles)-(2*nminall)))/((2*nalleles)*((2*nalleles)-1)))
				#eterm=hhat/(2*nalleles)
				frqline.append(str(int(nminall)))
				frqline.append(str(int(nalleles-nminall)))
				frqline="\t".join(frqline)+"\n"
			elif nalleles==1:
				nminall=mincop(alleles)
				frqline=[str(nminall/nalleles)]
				#hhat="None"
				#eterm="None"
				frqline.append(str(int(nminall)))
				frqline.append(str(int(nalleles-nminall)))
				frqline="\t".join(frqline)+"\n"
			else:
				frqline="N\tN\tN\n"
			if pyver==2:
				O.stdin.write(frqline)
			else:
				O.stdin.write((frqline).encode())
			#On.stdin.write(str(nalleles*2)+"\n")
			#print "\t".join(oriall)+"\t"+frqline
			#print "\n"
	else:
		#print "entre con "+i+"\n"
		for line in handles[0]:
			if pyver==2:
				newall=str(line.split()[0])
			else:
				newall=str(line.decode().split()[0])
			if newall=="0":
				frqline=newall+"\t0\t1\n"
			elif newall=="1":
				frqline=newall+"\t1\t0\n"
			else:
				frqline="N\tN\tN\n"
			if pyver==2:
				O.stdin.write(frqline)
			else:
				O.stdin.write((frqline).encode())
			#if newall != "None":
			#	On.stdin.write("2\n")
			#else:
			#	On.stdin.write("0\n")
	O.communicate()
	#On.communicate()

#create new freq file with extra cols
#open a handle for writing new freqs file
O=subprocess.Popen("gzip -c > "+newpref+"_freqs.gz", shell=True, stdin=subprocess.PIPE)
##Not necessary any more
##open a handle for writing new ncop file
#Onc=subprocess.Popen("gzip -c > "+newpref+"_ncop.gz", shell=True, stdin=subprocess.PIPE)

#create handles for reading pop freq files
handles=[]
for j in freqfilenames:
	handles.append(subprocess.Popen(["zcat",j],stdout=subprocess.PIPE).stdout)

##create handles for reading pop ncop files
#handlesnc=[]
#for j in freqfilenamesnc:
#	handlesnc.append(subprocess.Popen(["zcat",j],stdout=subprocess.PIPE).stdout)

##not necessary any more
##open a handle for reading original ncop file
#fnc=subprocess.Popen(["zcat",freqpref+"_ncop.gz"],stdout=subprocess.PIPE).stdout
#read through the original freqs file
with subprocess.Popen(["zcat",freqpref+"_freqs.gz"],stdout=subprocess.PIPE).stdout as f:
	for line in f:
		#split line from freqs file
		if pyver==2:
			fline=line.rstrip().split()
		else:
			fline=line.decode().rstrip().split()
		#get two last elements into an array
		ncops=[int(fline[len(fline)-2]), int(fline[len(fline)-1])]
		#rejoin line without the last two elements
		fline="\t".join(fline[:-2])
		##not necessary any more
		##read one line from the original ncop file and split it into an
		##array so that we can sum the counts for every population
		#ncops=map(int, fnc.readline().rstrip().split())
		
		#for each of the lines in the original freqs file, read one line for
		#each of the pop frequency files (kept in handles) and append to original
		#and keep the counts of minor and major alleles
		newcols=[]
		for j in handles:
			if pyver==2:
				popline=j.readline().rstrip()
			else:
				popline=j.readline().decode().rstrip()
			newcols.append(popline)
			pl=popline.split()
			if pl[1]!="N":
				ncops[0]+=int(pl[1])
				ncops[1]+=int(pl[2])
		
		#O.stdin.write("\t".join([line.rstrip(), "\t".join(newcols)])+"\n")
		if pyver==2:
			O.stdin.write("\t".join([fline, "\t".join(newcols)])+"\t")
		else:
			O.stdin.write(("\t".join([fline, "\t".join(newcols)])+"\t").encode())
		#write new counts
		#Onc.stdin.write(str(ncops[0])+"\t"+str(ncops[1])+"\n")
		if pyver==2:
			O.stdin.write(str(ncops[0])+"\t"+str(ncops[1])+"\n")
		else:
			O.stdin.write((str(ncops[0])+"\t"+str(ncops[1])+"\n").encode())

O.communicate()
#Onc.communicate()
			

#create new pop file with extra rows

l=open(freqpref+"_pop", "r").readlines()
nl=[]
ploidy=[]
for i in l:
	nl.append(i.split()[0])
	ploidy.append(i.split()[1])

for i in poporder:
	ploidy.append("1")
	if i not in nl:
		nl.append(i)
	else:
		nl.append(i+"_WGS")

Oh=open(newpref+"_pop", 'w')
for i in range(len(nl)):
	Oh.write(nl[i]+"\t"+ploidy[i]+"\n")

Oh.close()

#copy regions and chrs files from original, but with new names
makenewregs="cp "+freqpref+"_regions "+newpref+"_regions"
r=os.system(makenewregs)

makenewchrs="cp "+freqpref+"_chrs "+newpref+"_chrs"
r=os.system(makenewchrs)




