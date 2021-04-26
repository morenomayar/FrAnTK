
from __future__ import division
import sys
import os
import subprocess

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["bamfile", "freqpref", "trim", "MinMQ", "MinBQ", "indname", "popname", "regionsfile"]

argvect=sys.argv[1:]
arghash={}
for i in argvect:
	if "=" not in i:
		print(i+" is not a valid arg sntx\n")
		sys.exit(1)
	if len(i.split("="))!=2:
		print(i+" is not a valid arg, use arg=val")
		sys.exit(1)
	argname=i.split("=")[0].split()[0]
	if argname not in allowedpars:
		print(i+" is not a valid arg\n")
		sys.exit(1)
	arghash[argname]=i.split("=")[1].split()[0]

try:
	regionsfile=arghash["regionsfile"]
except:
	arghash["regionsfile"]=arghash["freqpref"]+"_regions"
	regionsfile=arghash["freqpref"]+"_regions"

for i in allowedpars:
	if i not in arghash:
		print(i+" is missing\n")
		sys.exit(1)

bamfile=arghash["bamfile"]
freqpref=arghash["freqpref"]
trim=arghash["trim"]
MinMQ=arghash["MinMQ"]
MinBQ=arghash["MinBQ"]
indname=arghash["indname"]
popname=arghash["popname"]

try:
	trim=int(trim)
	MinMQ=int(MinMQ)
	MinBQ=int(MinBQ)
except:
	print("trim, minmq and minbq should be int\n")

freqprefbase=freqpref.split("/")[len(freqpref.split("/"))-1]

cmd="EndQTrimPipe.pl -i "+bamfile+" -trim "+str(trim)+" -bed "+regionsfile+" - | samtools mpileup -Bl "+regionsfile+" -q "+str(MinMQ)+" -Q "+str(MinBQ)+" - | RandomizePU.pl - "
#cmd="perl EndQTrimPipe.pl -i "+bamfile+" -trim "+str(trim)+" -bed "+freqpref+"_regions - | samtools mpileup -Bl "+freqpref+"_regions -q "+str(MinMQ)+" -Q "+str(MinBQ)+" - | perl RandomizePU.pl - "

O=subprocess.Popen("gzip -c > "+ indname+"_"+popname+"_"+freqprefbase+"_spu.gz", shell=True, stdin=subprocess.PIPE)
b=open(regionsfile, "r")
bl=b.readline().split()
bn=str(bl[0])+"_"+str(bl[2])
a1=bl[3]
a2=bl[4]
nsites=0
with subprocess.Popen(cmd,stdout=subprocess.PIPE, shell=True).stdout as f:
	for line in f:
		if pyver==2:
			l=line.split()
		else:
			l=line.decode().split()
		ln=l[0]
		la=l[1]
		while ln != bn:
			if pyver==2:
				O.stdin.write("None\n")
			else:
				O.stdin.write(("None\n").encode())
			nsites+=1
			bl=b.readline().split()
			bn=str(bl[0])+"_"+str(bl[2])
			a1=bl[3]
			a2=bl[4]
		if la==a1:
			if pyver==2:
				O.stdin.write("1\n")
			else:
				O.stdin.write(("1\n").encode())
			nsites+=1
		elif la==a2:
			if pyver==2:
				O.stdin.write("0\n")
			else:
				O.stdin.write(("0\n").encode())
			nsites+=1
		else:
			if pyver==2:
				O.stdin.write("None\n")
			else:
				O.stdin.write(("None\n").encode())
			nsites+=1
		try:
			bl=b.readline().split()
			bn=str(bl[0])+"_"+str(bl[2])			
			a1=bl[3]
			a2=bl[4]
		except:
			#print "entre"
			break

while bl != []:
	if pyver==2:
		O.stdin.write("None\n")
	else:
		O.stdin.write(("None\n").encode())
	nsites+=1
	bl=b.readline().split()

b.close()
O.communicate()

