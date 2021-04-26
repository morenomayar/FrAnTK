
from __future__ import division
import sys
import os
import subprocess

pyver=2
if sys.version_info[0] == 3:
	xrange = range
	pyver=3

allowedpars=["bamfile", "plinkpref", "trim", "MinMQ", "MinBQ", "indname", "popname", "doCns", "regionsfile"]

argvect=sys.argv[1:]
arghash={}
arghash["doCns"]="F"
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
	arghash["regionsfile"]=arghash["plinkpref"]+"_regions"
	regionsfile=arghash["plinkpref"]+"_regions"

for i in allowedpars:
	if i not in arghash:
		print(i+" is missing\n")
		sys.exit(1)

bamfile=arghash["bamfile"]
plinkpref=arghash["plinkpref"]
trim=arghash["trim"]
MinMQ=arghash["MinMQ"]
MinBQ=arghash["MinBQ"]
indname=arghash["indname"]
popname=arghash["popname"]
doCns=arghash["doCns"]

try:
	trim=int(trim)
	MinMQ=int(MinMQ)
	MinBQ=int(MinBQ)
except:
	print("trim, minmq and minbq should be int\n")

plinkprefbase=plinkpref.split("/")[len(plinkpref.split("/"))-1]

if doCns=="T":
	cmd="EndQTrimPipe.pl -i "+bamfile+" -trim "+str(trim)+" -bed "+regionsfile+" - | samtools mpileup -RBl "+regionsfile+" -q "+str(MinMQ)+" -Q "+str(MinBQ)+" - | ConsensusPU.pl - "
else:
	cmd="EndQTrimPipe.pl -i "+bamfile+" -trim "+str(trim)+" -bed "+regionsfile+" - | samtools mpileup -RBl "+regionsfile+" -q "+str(MinMQ)+" -Q "+str(MinBQ)+" - | RandomizePU.pl - "

#cmd="perl EndQTrimPipe.pl -i "+bamfile+" -trim "+str(trim)+" -bed "+plinkpref+"_regions - | samtools mpileup -Bl "+plinkpref+"_regions -q "+str(MinMQ)+" -Q "+str(MinBQ)+" - | perl RandomizePU.pl - "

#O=subprocess.Popen("gzip -c > "+ indname+"_"+popname+"_"+plinkprefbase+"_spu.gz", shell=True, stdin=subprocess.PIPE)
O=open(indname+"_"+popname+"_"+plinkprefbase+".tped", 'w')
b=open(regionsfile, "r")
p=open(plinkpref+".bim", "r")
bl=b.readline().split()
pl=p.readline().split()
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
			#O.stdin.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			#O.stdin.write("0\t0\n")
			O.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			O.write("0\t0\n")
			nsites+=1
			bl=b.readline().split()
			pl=p.readline().split()
			bn=str(bl[0])+"_"+str(bl[2])
			a1=bl[3]
			a2=bl[4]
		if la==a1:
			#O.stdin.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			#O.stdin.write(a1+"\t"+a1+"\n")
			O.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			O.write(a1+"\t"+a1+"\n")
			nsites+=1
		elif la==a2:
			#O.stdin.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			#O.stdin.write(a2+"\t"+a2+"\n")
			O.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			O.write(a2+"\t"+a2+"\n")
			nsites+=1
		else:
			#O.stdin.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			#O.stdin.write("0\t0\n")
			O.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
			O.write("0\t0\n")
			nsites+=1
		try:
			bl=b.readline().split()
			pl=p.readline().split()
			bn=str(bl[0])+"_"+str(bl[2])			
			a1=bl[3]
			a2=bl[4]
		except:
			#print "entre"
			break

while bl != []:
	#O.stdin.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
	#O.stdin.write("0\t0\n")
	O.write(pl[0]+"\t"+pl[1]+"\t"+pl[2]+"\t"+pl[3]+"\t")
	O.write("0\t0\n")
	nsites+=1
	bl=b.readline().split()
	pl=p.readline().split()

b.close()
#O.communicate()
O.close()

O=open(indname+"_"+popname+"_"+plinkprefbase+".tfam", 'w')
O.write(popname+"\t"+indname+"\t0\t0\t0\t1\n")
O.close()

cmd="plink --allow-extra-chr --tfile "+indname+"_"+popname+"_"+plinkprefbase+" --make-bed --out "+indname+"_"+popname+"_"+plinkprefbase
r=os.system(cmd)
cmd="rm "+indname+"_"+popname+"_"+plinkprefbase+".tped"
r=os.system(cmd)
cmd="rm "+indname+"_"+popname+"_"+plinkprefbase+".tfam"
r=os.system(cmd)
