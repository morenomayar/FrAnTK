
Args<-commandArgs(T)
freqpref=Args[1]
bamfile=Args[2]
flag=Args[3]

#pars for testing
#freqpref="test2frq"
#bamfile<-"Pd1.bam"

bamcontigs<-gsub("^SN:", "", read.table(pipe(paste0("samtools view -H ", bamfile, " | grep -P \"^@SQ\t\"")), as.is=T)[,2])
chrs<-read.table(paste0(freqpref, "_chrs"), as.is=T)
freqscontigs<-chrs[,1]

#
#bad order example
#bamcontigs<-sort(bamcontigs)
#bamcontigs<-sort(bamcontigs[-2])
#

bamcontigsin<-as.character(bamcontigs[bamcontigs %in% freqscontigs])
freqscontigsin<-as.character(freqscontigs[freqscontigs %in% bamcontigs])
freqscontigsout<-as.character(freqscontigs[!(freqscontigs %in% bamcontigs)])

#check
if(flag=="check"){
	if(length(bamcontigsin)==0 & length(freqscontigsin)==0){
		# check at least one overlap
		message("There is no overlap between bam and freqs contigs. Make sure contigs follow the same nomenclature in both files. ")
		quit(save="no", status=1)
	}else if(identical(bamcontigsin, freqscontigsin)){
		# order ok
		message("Bam and freqs contigs are in the same order. ")
		quit(save="no", status=0)
	}else{
		message(paste0("Bam and freqs contigs are NOT in the same order. Try running:\nCheckContigOrder.R ", freqpref, " ", bamfile, " sort"))
		print(data.frame(freqscontigsin, bamcontigsin))
		quit(save="no", status=2)
	}
}else if(flag=="sort"){
	getRandString<-function(len=10){
		return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
	}
	RandExt<-getRandString()
	newfreqpref=paste0(freqpref, "_", RandExt)

	###
	#remove after testing
	#system(paste0("rm ", newfreqpref, "_freqs.gz"))
	#system(paste0("rm ", newfreqpref, "_regions"))
	#system(paste0("rm ", newfreqpref, "_chrs"))
	#system(paste0("rm ", newfreqpref, "_pop"))
	###

	#for each contig (in the bam order)
	for(i in c(bamcontigsin, freqscontigsout)){
		#get start and end pos in freqs
		st<-chrs[chrs[,1]==i, 4]
		en<-chrs[chrs[,1]==i, 5]
		#subset freqs file and append to new
		subsetcmd<-paste0("zcat ", freqpref, "_freqs.gz | sed -n '", st, ",", en, "p;", en, "q' | gzip -c >> ", newfreqpref, "_freqs.gz")
		system(subsetcmd)
		#subset regions file and append to new
		subsetcmd<-paste0("sed -n '", st, ",", en, "p;", en, "q' ", freqpref, "_regions >> ", newfreqpref, "_regions")
		system(subsetcmd)
	}

	#create new chrs
	#get number of snps per contig and lengths and firstpos
	nsnps<-chrs[,5]-chrs[,4]+1
	chrls<-chrs[,3]
	firstpos<-chrs[,2]

	#sort nsnps according to new contig order
	sortednsnps<-NULL
	sortedchrls<-NULL
	sortedfirstpos<-NULL
	for(i in c(bamcontigsin, freqscontigsout)){
		sortednsnps<-c(sortednsnps, nsnps[chrs[,1]==i])
		sortedchrls<-c(sortedchrls, chrs[chrs[,1]==i,3])
		sortedfirstpos<-c(sortedfirstpos, chrs[chrs[,1]==i,2])
	}

	#create coordinates
	newends<-cumsum(sortednsnps)
	newstarts<-c(1, newends+1)
	newstarts<-newstarts[-length(newstarts)]

	#write new chrs
	newchrs<-data.frame(c(bamcontigsin, freqscontigsout), sortedfirstpos, sortedchrls, newstarts, newends, stringsAsFactors=F)
	write.table(newchrs, file=paste0(newfreqpref, "_chrs"), quote=F, sep="\t", col.names=F, row.names=F)

	#copy _pop file
	system(paste0("cp ", freqpref, "_pop ", newfreqpref, "_pop"))

	message(paste0("New files are named ", newfreqpref, "*"))

	message("New order is: ")
	print(bamcontigsin)

	if(length(freqscontigsout)>0){
		message("These contigs are not present in bam, so we sent them to the end of the file")
		print(freqscontigsout)
	}

}else{
	stop("\nUSAGE: \nCheckContigOrder.R freqpref bamfile check\nCheckContigOrder.R freqpref bamfile sort")
}


