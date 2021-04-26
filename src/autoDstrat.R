
#####ArgParse#####

args<-commandArgs(TRUE)

options(scipen=999)

#names of recognized arguments
PossibleArgNames<-c("freqpref", "h1", "h2", "h3", "h4", "anc", "nthr", "bsize", "resfile", "rmtrns", "minder", "maxder", "dstep", "mincops", "catfile", "legfile", "minD", "maxD", "strattype")

USAGE<-("
	Compute D statistics of the form D(h1, h2; h3, h4) on sites with a given number of derived alleles, using (allele frequency) definition and normalization in Patterson et al, 2012 (Genetics). We compute D for sites with [minder, maxder] (step size=dstep) derived alleles, where h1,h2,h3,h4 have mincops1,mincops2,mincops3,mincops4 copies. 
	Arguments should be passed as: arg=value. Make sure that there is no spaces between arg, =, and value.
	freqpref\tSTR\tPrefix from prefix_freqs.gz file created by BuildFreqs.py. Required.
	h1\tSTR\tName of population to set in h1. Required
	h2\tSTR\tName of population to set in h2. Required
	h3\tSTR\tName of population to set in h3. Required
	h4\tSTR\tName of population to set in h4. Required
	anc\tSTR\tName of population with the ancestral allele (should be a single hap ind and should be present in prefix_pop). Required. 
	minder\tInt. Compute D on sites with >=minder derived alleles. Default=2.
	maxder\tInt. Compute D on sites with <=maxder derived alleles. Default=10.
	dstep\tInt. Step size for range (minder, maxder]. Default=1.
	mincops\tString. Comma-separated (no spaces) list of minimum number of total copies in h1,h2,h3,h4, e.g., 1,1,1,1 for no filter. Default=1,1,1,1. 
	strattype\t(1 | 2). Count derived alleles in h1,h2,h3,h4 (1) or in the complete panel (2). Default=1. 
	minD\tFloat. Lower ylim. Should be used together with maxD. Default=NA.
	maxD\tFloat. Upper ylim. Should be used together with minD. Default=NA.
	catfile\tSTR\tTab separated file specifying a category for each population (used for plotting). Default category 'Other' is assigned to populations not present in this file. Contents in catfile will override those in the Cat column of resfile. Format: population category. IMPORTANT: Since we are not looping over different populations, we use the category of H2 to assign colors and symbols.  
	legfile\tSTR\tTab separated file specifying color and pch for each category (used for plotting and legend). Random colors and pch are chosen for categories not present in the file. Format: category outercolor innercolor pch. (Use R color names and pch numbers) .
	nthr\tINT\tNumber of threads to use. Each test is sent to a different cpu. Use many cores! =9 Default: 1
	bsize\tINT\tSize of each block used for weighted block jacknife (in nucleotides). Default: 5000000
	resfile\tSTR\tThis program creates a txt output file. Edit that file and feed it here for replotting after editing txt output, catfile or legfile. This option overrides all others. 
	rmtrns\t(0 | 1)\tInclude or remove transition SNPs. Default: 0
	IMPORTANT: Make sure that pseudo-haploid populations were specified correctly in the _clust file before running BuildFreqs.py. Otherwise, derived alleles in pseudo-haploid pops will be counted twice. 
	Sample call: autof4subtr.R freqpref=prefix h1=popnameh1 h2=popnameh2 h3=popnameh3 h4=popnameh4 x=popnamex nthr=50 bsize=3500000 catfile=categories legfile=legendvalues rmtrns=1
")

#fill hash with NAs
ArgHash<<-new.env()
for(i in PossibleArgNames){
	ArgHash[[i]]<-NA
}

#default values for bsize, nthreads, derived alleles and missingness
ArgHash[["bsize"]]<-5000000
ArgHash[["nthr"]]<-1
ArgHash[["rmtrns"]]<-0
ArgHash[["strattype"]]<-1
ArgHash[["minder"]]<-2
ArgHash[["maxder"]]<-10
ArgHash[["dstep"]]<-1
ArgHash[["mincops"]]<-"1,1,1,1"

#remove spaces from argument names and vals
argnames<-gsub(" ", "", unlist(lapply(strsplit(args, "="), "[[", 1)))
#message(argnames)
argvals<-gsub(" ", "", unlist(lapply(strsplit(args, "="), "[[", 2)))
#message(argvals)

#check that all argnames are recognized
if(sum(!(argnames%in%PossibleArgNames))>0){
	message(USAGE)
	stop(paste(collapse=" ", c("Unrecognized arguments", argnames[!(argnames%in%PossibleArgNames)])))
}

#assign each arg val to an arg name in the arg hash
for(i in seq_along(argnames)){
	ArgHash[[argnames[i]]]<-argvals[i]
}

#all pops should be specified unless we have a resfile
if(sum(!is.na(ArgHash[["h1"]]), !is.na(ArgHash[["h2"]]), !is.na(ArgHash[["h3"]]), !is.na(ArgHash[["h4"]]), !is.na(ArgHash[["anc"]]))!=5 & is.na(ArgHash[["resfile"]])){
	message(USAGE)
	stop("h1, h2, h3, h4 and anc are required.")
}

#retrieve arg vals from arghash, so they are used in an indep variable
freqpref<-ArgHash[["freqpref"]]
h1<-ArgHash[["h1"]]
h2<-ArgHash[["h2"]]
h3<-ArgHash[["h3"]]
h4<-ArgHash[["h4"]]
anc<-ArgHash[["anc"]]
clustfile<-ArgHash[["clustfile"]]
catfile<-ArgHash[["catfile"]]
legfile<-ArgHash[["legfile"]]
nthr<-as.integer(ArgHash[["nthr"]])
bsize<-as.integer(ArgHash[["bsize"]])
resfile<-ArgHash[["resfile"]]
rmtrns<-as.integer(ArgHash[["rmtrns"]])
strattype <-as.integer(ArgHash[["strattype"]])

if(!is.na(ArgHash[["minder"]]) & !is.na(ArgHash[["maxder"]]) & !is.na(ArgHash[["dstep"]])){
	minder<-as.integer(ArgHash[["minder"]])
	maxder<-as.integer(ArgHash[["maxder"]])
	dstep<-as.numeric(ArgHash[["dstep"]])
	#check if min maxder could be cast as int
	if(is.na(minder) | is.na(maxder) | is.na(dstep)){
		message(USAGE)
		stop("minder, maxder and dstep should be numeric. ")
	}

	#check minder and maxder lower bounds
	if(minder<1 | maxder<2 | dstep<1){
		message(USAGE)
		stop("minder should be >=1, maxder should be >=2 and dstep should be >=1. ")
	}

	#check minder<maxder
	if(minder>=maxder){
		message(USAGE)
		stop("minder should be <maxder. ")
	}

	#check dstep is within the range
	pdiff<-maxder-minder
	if(dstep>pdiff){
		message(USAGE)
		stop("dstep should be <(maxder-minder). ")
	}

	#create range
	st<-seq(minder, maxder, dstep)
	en<-c(st[-1], st[length(st)]+dstep)-1
	st<-st-1
	en[length(en)]<-maxder

	minderrange<-st
	maxderrange<-en

}else{
	message(USAGE)
	stop("minder, maxder and dstep are required.")
}

	#retrieve mincopps array
	mincops<-ArgHash[["mincops"]]
	mincopsv<-strsplit(mincops, ",")[[1]]

	#check mincopsv has 4 elements
	if(length(mincopsv)!=4){
		message(USAGE)
		stop("mincops should be a comma-separated (no spaces) list of minimum number of total copies in h1,h2,h3,h4, e.g., 1,1,1,1 for no filter. ")
	}

	#check if mincopsv has floats
	mincopsv<-as.integer(mincopsv)

	if(any(is.na(mincopsv))){
		message(USAGE)
		stop("mincops should be a comma-separated (no spaces) list of minimum number of total copies in h1,h2,h3,h4, e.g., 1,1,1,1 for no filter. ")
	}

	#check if mincopsv are >0
	if(any(mincopsv<=0)){
		message(USAGE)
		stop("mincops should be >0. ")
	}

flt<-F
if(!is.na(ArgHash[["minD"]]) | !is.na(ArgHash[["maxD"]])){
	if(!is.na(ArgHash[["minD"]]) & !is.na(ArgHash[["maxD"]])){
		minD<-as.numeric(ArgHash[["minD"]])
		maxD<-as.numeric(ArgHash[["maxD"]])
		if(is.na(minD) | is.na(maxD)){
			message(USAGE)
			stop("minD and maxD should be floats.")
		}
		flt<-T
	}else{
		message(USAGE)
		stop("minD and maxD should be used together.")
	}
}

#check if we have required args or resfile
if((is.na(freqpref) | is.na(h1) | is.na(h2) | is.na(h3) | is.na(h4) | is.na(anc)) & is.na(ArgHash[["resfile"]])){
	message(USAGE)
	stop("freqpref, h1, h2, h3, h4, anc are required")
}

if(!(rmtrns==1 | rmtrns==0)){
	message(USAGE)
	stop("rmtrns can only be 0 or 1")
}

if(!(strattype==1 | strattype==2)){
	message(USAGE)
	stop("strattype can only be 1 or 2")
}


#only do this if resfile was not specified, otherwise just go below and do the plotting
if(is.na(resfile)){
	#read in pop names, based on the freq prefix
	pops<-read.table(paste(freqpref, "_pop", sep=""), as.is=T)[,1]
	#initialize the default categories for each variable pop
	cats<-rep("Other", length(pops))

	#check that h1 is in pops
	if(!(h1 %in% pops) & !is.na(h1)){
		message(USAGE)
		stop("h1 not in pops")
	}

	#check that h2 is in pops
	if(!(h2 %in% pops) & !is.na(h2)){
		message(USAGE)
		stop("h2 not in pops")
	}

	#check that h3 is in popfile
	if(!(h3 %in% pops) & !is.na(h3)){
		message(USAGE)
		stop("h3 not in pops")
	}

	#check that h4 is in popfile
	if(!(h4 %in% pops) & !is.na(h4)){
		message(USAGE)
		stop("h4 not in pops")
	}

	#check that anc is in popfile
	if(!(anc %in% pops) & !is.na(anc)){
		message(USAGE)
		stop("anc not in pops")
	}

	#check that nthr is an int
	if(is.na(nthr)){
		message(USAGE)
		stop("nthr should be int")
	}

	#check that bszie is an int
	if(is.na(bsize)){
		message(USAGE)
		stop("bsize should be int")
	}

	#read category file to give a category to the outputs, otherwise, all outputs will have "Other" in category
	if(!is.na(catfile)){
		a<-read.table(catfile, as.is=T)
		for(i in seq_along(pops)){
			if(pops[i] %in% a[,1]){
				cats[i]<-a[which(a[,1]==pops[i])[1],2]
			}
		}
	}

	#fix weird symbols in cats
	cats<-gsub(">", "\\\\>", cats)
	cats<-gsub("<", "\\\\<", cats)

	#build cmds
	#use random string to differentiate runs
	getRandString<-function(len=6){
		return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
	}
	RandExt<-getRandString()

	if(strattype==1){
		Dpref<-"Dstrat_"
		Dpref2<-"D_"
		if(rmtrns==1){
			Dpref<-"Dstrat_"
			Dpref2<-"D_"
		}

		Dprog<-"getDstrat.py"
		Dprog2<-"getD.py"
		if(rmtrns==1){
			Dprog<-"getDstratNT.py"
			Dprog2<-"getDNT.py"
		}
	}else if(strattype==2){
		Dpref<-"Dstrat2_"
		Dpref2<-"D_"
		if(rmtrns==1){
			Dpref<-"Dstrat2_"
			Dpref2<-"D_"
		}

		Dprog<-"getDstrat2.py"
		Dprog2<-"getD.py"
		if(rmtrns==1){
			Dprog<-"getDstrat2NT.py"
			Dprog2<-"getDNT.py"
		}
	}else{
		message(USAGE)
		stop("strattype can only be 1 or 2")
	}

	freqprefbase<-strsplit(freqpref, "/")[[1]][length(strsplit(freqpref, "/")[[1]])]

	#build cmds
	cats<-cats[pops==h2]
	resdir<-paste(sep="_", Dpref, h1, h2, h3, h4, anc, as.character(minder), as.character(maxder), as.character(mincops), freqprefbase, RandExt)
	cmds<-paste(Dprog, " h1=", h1, " h2=", h2, " h3=", h3, " h4=", h4, " anc=", anc, " minder=", as.character(minderrange), " maxder=", as.character(maxderrange), " mincops=", as.character(mincops), " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep="")
	cmds<-c(cmds, paste(Dprog2, " h1=", h1, " h2=", h2, " h3=", h3, " h4=", h4, " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep=""))
	#filenames will have the names of the output files that are produced
	filenames<-paste(sep="", resdir, "/", Dpref, freqprefbase, "_", h1, "_", h2, "_", h3, "_", h4, "_", anc, "_", as.character(minderrange), "_", as.character(maxderrange), "_", as.character(mincops), ".txt")
	filenames<-c(filenames, paste(sep="", resdir, "/", Dpref2, freqprefbase, "_", h1, "_", h2, "_", h3, "_", h4, ".txt"))
	#print(filenames)

	message("Running tests... ")
	system(paste0("mkdir ", resdir))
	#run cmds
	library(doParallel)
	registerDoParallel(cores=nthr)

	foreach(i=cmds) %dopar%{
		system(i)
	}

	#write results
	#read in all produced output files
	restab<-read.table(pipe(paste("cat", paste(filenames, collapse=" "))), as.is=T)
	#give names to data frame
	names(restab)<-c("H1pop", "H2pop", "H3pop", "H4pop", "D", "SE", "Z", "nSNPs", "nBlocks", "Plot", "Cat", "nABBA", "nBABA")
	#get adm props
	maxders<-c(maxderrange, "ALLSITES")
	restab<-cbind.data.frame(restab, maxders, stringsAsFactors=F)
	#resfile has the name of the concatenated results, will be used below, to enter the plotting section
	resfile<-paste(sep="", "Results_", resdir, ".txt")
	write.table(restab, file=resfile, quote=F, sep="\t", col.names=T, row.names=F)
}

#jump to here if resfile is specified
if(!is.na(resfile)){
	message("Plotting from file... ")
	Results<-read.table(resfile, as.is=T, h=T)

	if(length(unique(Results$H1pop))!=1 | length(unique(Results$H2pop))!=1 | length(unique(Results$H3pop))!=1 | length(unique(Results$H4pop))!=1){
		message(USAGE)
		stop("All columns H1pop, H2pop, H3pop, H4pop should be fixed. ")
	}

	Results<-Results[Results$Plot,]
	Results$D<-as.numeric(Results$D)
	baseD<-NA
	if("ALLSITES" %in% Results$maxders){
		baseD<-Results[Results$maxders=="ALLSITES", ]$D
		Results<-Results[Results$maxders!="ALLSITES",]
	}
	Results$SE<-as.numeric(Results$SE)
	Results$Z<-as.numeric(Results$Z)
	Results$nSNPs<-as.numeric(Results$nSNPs)
	Results$nBlocks<-as.numeric(Results$nBlocks)
	Results$Plot<-as.logical(Results$Plot)
	Results$maxders<-as.numeric(Results$maxders)
	Results<-Results[rowSums(is.na(Results))==0,]
	Upper<-Results$D+Results$SE
	Lower<-Results$D-Results$SE
	#This is a 99.9ci
	Upper3<-Results$D+qnorm(1-(.001)/2)*Results$SE
	Lower3<-Results$D-qnorm(1-(.001)/2)*Results$SE
	Results<-cbind(Results, Upper, Lower, Upper3, Lower3)

	MainTextY<-paste(sep="", "Dstrat(", Results$H1pop[1], ",", Results$H2pop[1], ";", Results$H3pop[1], ",", Results$H4pop[1], ")\n(", Results$H2pop[1], ",", Results$H3pop[1], ")<--->(", Results$H1pop[1], ",", Results$H3pop[1], ")")

	if(strattype==1){
		MainText<-"Number of derived alleles"
	}else{
		MainText<-"Number of derived alleles in panel"
	}

	#read category file to give a category to the results, otherwise, all outputs will have "Other" in category
	if(!is.na(catfile)){
		a<-read.table(catfile, as.is=T)
		if(Results$H2pop[1] %in% a[,1]){
			Results$Cat<-a[which(a[,1]==Results$H2pop[1])[1],2]
		}
	}

	#read legend file or do random colors
	if(!is.na(legfile)){
		CatSect<-read.table(legfile, as.is=T)
		CatSect<-CatSect[,1:4]
		CatSect<-rbind(CatSect, c("Other", "purple1", "purple1", "1"))
		CatSect<-CatSect[!duplicated(CatSect[,1]),]
		names(CatSect)<-c("Cat", "OuterColor", "InnerColor", "Pch")
		CatSect$Pch<-as.numeric(CatSect$Pch)
	}else{
		CatSect<-(data.frame(Cat="Other", OuterColor="purple1", InnerColor="purple1", Pch=1, stringsAsFactors=F))
	}

	#fill up random colors, etc
	for(i in unique(Results$Cat)){
		if(!(i %in% CatSect$Cat)){
			randomcolor<-sample(colors()[c(32:140, 365:657)], 1)
			NewRow<-c(i, randomcolor, randomcolor, 1)
			CatSect<-rbind(CatSect, NewRow)
			message(paste(sep="", "Added auto row to categories: ", paste(NewRow, collapse=" ")))
		}
	}

	for(i in seq_along(CatSect$OuterColor)){
		if(!(CatSect$OuterColor[i] %in% colors())){
			message(sep="", CatSect$OuterColor[i], " is not a valid color. Picking a random color... ")
			CatSect$OuterColor[i]<-sample(colors()[c(32:140, 365:657)], 1)
		}
	}

	for(i in seq_along(CatSect$InnerColor)){
		if(!(CatSect$InnerColor[i] %in% colors())){
			message(sep="", CatSect$InnerColor[i], " is not a valid color. Picking a random color... ")
			CatSect$InnerColor[i]<-sample(colors()[c(32:140, 365:657)], 1)
		}
	}

	CatSect$Pch<-as.numeric(CatSect$Pch)

	for(i in seq_along(CatSect$Pch)){
		if(CatSect$Pch[i] > 25 | CatSect$Pch[i] < 0 | is.na(CatSect$Pch[i])){
			message(paste(sep="", CatSect$Pch[i], " is not a valid pch. Setting to default... "))
			CatSect$Pch[i]<-1
		}
	}

	#this is a function for getting transparent colors
	transcol<-Vectorize(function(x, a){
		col<-col2rgb(x)
		newcol<-rgb(col[1]/255, col[2]/255, col[3]/255, a)
		return(newcol)
	})

	#if pch is a filled symbol, and inner and outer are the same, give inner alpha .5
	totransp<-CatSect$Pch>=21 & CatSect$Pch<=25 & CatSect$InnerColor==CatSect$OuterColor
	if(sum(totransp)>0){
		CatSect$InnerColor[totransp]<-unname(transcol(CatSect$InnerColor[totransp], .5))
	}

	#if pch is a filled symbol, but inner and outer are different, give inner alpha .65
	totransp<-CatSect$Pch>=21 & CatSect$Pch<=25 & CatSect$InnerColor!=CatSect$OuterColor
	if(sum(totransp)>0){
		CatSect$InnerColor[totransp]<-unname(transcol(CatSect$InnerColor[totransp], .65))
	}

	OuterColVals<-CatSect$OuterColor
	names(OuterColVals)<-CatSect$Cat
	InnerColVals<-CatSect$InnerColor
	names(InnerColVals)<-CatSect$Cat
	PchVals<-CatSect$Pch
	names(PchVals)<-CatSect$Cat

	#Do the plot
	library("ggplot2")
	library(grid)

	#from the categories to be plotted, get them in order from CatSect
	CatsInPlot<-CatSect[CatSect$Cat %in% Results$Cat,1]
	#change Results$Cat order as a factor
	Results$Cat<-factor(Results$Cat, levels=CatsInPlot)

	#This is an additional column with a rounded Z-score, 
	#which will be printed in the plot
	Results$z<-round(abs(Results$Z),2)

	F4PLOT<-ggplot(Results, aes(x=maxders, y=D, color=Cat, fill=Cat, shape=Cat))+
	geom_hline(yintercept=0, colour="gray", linetype = "longdash")+
	#height used to be 0.6 and 0.3, but now they are 0
	#geom_errorbarh(aes(xmin=Lower, xmax=Upper), show.legend=F, height=0)+
	geom_errorbar(aes(ymin=Lower3, ymax=Upper3), show.legend=F, width=0)+
	geom_point(aes(shape=Cat), size=2.7, show.legend=ifelse(length(table(Results$Cat))==1, F, T))+
	theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.x=element_text(size=7.5), axis.title.y=element_text(size=7.5), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.key=element_rect(fill="white", colour="white"))+
	#guides(col = guide_legend(ncol = 1))+
	xlab(MainText)+
	ylab(MainTextY)+
	#scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=ceiling(length(CatsInPlot)/21)))+
	scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=1))+
	scale_color_manual(values=OuterColVals)+
	scale_shape_manual(values=PchVals)+
	geom_text(data=Results[-Results$D<=0, ], aes(y=Lower3, x=maxders, label=z), hjust=1.05, show.legend=F, size=2.7, angle=90)+
	geom_text(data=Results[-Results$D>0, ], aes(y=Upper3, x=maxders, label=z), hjust=-.05, show.legend=F, size=2.7, angle=90)

	if(!is.na(baseD)){
		F4PLOT<-F4PLOT+geom_hline(yintercept=baseD, linetype="dashed")
	}

	if(flt==T){
		F4PLOT<-F4PLOT+coord_cartesian(ylim=c(minD, maxD))
	}else{
		F4PLOT<-F4PLOT+coord_cartesian(ylim=c(-(max(abs(Results$Lower3))+.1*max(abs(Results$Lower3))), (max(abs(Results$Upper3))+.1*max(abs(Results$Upper3)))))
	}

	pdf(paste(sep="", resfile, ".pdf"), height=3.3, width=6, useDingbats=F)
	print(F4PLOT, newpage=FALSE)
	dev.off()
	message(paste(sep="", "Plot located in: ", paste(sep="", resfile, ".pdf")))
}








