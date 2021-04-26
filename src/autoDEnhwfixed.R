
#####ArgParse#####

args<-commandArgs(TRUE)

#names of recognized arguments
PossibleArgNames<-c("freqpref", "h1", "h2", "h3", "h4", "catfile", "legfile", "nthr", "bsize", "resfile", "rmtrns")

USAGE<-("
	Compute D-statistics (or f4-statistics) of the form D(h1, h2; h3, h4), using (allele frequency) definition and normalization in Patterson et al, 2012 (Genetics). Fix three of h1,h2,h3,h4 and loop over the remaining one. 
	Arguments should be passed as: arg=value. Make sure that there is no spaces between arg, =, and value.
	freqpref\tSTR\tPrefix from prefix_freqs.gz file created by BuildFreqs.py. Required.
	h1\tSTR\tName of population to set in h1.
	h2\tSTR\tName of population to set in h2.
	h3\tSTR\tName of population to set in h3.
	h4\tSTR\tName of population to set in h4. Required
	NOTE: Exactly two of h1,h2,h3 are required. 
	catfile\tSTR\tTab separated file specifying a category for each population (used for plotting). Default category 'Other' is assigned to populations not present in this file. Contents in catfile will override those in the Cat column of resfile. Format: population category.  
	legfile\tSTR\tTab separated file specifying color and pch for each category (used for plotting and legend). Random colors and pch are chosen for categories not present in the file. Format: category outercolor innercolor pch. (Use R color names and pch numbers) .
	nthr\tINT\tNumber of threads to use. Each test is sent to a different cpu. Use many cores! =9 Default: 1
	bsize\tINT\tSize of each block used for weighted block jacknife (in nucleotides). Default: 5000000
	resfile\tSTR\tThis program creates a txt output file. Edit that file and feed it here for replotting after editing txt output, catfile or legfile. This option overrides all others. 
	rmtrns\t(0 | 1)\tInclude or remove transition SNPs. Default: 0
	f4\t(0 | 1)\tCompute f4-statistics (1) instead of D-statistics (0). Default: 0
	Sample call: autoDwfixed.R freqpref=prefix h1=popnameh1 h2=popnameh2 nthr=50 bsize=3500000 catfile=categories legfile=legendvalues rmtrns=1 f4=0
")

#fill hash with NAs
ArgHash<<-new.env()
for(i in PossibleArgNames){
	ArgHash[[i]]<-NA
}

#default values for bsize and nthreads
ArgHash[["bsize"]]<-5000000
ArgHash[["nthr"]]<-1
ArgHash[["rmtrns"]]<-0

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

#there should be at least one non-specified pop, which will be variable
if(sum(!is.na(ArgHash[["h1"]]), !is.na(ArgHash[["h2"]]), !is.na(ArgHash[["h3"]]))!=2 & is.na(ArgHash[["resfile"]])){
	message(USAGE)
	stop("exactly 2 pops should be specified from h1, h2, h3.")
}

if(is.na(ArgHash[["h4"]]) & is.na(ArgHash[["resfile"]])){
	message(USAGE)
	stop("h4 is required. Invert the ordering of the test (h1,h2 <-> h3,h4) to get the equivalent of varying h4.")
}

#retrieve arg vals from arghash, so they are used in an indep variable
freqpref<-ArgHash[["freqpref"]]
h1<-ArgHash[["h1"]]
h2<-ArgHash[["h2"]]
h3<-ArgHash[["h3"]]
h4<-ArgHash[["h4"]]
catfile<-ArgHash[["catfile"]]
legfile<-ArgHash[["legfile"]]
nthr<-as.integer(ArgHash[["nthr"]])
bsize<-as.integer(ArgHash[["bsize"]])
resfile<-ArgHash[["resfile"]]
rmtrns<-as.integer(ArgHash[["rmtrns"]])

if((is.na(freqpref) | is.na(h4)) & is.na(ArgHash[["resfile"]])){
	message(USAGE)
	stop("freqpref and h4 are required")
}

if(!(rmtrns==1 | rmtrns==0)){
	message(USAGE)
	stop("rmtrns can only be 0 or 1")
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
	seph4s<-strsplit(h4, ",")[[1]]
	if(sum(seph4s %in% pops)!=length(seph4s) & !is.na(h4)){
		message(USAGE)
		stop("h4 not in pops")
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

	Dpref<-"D_"

	Dprog<-"getEnhD.py"
	if(rmtrns==1){
		Dprog<-"getEnhDNT.py"
	}

	freqprefbase<-strsplit(freqpref, "/")[[1]][length(strsplit(freqpref, "/")[[1]])]

	#build cmds based on what is variable
	#if h1 is variable, h2 h3 and h4 are fixed
	if(is.na(h1)){
		a<-pops
		a<-a[pops!=h2 & pops!=h3 & !(pops %in% seph4s)]
		cats<-cats[pops!=h2 & pops!=h3 & !(pops %in% seph4s)]
		resdir<-paste(sep="_", Dpref, "h1", h2, h3, h4, freqprefbase, RandExt)
		cmds<-paste(Dprog, " h1=", a, " h2=", h2, " h3=", h3, " h4=", h4, " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep="")
		#filenames will have the names of the output files that are produced
		filenames<-paste(sep="", resdir, "/", Dpref, freqprefbase, "_", a, "_", h2, "_", h3, "_", h4, ".txt")
	}

	#if h2 is variable, h1 h3 and h4 are fixed
	if(is.na(h2)){
		a<-pops
		a<-a[pops!=h1 & pops!=h3 & !(pops %in% seph4s)]
		cats<-cats[pops!=h1 & pops!=h3 & !(pops %in% seph4s)]
		resdir<-paste(sep="_", Dpref, "h2", h1, h3, h4, freqprefbase, RandExt)
		cmds<-paste(Dprog, " h2=", a, " h1=", h1, " h3=", h3, " h4=", h4, " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep="")
		#filenames will have the names of the output files that are produced
		filenames<-paste(sep="", resdir, "/", Dpref, freqprefbase, "_", h1, "_", a, "_", h3, "_", h4, ".txt")
	}

	#if h3 is variable, h1 h2 and h4 are fixed
	if(is.na(h3)){
		a<-pops
		a<-a[pops!=h1 & pops!=h2 & !(pops %in% seph4s)]
		cats<-cats[pops!=h1 & pops!=h2 & !(pops %in% seph4s)]
		resdir<-paste(sep="_", Dpref, "h3", h1, h2, h4, freqprefbase, RandExt)
		cmds<-paste(Dprog, " h1=", h1, " h2=", h2, " h3=", a, " h4=", h4, " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep="")
		#filenames will have the names of the output files that are produced
		filenames<-paste(sep="", resdir, "/", Dpref, freqprefbase, "_", h1, "_", h2, "_", a, "_", h4, ".txt")
	}

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
	#resfile has the name of the concatenated results, will be used below, to enter the plotting section
	resfile<-paste(sep="", "Results_", resdir, ".txt")
	write.table(restab, file=resfile, quote=F, sep="\t", col.names=T, row.names=F)
}

#jump to here if resfile is specified
if(!is.na(resfile)){
	message("Plotting from file... ")
	Results<-read.table(resfile, as.is=T, h=T)
	Results<-Results[Results$Plot,]
	Results$D<-as.numeric(Results$D)
	Results$SE<-as.numeric(Results$SE)
	Results$Z<-as.numeric(Results$Z)
	Results$nSNPs<-as.numeric(Results$nSNPs)
	Results$nBlocks<-as.numeric(Results$nBlocks)
	Results$Plot<-as.logical(Results$Plot)
	Results<-Results[rowSums(is.na(Results))==0,]
	Fixed<-unname(unlist(lapply(apply(Results[,1:4], 2, table), length))==1)
	Upper<-Results$D+Results$SE
	Lower<-Results$D-Results$SE
	#This is a 99.9ci
	Upper3<-Results$D+qnorm(1-(.001)/2)*Results$SE
	Lower3<-Results$D-qnorm(1-(.001)/2)*Results$SE
	Results<-cbind(Results, Upper, Lower, Upper3, Lower3)

	#find out if there is a fixed column
	if(identical(Fixed, c(F,T,T,T))){
		Fixed1<-Results[1,2]
		Fixed2<-Results[1,3]
		Fixed3<-Results[1,4]
		MainText<-paste(sep="", "D(H1, ", Fixed1, "; ", Fixed2, ", ",  "Pool", ")\n(", Fixed1, ",", Fixed2, ") <---> (H1,", Fixed2, ")")
		variablecol<-Results$H1pop
		xlabtext<-"H1"
		TreeLabs<-gsub("_", "", c("H1", Fixed1, Fixed2, Fixed3))
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$H1, levels=Results$H1[sort.list(Results$D)])
		names(Results)<-gsub("H1pop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$D)])
	}else if(identical(Fixed, c(T,F,T,T))){
		Fixed1<-Results[1,1]
		Fixed2<-Results[1,3]
		Fixed3<-Results[1,4]
		MainText<-paste(sep="", "D(", Fixed1, ", H2; ", Fixed2, ", ",  "Pool", ")\n(", "H2", ",", Fixed2, ") <---> (", Fixed1, ",", Fixed2, ")")
		variablecol<-Results$H2pop
		xlabtext<-"H2"
		TreeLabs<-gsub("_", "", c(Fixed1, "H2", Fixed2, Fixed3))
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$H2, levels=Results$H2[sort.list(Results$D)])
		names(Results)<-gsub("H2pop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$D)])
	}else if(identical(Fixed, c(T,T,F,T))){
		Fixed1<-Results[1,1]
		Fixed2<-Results[1,2]
		Fixed3<-Results[1,4]
		MainText<-paste(sep="", "D(", Fixed1, ", ", Fixed2, "; H3, ", "Pool", ")\n(", Fixed2, ",", "H3", ") <---> (", Fixed1, ",", "H3", ")")
		variablecol<-Results$H3pop
		xlabtext<-"H3"
		TreeLabs<-gsub("_", "", c(Fixed1, Fixed2, "H3", Fixed3))
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$H3, levels=Results$H3[sort.list(Results$D)])
		names(Results)<-gsub("H3pop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$D)])
	}else{
		message(USAGE)
		stop("There are less than 2 fixed columns in the results, or variable outgroup. ")
	}

	#read category file to give a category to the results, otherwise, all outputs will have "Other" in category
	if(!is.na(catfile)){
		a<-read.table(catfile, as.is=T)
		for(i in seq_along(variablecol)){
			if(variablecol[i] %in% a[,1]){
				Results$Cat[i]<-a[which(a[,1]==variablecol[i])[1],2]
			}
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

	F4PLOT<-ggplot(Results, aes(x=D, y=X, color=Cat, fill=Cat, shape=Cat))+
	geom_vline(xintercept=0, colour="gray", linetype = "longdash")+
	#height used to be 0.6 and 0.3, but now they are 0
	#geom_errorbarh(aes(xmin=Lower, xmax=Upper), show.legend=F, height=0)+
	geom_errorbarh(aes(xmin=Lower3, xmax=Upper3), show.legend=F, height=0)+
	geom_point(aes(shape=Cat), size=2.7, show.legend=ifelse(length(table(Results$Cat))==1, F, T))+
	theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.key=element_rect(fill="white", colour="white"))+
	#guides(col = guide_legend(ncol = 1))+
	xlab(MainText)+
	ylab(xlabtext)+
	xlim(-(max(abs(Results$Lower3))+.1*max(abs(Results$Lower3))), (max(abs(Results$Upper3))+.1*max(abs(Results$Upper3))))+
	#scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=ceiling(length(CatsInPlot)/21)))+
	scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=1))+
	scale_color_manual(values=OuterColVals)+
	scale_shape_manual(values=PchVals)+
	geom_text(data=Results[-Results$D<=0, ], aes(x=Lower3, y=X, label=z), hjust=1.05, show.legend=F, size=2.7)+
	geom_text(data=Results[-Results$D>0, ], aes(x=Upper3, y=X, label=z), hjust=-.05, show.legend=F, size=2.7)

	pdf(paste(sep="", resfile, ".pdf"), height=ifelse(dim(Results)[1]<=30, dim(Results)[1]*8/52, dim(Results)[1]*6.5/52), width=ifelse(length(table(Results$Cat))==1, 5.5, 7.5), useDingbats=F)
	print(F4PLOT)
	dev.off()
	message(paste(sep="", "Plot located in: ", paste(sep="", resfile, ".pdf")))
}

