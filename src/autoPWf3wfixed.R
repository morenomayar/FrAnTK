
#####ArgParse#####

args<-commandArgs(TRUE)

#names of recognized arguments
PossibleArgNames<-c("freqpref", "x", "y", "target", "catfile", "legfile", "nthr", "bsize", "rmtrns")

USAGE<-("
	Compute f3 statistics of the form f3(h1, h2; target), using definition and normalization in Patterson et al, 2012 (Genetics), for two samples and get pairwise comparisons. x and y will be fixed in h1 and h2 will be all other pops in the panel, target will also be fixed. 
	Arguments should be passed as: arg=value. Make sure that there is no spaces between arg, =, and value.
	freqpref\tSTR\tPrefix from prefix_freqs.gz file created by BuildFreqs.py. Required.
	x\tSTR\tName of population to set in h1 in the x-axis.
	y\tSTR\tName of population to set in h1 in the y-axis.
	target\tSTR\tName of population to set in target. Fixed
	catfile\tSTR\tTab separated file specifying a category for each population (used for plotting). Default category 'Other' is assigned to populations not present in this file. Contents in catfile will override those in the Cat column of resfile. Format: population category.  
	legfile\tSTR\tTab separated file specifying color and pch for each category (used for plotting and legend). Random colors and pch are chosen for categories not present in the file. Format: category outercolor innercolor pch. (Use R color names and pch numbers) .
	nthr\tINT\tNumber of threads to use. Each test is sent to a different cpu. Use many cores! =9 Default: 1
	bsize\tINT\tSize of each block used for weighted block jacknife (in nucleotides). Default: 5000000
	rmtrns\t(0 | 1)\tInclude or remove transition SNPs. Default: 0
	Sample call: autoPWf3wfixed.R freqpref=prefix x=firstpopnameh1 y=secondpopnameh1 nthr=50 bsize=3500000 catfile=categories legfile=legendvalues rmtrns=1
")

if(length(args)==0){
	message(USAGE)
	q("no")
}

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

#retrieve arg vals from arghash, so they are used in an indep variable
freqpref<-ArgHash[["freqpref"]]
x<-ArgHash[["x"]]
y<-ArgHash[["y"]]
target<-ArgHash[["target"]]
catfile<-ArgHash[["catfile"]]
legfile<-ArgHash[["legfile"]]
nthr<-as.integer(ArgHash[["nthr"]])
bsize<-as.integer(ArgHash[["bsize"]])
rmtrns<-as.integer(ArgHash[["rmtrns"]])

if(is.na(freqpref)){
	message(USAGE)
	stop("freqpref is required")
}

if(!(rmtrns==1 | rmtrns==0)){
	message(USAGE)
	stop("rmtrns can only be 0 or 1")
}


####

#Important to go into h1 conditions
h2<-NA
	#read in pop names, based on the freq prefix
	pops<-read.table(paste(freqpref, "_pop", sep=""), as.is=T)[,1]
	#initialize the default categories for each variable pop
	cats<-rep("Other", length(pops))

	#check that x is in pops
	if(!(x %in% pops) & !is.na(x)){
		message(USAGE)
		stop("x not in pops")
	}

	#check that y is in pops
	if(!(y %in% pops) & !is.na(y)){
		message(USAGE)
		stop("y not in pops")
	}

	#check that target is in popfile
	if(!(target %in% pops) & !is.na(target)){
		message(USAGE)
		stop("target not in pops")
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


#############
#Do x f3

h1<-x

	#build cmds
	#use random string to differentiate runs
	getRandString<-function(len=6){
		return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
	}
	RandExt<-getRandString()

	Dpref<-"f3_"

	Dprog<-"getf3.py"
	if(rmtrns==1){
		Dprog<-"getf3NT.py"
	}

	freqprefbase<-strsplit(freqpref, "/")[[1]][length(strsplit(freqpref, "/")[[1]])]

	#build cmds based on what is variable
	#if h1 is variable, h2 and target are fixed
	oricats<-cats
	if(is.na(h2)){
		a<-pops
		a<-a[pops!=h1 & pops!=target]
		cats<-cats[pops!=h1 & pops!=target]
		resdir<-paste(sep="_", Dpref, "h2", h1, target, freqprefbase, RandExt)
		cmds<-paste(Dprog, " h2=", a, " h1=", h1, " target=", target, " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep="")
		#filenames will have the names of the output files that are produced
		filenames<-paste(sep="", resdir, "/f3_", freqprefbase, "_", h1, "_", a, "_", target, ".txt")
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
	names(restab)<-c("Apop", "Bpop", "Opop", "f3", "SE", "Z", "nSNPs", "nBlocks", "Plot", "Cat")
	#resfile has the name of the concatenated results, will be used below, to enter the plotting section
	resfile<-paste(sep="", "Results_", resdir, ".txt")
	resfilex<-resfile
	write.table(restab, file=resfile, quote=F, sep="\t", col.names=T, row.names=F)

#jump to here if resfile is specified
if(!is.na(resfile)){
	message("Plotting from file... ")
	Results<-read.table(resfile, as.is=T, h=T)
	Results<-Results[Results$Plot,]
	Results$f3<-as.numeric(Results$f3)
	Results$SE<-as.numeric(Results$SE)
	Results$Z<-as.numeric(Results$Z)
	Results$nSNPs<-as.numeric(Results$nSNPs)
	Results$nBlocks<-as.numeric(Results$nBlocks)
	Results$Plot<-as.logical(Results$Plot)
	Results<-Results[rowSums(is.na(Results))==0,]
	Fixed<-unname(unlist(lapply(apply(Results[,1:3], 2, table), length))==1)
	Upper<-Results$f3+qnorm(1-(.05)/2)*Results$SE
	Lower<-Results$f3-qnorm(1-(.05)/2)*Results$SE
	Results<-cbind(Results, Upper, Lower)

	#find out if there is a fixed column
	if(identical(Fixed, c(F,T,T))){
		Fixed1<-Results[1,2]
		Fixed2<-Results[1,3]
		MainText<-paste(sep="", "f3(", Fixed2, "; X, ", Fixed1, ")")
		variablecol<-Results$Apop
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$Apop, levels=Results$Apop[sort.list(Results$f3)])
		names(Results)<-gsub("Apop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$f3)])
	}else if(identical(Fixed, c(T,F,T))){
		Fixed1<-Results[1,1]
		Fixed2<-Results[1,3]
		MainText<-paste(sep="", "f3(", Fixed2, "; ", Fixed1, ", X)")
		variablecol<-Results$Bpop
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$Bpop, levels=Results$Bpop[sort.list(Results$f3)])
		names(Results)<-gsub("Bpop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$f3)])
	}else if(identical(Fixed, c(T,T,F))){
		Fixed1<-Results[1,1]
		Fixed2<-Results[1,2]
		MainText<-paste(sep="", "f3(X; ", Fixed1, ", ", Fixed2, ")")
		variablecol<-Results$Opop
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$Opop, levels=Results$Opop[sort.list(Results$f3)])
		names(Results)<-gsub("Opop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$f3)])
	}else{
		message(USAGE)
		stop("There are less than 2 fixed columns in the results. ")
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

	#from the categories to be plotted, get them in order from CatSect
	CatsInPlot<-CatSect[CatSect$Cat %in% Results$Cat,1]
	#change Results$Cat order as a factor
	Results$Cat<-factor(Results$Cat, levels=CatsInPlot)

	F3PLOT<-ggplot(Results, aes(x=f3, y=X, color=Cat, fill=Cat, shape=Cat))+
	geom_errorbarh(aes(xmin=Lower, xmax=Upper, colour=Cat), show.legend=F, height=0)+
	geom_point(aes(shape=Cat), size=2.7, show.legend=ifelse(length(table(Results$Cat))==1, F, T))+
	theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.key=element_rect(fill="white", colour="white"))+
	xlab(MainText)+
	ylab("X")+
	scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=1))+
	scale_color_manual(values=OuterColVals)+
	scale_shape_manual(values=PchVals)

	#print the plot
	pdf(paste(sep="", resfile, ".pdf"), height=dim(Results)[1]*6.5/52, width=ifelse(length(table(Results$Cat))==1, 4, 6), useDingbats=F)
	print(F3PLOT)
	dev.off()
	message(paste(sep="", "Plot located in: ", paste(sep="", resfile, ".pdf")))
}

#############
#Do y f3

h1<-y

	#build cmds
	#keep x random string to differentiate runs

	Dpref<-"f3_"

	Dprog<-"getf3.py"
	if(rmtrns==1){
		Dprog<-"getf3NT.py"
	}

	freqprefbase<-strsplit(freqpref, "/")[[1]][length(strsplit(freqpref, "/")[[1]])]

	#build cmds based on what is variable
	#if h1 is variable, h2 and target are fixed
	cats<-oricats
	if(is.na(h2)){
		a<-pops
		a<-a[pops!=h1 & pops!=target]
		cats<-cats[pops!=h1 & pops!=target]
		resdir<-paste(sep="_", Dpref, "h2", h1, target, freqprefbase, RandExt)
		cmds<-paste(Dprog, " h2=", a, " h1=", h1, " target=", target, " freqpref=", freqpref, " resdir=", resdir, " bsize=", as.character(bsize), " category=", cats, sep="")
		#filenames will have the names of the output files that are produced
		filenames<-paste(sep="", resdir, "/f3_", freqprefbase, "_", h1, "_", a, "_", target, ".txt")
	}

	message("Running tests... ")
	system(paste0("mkdir ", resdir))
	#run cmds

	foreach(i=cmds) %dopar%{
		system(i)
	}

	#write results
	#read in all produced output files
	restab<-read.table(pipe(paste("cat", paste(filenames, collapse=" "))), as.is=T)
	#give names to data frame
	names(restab)<-c("Apop", "Bpop", "Opop", "f3", "SE", "Z", "nSNPs", "nBlocks", "Plot", "Cat")
	#resfile has the name of the concatenated results, will be used below, to enter the plotting section
	resfile<-paste(sep="", "Results_", resdir, ".txt")
	resfiley<-resfile
	write.table(restab, file=resfile, quote=F, sep="\t", col.names=T, row.names=F)

#jump to here if resfile is specified
if(!is.na(resfile)){
	message("Plotting from file... ")
	Results<-read.table(resfile, as.is=T, h=T)
	Results<-Results[Results$Plot,]
	Results$f3<-as.numeric(Results$f3)
	Results$SE<-as.numeric(Results$SE)
	Results$Z<-as.numeric(Results$Z)
	Results$nSNPs<-as.numeric(Results$nSNPs)
	Results$nBlocks<-as.numeric(Results$nBlocks)
	Results$Plot<-as.logical(Results$Plot)
	Results<-Results[rowSums(is.na(Results))==0,]
	Fixed<-unname(unlist(lapply(apply(Results[,1:3], 2, table), length))==1)
	Upper<-Results$f3+qnorm(1-(.05)/2)*Results$SE
	Lower<-Results$f3-qnorm(1-(.05)/2)*Results$SE
	Results<-cbind(Results, Upper, Lower)

	#find out if there is a fixed column
	if(identical(Fixed, c(F,T,T))){
		Fixed1<-Results[1,2]
		Fixed2<-Results[1,3]
		MainText<-paste(sep="", "f3(", Fixed2, "; X, ", Fixed1, ")")
		variablecol<-Results$Apop
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$Apop, levels=Results$Apop[sort.list(Results$f3)])
		names(Results)<-gsub("Apop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$f3)])
	}else if(identical(Fixed, c(T,F,T))){
		Fixed1<-Results[1,1]
		Fixed2<-Results[1,3]
		MainText<-paste(sep="", "f3(", Fixed2, "; ", Fixed1, ", X)")
		variablecol<-Results$Bpop
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$Bpop, levels=Results$Bpop[sort.list(Results$f3)])
		names(Results)<-gsub("Bpop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$f3)])
	}else if(identical(Fixed, c(T,T,F))){
		Fixed1<-Results[1,1]
		Fixed2<-Results[1,2]
		MainText<-paste(sep="", "f3(X; ", Fixed1, ", ", Fixed2, ")")
		variablecol<-Results$Opop
		#this was the previous way to represent the fixed column, 
		#now we just change col names and reorder Results accordingly
		#X<<-factor(Results$Opop, levels=Results$Opop[sort.list(Results$f3)])
		names(Results)<-gsub("Opop", "X", names(Results))
		Results$X<-factor(Results$X, levels=Results$X[sort.list(Results$f3)])
	}else{
		message(USAGE)
		stop("There are less than 2 fixed columns in the results. ")
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

	#from the categories to be plotted, get them in order from CatSect
	CatsInPlot<-CatSect[CatSect$Cat %in% Results$Cat,1]
	#change Results$Cat order as a factor
	Results$Cat<-factor(Results$Cat, levels=CatsInPlot)

	F3PLOT<-ggplot(Results, aes(x=f3, y=X, color=Cat, fill=Cat, shape=Cat))+
	geom_errorbarh(aes(xmin=Lower, xmax=Upper, colour=Cat), show.legend=F, height=0)+
	geom_point(aes(shape=Cat), size=2.7, show.legend=ifelse(length(table(Results$Cat))==1, F, T))+
	theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.key=element_rect(fill="white", colour="white"))+
	xlab(MainText)+
	ylab("X")+
	scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=1))+
	scale_color_manual(values=OuterColVals)+
	scale_shape_manual(values=PchVals)

	#print the plot
	pdf(paste(sep="", resfile, ".pdf"), height=dim(Results)[1]*6.5/52, width=ifelse(length(table(Results$Cat))==1, 4, 6), useDingbats=F)
	print(F3PLOT)
	dev.off()
	message(paste(sep="", "Plot located in: ", paste(sep="", resfile, ".pdf")))
}

##########
#Call PlotPairwiseOutF3.R

#Build cmd

if(!is.na(legfile)){
	plotcmd<-paste(sep="", "PlotPairwiseOutF3.R x=", x, " y=", y, " xfile=", resfilex, " yfile=", resfiley, " target=", target, " pref=", freqprefbase, "_", RandExt, " legfile=", legfile)
}else{
	plotcmd<-paste(sep="", "PlotPairwiseOutF3.R x=", x, " y=", y, " xfile=", resfilex, " yfile=", resfiley, " target=", target, " pref=", freqprefbase, "_", RandExt)
}

#Plot

system(plotcmd)