
#####ArgParse#####

args<-commandArgs(TRUE)

#names of recognized arguments
PossibleArgNames<-c("xfile", "yfile", "x", "y", "target", "legfile", "pref", "minf", "maxf")

USAGE<-("
	Plot pairwise comparisons of outgroup f3 statistics of the form f3(h1, h2; target) and standardized residuals from lm.  
	Arguments should be passed as: arg=value. Make sure that there is no spaces between arg, =, and value.
	pref\tSTR\tPrefix for outfiles. Default: Results
	xfile\tSTR\tFile name from autof3wfixed.R. Required.
	yfile\tSTR\tFile name from autof3wfixed.R. Required.
	x\tSTR\tName of population to set in x-axis.
	y\tSTR\tName of population to set in y-axis.
	minf\tMinumum f3 to consider in pairwise comparisons. Should be used together with maxf. Default=NA.
	maxf\tMaximum f3 to consider in pairwise comparisons. Should be used together with minf. Default=NA.
	target\tSTR\tName of population that was used for f3.
	legfile\tSTR\tTab separated file specifying color and pch for each category (used for plotting and legend). Random colors and pch are chosen for categories not present in the file. Format: category outercolor innercolor pch. (Use R color names and pch numbers) .
	Sample call: PlotPairwiseOutF3.R pref=Panel x=popnameh1 y=popnameh2 xfile=Results_f3__h2_x_target_Panel_PhDkwg.txt yfile=Results_f3__h2_y_target_Panel_geKEGr.txt legfile=legendvalues 
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
ArgHash[["pref"]]<-"Results"

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
x<-ArgHash[["x"]]
y<-ArgHash[["y"]]
xfile<-ArgHash[["xfile"]]
yfile<-ArgHash[["yfile"]]
target<-ArgHash[["target"]]
legfile<-ArgHash[["legfile"]]
pref<-ArgHash[["pref"]]

flt<-F
if(!is.na(ArgHash[["minf"]]) || !is.na(ArgHash[["maxf"]])){
	if(!is.na(ArgHash[["minf"]]) && !is.na(ArgHash[["maxf"]])){
		minf<-as.numeric(ArgHash[["minf"]])
		maxf<-as.numeric(ArgHash[["maxf"]])
		flt<-T
	}else{
		message(USAGE)
		stop("minf and maxf should be used together.")
	}
}

#############

#Read files

a<-read.table(xfile, as.is=T, h=T)
a$f3<-as.numeric(a$f3)
a$SE<-as.numeric(a$SE)
a$Z<-as.numeric(a$Z)
a$nSNPs<-as.numeric(a$nSNPs)
a$nBlocks<-as.numeric(a$nBlocks)
a$Plot<-as.logical(a$Plot)
a<-a[rowSums(is.na(a))==0,]
a<-a[a$Plot==T,]
Upper<-a$f3+qnorm(1-(.05)/2)*a$SE
Lower<-a$f3-qnorm(1-(.05)/2)*a$SE
a<-cbind(a, Upper, Lower)
xl<-paste(sep="", "f3(", target, "; ", x, ", X)")
if(flt==T){
	a<-a[a$f3>=minf & a$f3<=maxf, ]
}


b<-read.table(yfile, as.is=T, h=T)
b$f3<-as.numeric(b$f3)
b$SE<-as.numeric(b$SE)
b$Z<-as.numeric(b$Z)
b$nSNPs<-as.numeric(b$nSNPs)
b$nBlocks<-as.numeric(b$nBlocks)
b$Plot<-as.logical(b$Plot)
b<-b[rowSums(is.na(b))==0,]
b<-b[b$Plot==T,]
Upper<-b$f3+qnorm(1-(.05)/2)*b$SE
Lower<-b$f3-qnorm(1-(.05)/2)*b$SE
b<-cbind(b, Upper, Lower)
yl<-paste(sep="", "f3(", target, "; ", y, ", X)")
if(flt==T){
	b<-b[b$f3>=minf & b$f3<=maxf, ]
}


#match b with a
ToRm<-NULL
for(i in a$Bpop){
	if(sum(b$Bpop==i)!=1){
		ToRm<-c(ToRm, i)
	}
}
ToRm<-unique(ToRm)
a<-a[!a$Bpop %in% ToRm, ]

bf3<-NULL
bU<-NULL
bL<-NULL
for(i in a$Bpop){
	bf3<-c(bf3, b[b$Bpop==i,4])
	bU <-c(bU, b[b$Bpop==i,11])
	bL <-c(bL, b[b$Bpop==i,12])
}

a<-cbind(a, bf3, bU, bL)

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
	for(i in unique(a$Cat)){
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

#get regression and residuals
reg<-(lm(a$f3~a$bf3))
stdres<-(reg$residuals-mean(reg$residuals))/sd(reg$residuals)

X<-factor(a$Bpop, levels=a$Bpop[sort.list(stdres)])
a<-cbind(a, stdres, X)

#Do the plots
library("ggplot2")

	#from the categories to be plotted, get them in order from CatSect
	CatsInPlot<-CatSect[CatSect$Cat %in% a$Cat,1]
	#change Results$Cat order as a factor
	a$Cat<-factor(a$Cat, levels=CatsInPlot)

F3PLOT<-ggplot(a, aes(x=f3, y=bf3))+
geom_abline(intercept=0, slope=1, color="black", linetype=2, alpha=.25)+
#geom_smooth(aes(x=f3, y=bf3), method="lm", se=F, fullrange=TRUE, show.legend=F, alpha=.35, size=.35)+
geom_smooth(method="lm", se=F, show.legend=F, fullrange=TRUE, alpha=.35, size=.35, color="black")+
geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=Cat), show.legend=F, height=0)+
geom_errorbar(aes(ymin=bU, ymax=bL, color=Cat), show.legend=F, width=0)+
geom_point(aes(color=Cat, fill=Cat, shape=Cat), size=2.3, show.legend=ifelse(length(table(a$Cat))==1, F, T))+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 0, hjust = .5, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.key=element_rect(fill="white", colour="white"))+
xlab(xl)+
ylab(yl)+
scale_fill_manual(values= InnerColVals)+
scale_color_manual(values=OuterColVals)+
scale_shape_manual(values=PchVals)

firstthr<-c(qnorm(1-(c(.05, .01, .001))/2), -qnorm(1-(c(.05, .01, .001))/2))
resPLOT<-ggplot(a, aes(x=stdres, y=X, color=Cat, fill=Cat, shape=Cat))+
geom_vline(xintercept=firstthr, colour="gray", linetype = "longdash")+
geom_point(aes(shape=Cat), size=2.7, show.legend=ifelse(length(table(a$Cat))==1, F, T))+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 0, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.key=element_rect(fill="white", colour="white"))+
xlab(paste("Standardized residuals\n", y, " <      > ", x))+
scale_fill_manual(values= InnerColVals, guide=guide_legend(ncol=1))+
scale_color_manual(values=OuterColVals)+
scale_shape_manual(values=PchVals)


pdf(paste(sep="", x, "_", y, "_", pref, "_PW_f3.pdf"), height=4.5, width=8.5, useDingbats=F)
print(F3PLOT)
dev.off()
message(paste(sep="", "Pairwise plot located in: ", paste(sep="", x, "_", y, "_", pref, "_PW_f3.pdf")))

pdf(paste(sep="", x, "_", y, "_", pref, "_res.pdf"), height=ifelse(dim(a)[1]<=30, dim(a)[1]*8/52, dim(a)[1]*6.5/52), width=ifelse(length(table(a$Cat))==1, 4, 6), useDingbats=F)
print(resPLOT)
dev.off()
message(paste(sep="", "Residuals plot located in: ", paste(sep="", x, "_", y, "_", pref, "_res.pdf")))



