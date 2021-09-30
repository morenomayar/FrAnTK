#!/usr/bin/Rscript

#####ArgParse#####

args<-commandArgs(TRUE)

options(scipen=999)

PossibleFirstArgNames<-c("help", "BuildFreqs", "getf3", "getD", "getF4", "getF4Ratio", "getF4subtr", "getPWdist", "getEnhD", "getDstrat", "getDstrat2", "getDtrip", "autof3wfixed", "autoPWf3wfixed", "PlotPairwiseOutF3", "autoDwfixed", "autof4subtr", "autoPWdistwfixed", "autoDEnhwfixed", "autoDstrat", "addBams")

USAGE<-("FrAnTK master wrapper. 

	Usage: FrAnTK <command> [command specific arguments]

	Command specific arguments should be passed as: arg=value. Make sure that there is no spaces between arg, =, and value.

	Command-specific help is displayed when no specific arguments are supplied. 

	Run FrAnTK without arguments to get this message.

	---------
	Commands:
	--------- 

	help\tGet this message. 

	-------------------------------
	Precomputing allele frequencies
	-------------------------------

	BuildFreqs\tGet allele frequencies from a plink file for frquency-based-statistics computation. (Call BuildFreqs.py)

	----------------------------
	Computing single statistics 
	----------------------------

	getf3\tCompute f3-statistics of the form f3(h1, h2; target).  (Call getf3.py) Wrapper: autof3wfixed and autoPWf3wfixed. 

	getD\tCompute D-statistics of the form D(h1, h2; h3, h4).  (Call getD.py) Wrapper: autoDwfixed. 

	getF4\tCompute f4-statistics of the form f4(h1, h2; h3, h4).  (Call getF4.py) Wrapper: autoDwfixed. 

	getF4Ratio\tCompute f4-ratio of the form f4(h1, h4; x, h3)/f4(h1, h4; h2, h3), using definition in Patterson et al, 2012 [f4(A, O; x, C)/f4(A, O; B, C)]. We assume x is admixed between pops related to h2 (B) and h3 (C); and h1 (A) forms a clade with h2.  (Call getF4Ratio.py) 

	getF4subtr\tCompute admixture subtracted F4 statistics of the form (f4(h1, h2; h3, h4)-padm*f4(h1, x; h3, h4))/(1-padm), using definition in Reich et al. 2012.  (Call getF4subtr.py) Wrapper: autof4subtr. 

	getPWdist\tCompute pairwise distance between 2 pops, defined as the average of p_h1*q_h2+q_h1*p_h2 over all sites.  (Call getPWdist.py) Wrapper: autoPWdistwfixed.  

	getEnhD\tCompute enhanced D-statistics of the form D(h1, h2; h3, h4) (Meyer et al, 2012). We restrict to sites where h4 is invariant for one allele. Works better with many sites and not too many individuals/populations with very different ancestries in h4.  (Call getEnhD.py) Wrapper: autoDEnhwfixed. 

	getDstrat\tCompute D-statistics of the form D(h1, h2; h3, h4), on sites where there are (minder, maxder] derived alleles in h1+h2+h3+h4. Similar to what is done in Prüfer et al, 2013, but we pick segregating sites in h1+h2+h3+h4.  (Call getDstrat.py) Wrapper: autoDstrat. 

	getDstrat2\tCompute D-statistics of the form D(h1, h2; h3, h4), on sites where there are (minder, maxder] derived alleles in the complete panel (freqpref). Similar to what is done in Prüfer et al, 2013, but we pick segregating sites in the complete panel.  (Call getDstrat2.py) Wrapper: autoDstrat. 
	
	getDtrip\tCompute the three possible arrangements of D(h1, h2; h3, h4).  (Call getDtrip.py) 

	-------------------
	Automated wrappers 
	-------------------

	autof3wfixed\tCompute f3 statistics of the form f3(h1, h2; target). Fix two of h1,h2,target and loop over the remaining one. (Call autof3wfixed.R) 

	autoPWf3wfixed\tCompute f3 statistics of the form f3(h1, h2; target), for two samples and get pairwise comparisons. x and y will be fixed in h1 and h2 will be all other pops in the panel, target will also be fixed.  (Call autof3wfixed.R) 

	PlotPairwiseOutF3\tPlot pairwise comparisons of outgroup f3 statistics of the form f3(h1, h2; target) and standardized residuals from lm.  (Call PlotPairwiseOutF3.R) 

	autoDwfixed\tCompute D statistics of the form D(h1, h2; h3, h4). Fix three of h1,h2,h3,h4 and loop over the remaining one. Run with f4=1 to compute f4 instead of D.  (Call autoDwfixed.R) 

	autof4subtr\tCompute admixture subtracted f4 statistics of the form (f4(h1, h2; h3, h4)-padm*f4(h1, x; h3, h4))/(1-padm), using (allele frequency) definition in Reich et al, 2012 (Nature). We compute f4_subtr for values of padm within the range [minp, maxp], with step size pstep.  (Call autof4subtr.R) 

	autoPWdistwfixed\tCompute pairwise distance between 2 pops, defined as the average of p_h1*q_h2+q_h1*p_h2 over all sites. Fix one of h1,h2 and loop over the remaining one.  (Call autoPWdistwfixed.R) 

	autoDEnhwfixed\tCompute enhanced D-statistics of the form D(h1, h2; h3, h4), restricting to sites where h4 is invariant for one allele (Meyer et al, 2012). Works better with many sites and not too many individuals/populations with very different ancestries in h4. Fix three of h1,h2,h3,h4 and loop over the remaining one.  (Call autoDEnhwfixed.R) 

	autoDstrat\tCompute D statistics of the form D(h1, h2; h3, h4) on sites with a given number of derived alleles. We compute D for sites with [minder, maxder] (step size=dstep) derived alleles, where h1,h2,h3,h4 have mincops1,mincops2,mincops3,mincops4 copies. Set strattype=1 or strattype=2 to count derived alleles in h1,h2,h3,h4 (1) or in the complete panel (2).  (Call autoDstrat.R) 

	----------------------
	Adding sequencing data
	----------------------

	addBams\tAdd sequencing data in one or more bam files to a precomputed frequency file (output from BuildFreqs.py). For each bam file, one random allele is sampled at each site present in the frequency file. (Call addBams.py). 

	See github.com/morenomayar/FrAnTK#adding-sequencing-data-2-bam2plinkpy and github.com/morenomayar/FrAnTK#adding-sequencing-data-without-a-genotype-reference-dataset-builddummyfreqspy for additional seq. data manipulation options. 

	----------------------
	Data format conversion
	----------------------

	Freqs2Treemix\tConvert freqs file to treemix input. 

	--------
	Examples
	--------

	FrAnTK BuildFreqs plinkpref=plinkfilename clustfile=popfilename npops=120 prefout=prefix
	FrAnTK getEnhD h1=popnameh1 h2=popnameh2 h3=popnameh3 h4=popnameh4_1,popnameh4_2,popnameh4_3 bsize=3500000 freqpref=prefix resdir=dirname
	FrAnTK autoDwfixed freqpref=prefix h1=popnameh1 h2=popnameh2 nthr=50 bsize=3500000 catfile=categories legfile=legendvalues rmtrns=1 f4=0
	FrAnTK autof3wfixed resfile=txtoutputfile catfile=categories legfile=legendvalues
	FrAnTK addBams listname=list.txt freqpref= prefix newpref=prefixforout nthr=2

")

#print(args)

if(length(args)==0){
	message(USAGE)
	opt <- options(show.error.messages = FALSE)
	on.exit(options(opt))
	stop()
}

if(!(args[1] %in% PossibleFirstArgNames)){
	message(USAGE)
	stop(paste0("Unrecognized FrAnTK command: ", args[1]))
	
}

if(args[1]=="help"){
	message(USAGE)
	opt <- options(show.error.messages = FALSE)
	on.exit(options(opt))
	stop()
}

	ShortCallHash<-new.env()
	ShortCallHash[["BuildFreqs"]]<-"BuildFreqs.py"
	ShortCallHash[["getf3"]]<-"getf3.py"
	ShortCallHash[["getD"]]<-"getD.py"
	ShortCallHash[["getF4"]]<-"getF4.py"
	ShortCallHash[["getF4Ratio"]]<-"getF4Ratio.py"
	ShortCallHash[["getF4subtr"]]<-"getF4subtr.py"
	ShortCallHash[["getPWdist"]]<-"autoPWdistwfixed"
	ShortCallHash[["getEnhD"]]<-"getEnhD.py"
	ShortCallHash[["getDstrat"]]<-"getDstrat.py"
	ShortCallHash[["getDstrat2"]]<-"getDstrat2.py"
	ShortCallHash[["getDtrip"]]<-"getDtrip.py"
	ShortCallHash[["autof3wfixed"]]<-"autof3wfixed.R"
	ShortCallHash[["autoPWf3wfixed"]]<-"autof3wfixed.R"
	ShortCallHash[["PlotPairwiseOutF3"]]<-"PlotPairwiseOutF3"
	ShortCallHash[["autoDwfixed"]]<-"autoDwfixed.R"
	ShortCallHash[["autof4subtr"]]<-"autof4subtr.R"
	ShortCallHash[["autoPWdistwfixed"]]<-"autoPWdistwfixed.R"
	ShortCallHash[["autoDEnhwfixed"]]<-"autoDEnhwfixed.R"
	ShortCallHash[["autoDstrat"]]<-"autoDstrat.R"
	ShortCallHash[["addBams"]]<-"addBams.py"
	ShortCallHash[["Freqs2Treemix"]]<-"Freqs2Treemix.py"

	#get tool name
	tool<-ShortCallHash[[args[1]]]
	
	if(length(args)==1){
		system(tool)
	}

	if(length(args)>1){
		#get cmd-specific arguments
		cmdargs<-args[-1]
		cmdargs<-paste0(cmdargs, collapse=" ")
		#build frantk cmd
		frantkcmd<-paste(tool, cmdargs)
		system(frantkcmd)
	}






