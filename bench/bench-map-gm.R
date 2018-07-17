#!/usr/bin/Rscript --vanilla

# if you want see commands in output file
# options(echo=TRUE)


############ INPUT ARGS ############

args <- commandArgs(TRUE)
# print(args)

input <- args[1]
output <- args[2]

if(length(args) != 2) { print("usage: input.data output.pdf"); q(); }

############ AUXILIARY ROUTINES ############

library(plyr)      # sudo apt-get install r-cran-plyr
library(ggplot2)   # sudo apt-get install r-cran-ggplot2
library(xtable) # install via R environment with install.packages("xtable"). latex in R library
require(tikzDevice)
library(systemfit)
require(graphics)

source('multiplot.R')


############ LOADING DATA ############

# load raw data
ALLDATA <- read.table(input)

options(max.print = 10000)

############ PROCESSING DATA ############

# set column names
# COLLECTION GRAPH NVTXS NEDGES AVGDEG MINDEG MAXDEG METHOD NPARTS NTHREADS SEED UB EDGECut NBISLANDS TIME CHECK
# method  collection   graph   topo  seeds  duration   maxCongestion    maxDilation    avgDilation  Cut   Co
# alma data
names(ALLDATA) <- c('FILE','COLLECTION','METHOD','TOPOLOGY','SEED','T','MAXCON','MAXDIL','AVGDIL','Cut','Co')

# remove COLLECTION column
ALLDATA$FILE <- paste(ALLDATA$COLLECTION, ALLDATA$FILE, sep="_")
ALLDATA$COLLECTION <- list(NULL) 

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
#METHODREF <- "TiMEr"
#METHODREF <- "SCOTCH"
#METHODREF <- "KAHIP"
cat("METHODREF = ", METHODREF, "\n")

ALLDATA[c('FILE','METHOD','TOPOLOGY','MAXCON','T')]


# normalize data: TIME and EDGECut relative to the reference method average
# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(FILE, METHOD, TOPOLOGY),
                  transform,
                  TMEAN=mean(T,na.rm=TRUE),
                  MAXCONMEAN=mean(MAXCON,na.rm=TRUE),
                  MAXDILMEAN=mean(MAXDIL,na.rm=TRUE),
                  AVGDILMEAN=mean(AVGDIL,na.rm=TRUE),
                  CutMEAN=mean(Cut,na.rm=TRUE),
                  CoMEAN=mean(Co,na.rm=TRUE),

		  minT=min(T,na.rm=TRUE),
                  MAXCONMIN=min(MAXCON,na.rm=TRUE),
                  MAXDILMIN=min(MAXDIL,na.rm=TRUE),
                  AVGDILMIN=min(AVGDIL,na.rm=TRUE),
                  minCut=min(Cut,na.rm=TRUE),
                  minCo=min(Co,na.rm=TRUE),

		  maxT=max(T,na.rm=TRUE),
                  MAXCONMAX=max(MAXCON,na.rm=TRUE),
                  MAXDILMAX=max(MAXDIL,na.rm=TRUE),
                  AVGDILMAX=max(AVGDIL,na.rm=TRUE),
                  maxCut=max(Cut,na.rm=TRUE),
                  maxCo=max(Co,na.rm=TRUE)
		  )

# here I have the arithmetic mean, min and max values over the seeds for each graphxtopo for both methods
NORMDATA[c('FILE','METHOD','TOPOLOGY','minCut','CutMEAN','maxCut','minT','TMEAN','maxT')]


# normalize just the mean value based on the reference mean value
NORMDATA <- ddply(NORMDATA, .(FILE, TOPOLOGY),
                  transform,
                  CutR=CutMEAN/CutMEAN[METHOD==METHODREF],
	          minCutR=minCut/minCut[METHOD==METHODREF],
		  maxCutR=maxCut/maxCut[METHOD==METHODREF],	

		  MAXCONR=MAXCONMEAN/MAXCONMEAN[METHOD==METHODREF],
		  MAXCONMINR=MAXCONMIN/MAXCONMIN[METHOD==METHODREF],
		  MAXCONMAXR=MAXCONMAX/MAXCONMAX[METHOD==METHODREF],	

		  MAXDILR=MAXDILMEAN/MAXDILMEAN[METHOD==METHODREF],
		  MAXDILMINR=MAXDILMIN/MAXDILMIN[METHOD==METHODREF],
		  MAXDILMAXR=MAXDILMAX/MAXDILMAX[METHOD==METHODREF],

		  AVGDILR=AVGDILMEAN/AVGDILMEAN[METHOD==METHODREF],
  		  AVGDILMINR=AVGDILMIN/AVGDILMIN[METHOD==METHODREF],
		  AVGDILMAXR=AVGDILMAX/AVGDILMAX[METHOD==METHODREF],

		  CoR=CoMEAN/CoMEAN[METHOD==METHODREF],
    		  minCoR=minCo/minCo[METHOD==METHODREF],
		  maxCoR=maxCo/maxCo[METHOD==METHODREF],

		  TR=TMEAN/TMEAN[METHOD==METHODREF],
		  minTR=minT/minT[METHOD==METHODREF],
		  maxTR=maxT/maxT[METHOD==METHODREF]
		  )

# after that ref method should be on 1 for all metrics
NORMDATA[c('FILE','METHOD','TOPOLOGY', 'minCutR','CutR','maxCutR','TR')] ##'minCutR','CutR','maxCutR','TR')]

# remove useless columns
NORMDATA[c("Cut","T","CutMEAN","TMEAN")] <- list(NULL) 

# compute geom mean
MEANDATA <- ddply(NORMDATA,
                  .(METHOD, TOPOLOGY), # compute separate mean for each METHOD x NPARTS
                  summarise,
                  Cut=exp(mean(log(CutR))),
		  minCut = exp(mean(log(minCutR))),
		  maxCut = exp(mean(log(maxCutR))),
		  CutSD=exp(sd(log(CutR))),
		  minCutSD=exp(sd(log(minCutR))),
		  maxCutSD=exp(sd(log(maxCutR))),

		  T=exp(mean(log(TR))),
                  minT=exp(mean(log(minTR))),
		  maxT=exp(mean(log(maxTR))),
		  TSD=exp(sd(log(TR))),
		  minTSD=exp(sd(log(minTR))),
		  maxTSD=exp(sd(log(maxTR))),

		  MAXCON=exp(mean(log(MAXCONR))),
		  MAXCONMIN=exp(mean(log(MAXCONMINR))),
		  MAXCONMAX=exp(mean(log(MAXCONMAXR))),
		  MAXCONSD=exp(sd(log(MAXCONR))),

		  MAXDIL=exp(mean(log(MAXDILR))),
		  MAXDILMIN=exp(mean(log(MAXDILMINR))),
		  MAXDILMAX=exp(mean(log(MAXDILMAXR))),
		  MAXDILSD=exp(sd(log(MAXDILR))),
		  
		  AVGDIL=exp(mean(log(AVGDILR))),
  		  AVGDILMIN=exp(mean(log(AVGDILMINR))),
		  AVGDILMAX=exp(mean(log(AVGDILMAXR))),
		  AVGDILSD=exp(sd(log(AVGDILR))),
		  
                  Co=exp(mean(log(CoR))),
    		  minCo=exp(mean(log(minCoR))),
		  maxCo=exp(mean(log(maxCoR))),
                  CoSD=exp(sd(log(CoR))),
		  minCoSD=exp(sd(log(minCoR))),
		  maxCoSD=exp(sd(log(maxCoR)))

                  )


# add TOPOLOGY column to use TOPO as a factor rather than numeric
# NORMDATA$NPARTSF <- factor(NORMDATA$NPARTS) 
# MEANDATA$NPARTSF <- factor(MEANDATA$NPARTS) 


############ REMOVE A METHOD ############

## MEANDATA <- MEANDATA[MEANDATA$METHOD != 'SCOTCH', ]

############ DUMP RESULTS ############




############ MY OWN COLORS ############

### configure my own palette based on brewer one

library(RColorBrewer)
PALSIZE <- max(length(unique(ALLDATA$METHOD)),3)
MYPAL <- brewer.pal(PALSIZE, "Blues")
MYPAL[1] <- "blue"
MYPAL[2] <- "lightblue"
MYPAL[3] <- "green"
MYPAL[4] <- "lightgreen"
MYPAL[5] <- "purple"
MYPAL[6] <- "pink"
MYPAL[7] <- "red"
MYPAL[8] <- "magenta"
MYPAL[9] <- "yellow"

############ BOXPLOT ############

 library(RColorBrewer)
 PALSIZE <- max(length(unique(ALLDATA$METHOD)),3)
 MYPAL <- brewer.pal(PALSIZE, "Set1")
 MYPAL  <- brewer.pal(PALSIZE, "Spectral")
d_colors_methods <- c(
    'SCOTCH'      = "orange2", #"coral2", #  "#E47658" "#F0936D" "#F4A989" "darkblue",
    'RB*'    = "palegreen4", #"yellow",
    'KGGGP_L'    = "magenta",
    'KGGGP'   = "cadetblue2", # red
    'KMETIS' = "#F8BFA5" #"green"
)

############ LINEPLOT ############

############ BARPLOT ############

# display graphics with lines and errorbar
   cutplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=Cut, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=Cut/CutSD, ymax=Cut*CutSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative cut") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  
# to remove diagonal lines in legend.key: + guides(fill = guide_legend(override.aes = list(colour = NULL)))


   durplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=T, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=T/TSD, ymax=T*TSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative dur") +
    coord_cartesian(ylim = c(0.0, 3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
    scale_y_continuous(breaks=c(0.50,1.00,1.50,5.00,50.00))



   maxconplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=MAXCON, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=MAXCON/MAXCONSD, ymax=MAXCON*MAXCONSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative maxcon") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  

   maxdilplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=MAXDIL, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=MAXDIL/MAXDILSD, ymax=MAXDIL*MAXDILSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative maxdil") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  

   avgdilplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=AVGDIL, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=AVGDIL/AVGDILSD, ymax=AVGDIL*AVGDILSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative avgdil") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  

   tcvplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=Co, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=Co/CoSD, ymax=Co*CoSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative tcv") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  

# reshape data to have all metrics in one graph
 library(reshape2)
# remove ref method
MEANDATA <- MEANDATA[MEANDATA$METHOD != METHODREF, ]
#MEANDATA <- MEANDATA[MEANDATA$METHOD != 'SCOTCH', ]

###### important dumping 
# dump meandata (geom mean) before melting
MEANDATA[c('METHOD','TOPOLOGY','T','Co','Cut')] ## 'minCut','Cut','maxCut','CutSD', 'T','TSD')]

cat("mean over topologies for co: ")
mean(MEANDATA[['Co']])

cat("mean over topologies for cut: ")
mean(MEANDATA[['Cut']])

cat("mean over topologies for time: ")
mean(MEANDATA[['T']])


MEANDATA[c("METHOD")] <- list(NULL)
sdfm <- MEANDATA[c('TOPOLOGY','CutSD','TSD','MAXCONSD','MAXDILSD','AVGDILSD','CoSD')]
## to exclude the following metrics
#sdfm[c("AVGDILSD")] <- list(NULL)
#sdfm[c("TSD")] <- list(NULL)


 durdata <- MEANDATA[c('TOPOLOGY','minT','T','maxT')]
 durdata
 durdata <- melt(durdata, id.vars = "TOPOLOGY")
 durdatasd <- MEANDATA[c('TOPOLOGY','minTSD','TSD','maxTSD')]
# durdatasd
 durdatasd <- melt(durdatasd, id.vars = "TOPOLOGY")
 

 # duration plot for mean max min with sd for each value
 durdata_plot <- ggplot(durdata, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/ durdatasd$value, ymax=value* durdatasd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative runtime") +
    coord_cartesian(ylim = c(0.0, 3)) +	  
    # for scotch
    # coord_cartesian(ylim = c(0.8, 100)) +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    scale_y_continuous(breaks=c(0.5,1.00,1.50,2.00))
    # for scotch
    # scale_y_continuous(breaks=c(1.00,5.00,10.00,20.00,50.00))	


mindfm <- MEANDATA[c('TOPOLOGY','minCut','minT','MAXCONMIN','MAXDILMIN','AVGDILMIN','minCo')]

maxdfm <- MEANDATA[c('TOPOLOGY','maxCut','maxT','MAXCONMAX','MAXDILMAX','AVGDILMAX','maxCo')]


  printdata <- MEANDATA[c('TOPOLOGY','minCo','Co','maxCo')]  # ,'MAXCON','MAXDIL'
  printdata

 fdfm <- MEANDATA[c('TOPOLOGY','minCut','maxCut','Cut','minCo','maxCo','Co')]  # ,'MAXCON','MAXDIL'
# for presentation: only meanCut and meanCo
# fdfm <- MEANDATA[c('TOPOLOGY','Cut','Co')] # only meanCut and meanCo  

# following two lines to change the names of columns in the fig
  colnames(fdfm) <- c('TOPOLOGY', "qCut^{gm}_{min}", "qCut^{gm}_{max}", "qCut^{gm}_{mean}","qCo^{gm}_{min}", "qCo^{gm}_{max}", "qCo^{gm}_{mean}")
  fdfm <- xtable(fdfm)

fdfm <- melt(fdfm, id.vars = "TOPOLOGY")

 fdfm


 fdfsd <- MEANDATA[c('TOPOLOGY','minCutSD','maxCutSD','CutSD','minCoSD','maxCoSD','CoSD')]  # ,'MAXCONSD','MAXDILSD'
 # for presentation: only meanCut and meanCo						    
#  fdfsd <- MEANDATA[c('TOPOLOGY','CutSD','CoSD')]					    
 fdfsd <- melt(fdfsd, id.vars = "TOPOLOGY")

 fdfsd


fplot <- ggplot(fdfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/ fdfsd$value, ymax=value* fdfsd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.5, 1.3)) +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    scale_y_continuous(breaks=c(0.6,0.8,1.00,1.2))


# remove all sd values
MEANDATA[c("CutSD","TSD","MAXCONSD","MAXDILSD","AVGDILSD","CoSD")] <- list(NULL)

# remove the extra duration sd
MEANDATA[c("maxTSD","minTSD")] <- list(NULL)

# remove the extra duration sd
MEANDATA[c("maxCutSD","minCutSD","maxCoSD","minCoSD")] <- list(NULL)


## do not include avgdil
#MEANDATA[c("AVGDIL")] <- list(NULL)
## do not include duration
#MEANDATA[c("T")] <- list(NULL)
## remove all min and max from meandata
MEANDATA[c('minCut','minT','MAXCONMIN','AVGDILMIN','MAXDILMIN','minCo','maxCut','maxT','MAXCONMAX','AVGDILMAX','MAXDILMAX','maxCo')] <- list(NULL)

# melt all the remaining metrics (excluding min max sd)
dfm <- melt(MEANDATA, id.vars = "TOPOLOGY")
# melt separately sd
sdfm <- melt(sdfm, id.vars = "TOPOLOGY")

# melt separately min max
mindfm <- melt(mindfm, id.vars = "TOPOLOGY")
maxdfm <- melt(maxdfm, id.vars = "TOPOLOGY")


# remove from the plots duration and avgdil
dfm <- dfm[dfm$variable != 'T', ]
sdfm <- sdfm[sdfm$variable != 'TSD', ]
mindfm <- mindfm[mindfm$variable != 'minT', ]
maxdfm <- maxdfm[maxdfm$variable != 'maxT', ]
# make sure all dfm sdfm mindfm maxdfm have the same size!

dfm <- dfm[dfm$variable != 'AVGDIL', ]
sdfm <- sdfm[sdfm$variable != 'AVGDILSD', ]
mindfm <- mindfm[mindfm$variable != 'AVGDILMIN', ]
maxdfm <- maxdfm[maxdfm$variable != 'AVGDILMAX', ]

dfm <- dfm[dfm$variable != 'MAXDIL', ]
sdfm <- sdfm[sdfm$variable != 'MAXDILSD', ]
mindfm <- mindfm[mindfm$variable != 'MAXDILMIN', ]
maxdfm <- maxdfm[maxdfm$variable != 'MAXDILMAX', ]

dfm <- dfm[dfm$variable != 'MAXCON', ]
sdfm <- sdfm[sdfm$variable != 'MAXCONSD', ]
mindfm <- mindfm[mindfm$variable != 'MAXCONMIN', ]
maxdfm <- maxdfm[maxdfm$variable != 'MAXCONMAX', ]

# dump
# cat("after melting dfm:\n")
# dfm
# cat("sd data after melting sdfm:\n")
# sdfm




min = dfm$value/sdfm$value
min <- melt(min)
#min

max = dfm$value*sdfm$value
max <- melt(max)
#max



PALSIZE <- max(length(unique(dfm$variable)),3)
MYPAL <- brewer.pal(PALSIZE, "Blues")
MYPAL[1] <- "purple" # "blue"
MYPAL[2] <- "lightblue"
MYPAL[3] <- "magenta" #"green"
MYPAL[4] <- "lightgreen"
MYPAL[5] <- "yellow"
MYPAL[6] <- "pink"
MYPAL[7] <- "red"
MYPAL[8] <- "green"
MYPAL[9] <- "yellow"








# this is the figure to compare min, mean and max relative values all (one by one) against the min mean and max of the reference method.
# so min mean and max are not compared among them but against 1
# bar plots are not the best for this representation!!! Try symbols
bplot_minmax <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_errorbar(aes(ymin=mindfm$value, ymax=maxdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    #geom_point(position=position_dodge(width=0.9), size = 2) +
    #geom_errorbar(aes(ymin=mindfm$value, ymax=maxdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +		       
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    scale_y_continuous(breaks=c(0.8,0.95,1.05,1.10,1.2))


# splot_minmax for symbols, not bars (same as the above though)
splot_minmax <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    #geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    #geom_errorbar(aes(ymin=mindfm$value, ymax=maxdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    geom_point(position=position_dodge(width=0.9), size = 2) +
    geom_errorbar(aes(ymin=mindfm$value, ymax=maxdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + 
   scale_y_continuous(breaks=c(0.8,0.95,1.05,1.10,1.2))		     


## for the plots with geometric sd do not include time because the results are not redable (large difference in time from one graph to another)
#dfm <- dfm[dfm$variable != 'T', ]
#sdfm <- dfm[sdfm$variable != 'TSD', ]

splot <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    #geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    #geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    geom_point(position=position_dodge(width=0.9), size = 2) +
    geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.2, 1.5)) +
#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
#+        		     scale_y_continuous(breaks=c(0.80,1.00,1.20))  

# this figure is good for the sd of the geometric mean over all graph
bplot <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    #geom_point(position=position_dodge(width=0.9), size = 2) +
    #geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.5, 1.3)) +

#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
   scale_y_continuous(breaks=c(0.7,0.8,1.10,1.2))



pdf(output)
### for now cutplot
### for now durplot
### for now durdata_plot
# maxconplot
# maxdilplot
# avgdilplot
### for now tcvplot
# sd for the geometric mean value (bar representation)
### for now bplot
#splot_minmax
fplot
cat("writing plot in ", output, "\n");

## export to eps for the fonts!!
#setEPS()
## almaIdentity-quality.eps
## almaLibtopo-quality.eps
## almaScotch-quality.eps
## almaGreedy-quality.eps
#postscript("almaGreedy-quality.eps")
#fplot
#dev.off()

### EOF
