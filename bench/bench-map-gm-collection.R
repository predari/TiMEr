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


source('multiplot.R')


############ LOADING DATA ############

# load raw data
ALLDATA <- read.table(input)

options(max.print = 10000)

############ PROCESSING DATA ############

# set column names
# COLLECTION GRAPH NVTXS NEDGES AVGDEG MINDEG MAXDEG METHOD NPARTS NTHREADS SEED UB EDGECut NBISLANDS TIME CHECK
# method  collection   graph   topo  seeds  duration   maxCongestion    maxDilation    avgDilation  Cut   TCV
# alma data
names(ALLDATA) <- c('FILE','COLLECTION','METHOD','TOPOLOGY','SEED','Time','MAXCON','MAXDIL','AVGDIL','Cut','TCV')

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
cat("METHODREF = ", METHODREF, "\n")

ALLDATA[c('FILE','METHOD','TOPOLOGY','Cut','Time')]


# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(COLLECTION, FILE, METHOD, TOPOLOGY),
                  transform,
                  TimeMEAN=mean(Time,na.rm=TRUE),
                  MAXCONMEAN=mean(MAXCON,na.rm=TRUE),
                  MAXDILMEAN=mean(MAXDIL,na.rm=TRUE),
                  AVGDILMEAN=mean(AVGDIL,na.rm=TRUE),
                  CutMEAN=mean(Cut,na.rm=TRUE),
                  TCVMEAN=mean(TCV,na.rm=TRUE),

		  minTime=min(Time,na.rm=TRUE),
                  MAXCONMIN=min(MAXCON,na.rm=TRUE),
                  MAXDILMIN=min(MAXDIL,na.rm=TRUE),
                  AVGDILMIN=min(AVGDIL,na.rm=TRUE),
                  minCut=min(Cut,na.rm=TRUE),
                  TCVMIN=min(TCV,na.rm=TRUE),

		  maxTime=max(Time,na.rm=TRUE),
                  MAXCONMAX=max(MAXCON,na.rm=TRUE),
                  MAXDILMAX=max(MAXDIL,na.rm=TRUE),
                  AVGDILMAX=max(AVGDIL,na.rm=TRUE),
                  maxCut=max(Cut,na.rm=TRUE),
                  TCVMAX=max(TCV,na.rm=TRUE)
		  )

# here I have the arithmetic mean, min and max values over the seeds for each graphxtopo for both methods
NORMDATA[c('FILE','METHOD','TOPOLOGY','minCut','CutMEAN','maxCut','minTime','TimeMEAN','maxTime')]


# normalize just the mean value based on the reference mean value
NORMDATA <- ddply(NORMDATA, .(COLLECTION, FILE, TOPOLOGY),
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

		  TCVR=TCVMEAN/TCVMEAN[METHOD==METHODREF],
    		  TCVMINR=TCVMIN/TCVMIN[METHOD==METHODREF],
		  TCVMAXR=TCVMAX/TCVMAX[METHOD==METHODREF],

		  TimeR=TimeMEAN/TimeMEAN[METHOD==METHODREF],
		  minTimeR=minTime/minTime[METHOD==METHODREF],
		  maxTimeR=maxTime/maxTime[METHOD==METHODREF]
		  )

# after that ref method should be on 1 for all metrics
NORMDATA[c('FILE','COLLECTION','METHOD','TOPOLOGY', 'minCutR','CutR','maxCutR','TimeR')] ##'minCutR','CutR','maxCutR','TimeR')]

# remove useless columns
NORMDATA[c("Cut","Time","CutMEAN","TimeMEAN")] <- list(NULL) 

# compute geom mean
MEANDATA <- ddply(NORMDATA,
                  .(COLLECTION, METHOD, TOPOLOGY), # compute separate mean for each METHOD x NPARTS
                  summarise,
                  Cut=exp(mean(log(CutR))),
		  minCut = exp(mean(log(minCutR))),
		  maxCut = exp(mean(log(maxCutR))),
		  CutSD=exp(sd(log(CutR))),
		  minCutSD=exp(sd(log(minCutR))),
		  maxCutSD=exp(sd(log(maxCutR))),

		  Time=exp(mean(log(TimeR))),
                  minTime=exp(mean(log(minTimeR))),
		  maxTime=exp(mean(log(maxTimeR))),
		  TimeSD=exp(sd(log(TimeR))),
		  minTimeSD=exp(sd(log(minTimeR))),
		  maxTimeSD=exp(sd(log(maxTimeR))),

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
		  
                  TCV=exp(mean(log(TCVR))),
    		  TCVMIN=exp(mean(log(TCVMINR))),
		  TCVMAX=exp(mean(log(TCVMAXR))),
                  TCVSD=exp(sd(log(TCVR))),
		  TCVMINSD=exp(sd(log(TCVMINR))),
		  TCVMAXSD=exp(sd(log(TCVMAXR)))

                  )




############ REMOVE A METHOD ############

## MEANDATA <- MEANDATA[MEANDATA$METHOD != 'SCOTCH', ]

############ DUMP RESULTS ############

# dump meandata (geom mean) before melting
MEANDATA[c('COLLECTION', 'METHOD','TOPOLOGY', 'minCut','Cut','maxCut','CutSD', 'Time','TimeSD')] ## 'minCut','Cut','maxCut','CutSD', 'Time','TimeSD')] 

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
    ## facet_grid(. ~ COLLECTION, scales="free")+
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") + 
    xlab("topologies") + ylab("relative cut") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  
# to remove diagonal lines in legend.key: + guides(fill = guide_legend(override.aes = list(colour = NULL)))


   durplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=Time, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=Time/TimeSD, ymax=Time*TimeSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") +
    ## facet_grid(. ~ COLLECTION, scales="free") +
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") + 
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
            ## facet_grid(. ~ COLLECTION, scales="free") +
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") + 
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
            ## facet_grid(. ~ COLLECTION, scales="free") +
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") + 
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
            ## facet_grid(. ~ COLLECTION, scales="free") +
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") + 
    xlab("topologies") + ylab("relative avgdil") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  

   tcvplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=TCV, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=TCV/TCVSD, ymax=TCV*TCVSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") +
            ## facet_grid(. ~ COLLECTION, scales="free") +
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") + 
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
MEANDATA[c("METHOD")] <- list(NULL)




N_MEANDATA <- MEANDATA[MEANDATA$COLLECTION == 'network', ]
D_MEANDATA <- MEANDATA[MEANDATA$COLLECTION == 'dimacs', ]

 durdata <- N_MEANDATA[c('TOPOLOGY','minTime','maxTime','Time')]
 durdata <- melt(durdata, id.vars = "TOPOLOGY")
 durdata
 durdatasd <- N_MEANDATA[c('TOPOLOGY','minTimeSD','maxTimeSD','TimeSD')]
 durdatasd <- melt(durdatasd, id.vars = "TOPOLOGY")
 durdatasd

 # duration plot for mean max min with sd for each value
 n_durdata_plot <- ggplot(durdata, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/ durdatasd$value, ymax=value* durdatasd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results - network collection") +
    #coord_cartesian(ylim = c(0.0, 5)) +	  
    # for scotch
     coord_cartesian(ylim = c(0.8, 100)) +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    #scale_y_continuous(breaks=c(0.5,1.00,1.50,2.00))
    # for scotch
    scale_y_continuous(breaks=c(1.00,5.00,10.00,20.00,50.00))	


 durdata <- D_MEANDATA[c('TOPOLOGY','minTime','maxTime','Time')]
 durdata <- melt(durdata, id.vars = "TOPOLOGY")
 durdata
 durdatasd <- D_MEANDATA[c('TOPOLOGY','minTimeSD','maxTimeSD','TimeSD')]
 durdatasd <- melt(durdatasd, id.vars = "TOPOLOGY")
 durdatasd

 # duration plot for mean max min with sd for each value
 d_durdata_plot <- ggplot(durdata, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/ durdatasd$value, ymax=value* durdatasd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results  - dimacs collection") +
    #coord_cartesian(ylim = c(0.0, 5)) +	  
    # for scotch
     coord_cartesian(ylim = c(0.8, 100)) +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    #scale_y_continuous(breaks=c(0.5,1.00,1.50,2.00))
    # for scotch
    scale_y_continuous(breaks=c(1.00,5.00,10.00,20.00,50.00))	




 fdfm <- N_MEANDATA[c('TOPOLOGY','minCut','maxCut','Cut','TCVMIN','TCVMAX','TCV')]
 fdfm <- melt(fdfm, id.vars = "TOPOLOGY")

 fdfsd <- N_MEANDATA[c('TOPOLOGY','CutSD','minCutSD','maxCutSD','TCVSD','TCVMINSD','TCVMAXSD')]   
 fdfsd <- melt(fdfsd, id.vars = "TOPOLOGY")
 
n_fplot <- ggplot(fdfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/fdfsd$value, ymax=value* fdfsd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results - network collection") +
    coord_cartesian(ylim = c(0.5, 1.3)) +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    scale_y_continuous(breaks=c(0.7,0.8,1.10,1.2))

 fdfm <- D_MEANDATA[c('TOPOLOGY','minCut','maxCut','Cut','TCVMIN','TCVMAX','TCV')]
 fdfm <- melt(fdfm, id.vars = "TOPOLOGY")
 fdfm

 fdfsd <- D_MEANDATA[c('TOPOLOGY','CutSD','minCutSD','maxCutSD','TCVSD','TCVMINSD','TCVMAXSD')]   
 fdfsd <- melt(fdfsd, id.vars = "TOPOLOGY")
 fdfm
 
d_fplot <- ggplot(fdfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/fdfsd$value, ymax=value* fdfsd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results - dimacs collection") +
    coord_cartesian(ylim = c(0.5, 1.3)) +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    scale_y_continuous(breaks=c(0.7,0.8,1.10,1.2))










 dfm <- D_MEANDATA[c('TOPOLOGY','Cut','TCV')]
 dfm <- melt(dfm, id.vars = "TOPOLOGY")

 sdfm <- D_MEANDATA[c('TOPOLOGY','CutSD','TCVSD')]   
 sdfm <- melt(sdfm, id.vars = "TOPOLOGY")
 cat("after melting dfm:\n")
 dfm
 cat("sd data after melting sdfm:\n")
 sdfm

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

# this figure is good for the sd of the geometric mean over all graph
d_bplot <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    #geom_point(position=position_dodge(width=0.9), size = 2) +
    #geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results - dimacs collection") +
    coord_cartesian(ylim = c(0.5, 1.3)) +

#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
   scale_y_continuous(breaks=c(0.7,0.8,1.10,1.2))


 dfm <- N_MEANDATA[c('TOPOLOGY','Cut','TCV')]
 dfm <- melt(dfm, id.vars = "TOPOLOGY")

 sdfm <- N_MEANDATA[c('TOPOLOGY','CutSD','TCVSD')]   
 sdfm <- melt(sdfm, id.vars = "TOPOLOGY")
 cat("after melting dfm:\n")
 dfm
 cat("sd data after melting sdfm:\n")
 sdfm



# this figure is good for the sd of the geometric mean over all graph
n_bplot <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    #geom_point(position=position_dodge(width=0.9), size = 2) +
    #geom_errorbar(aes(ymin=value/sdfm$value, ymax=value*sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results - network collection") +
    coord_cartesian(ylim = c(0.5, 1.3)) +

#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
   scale_y_continuous(breaks=c(0.7,0.8,1.10,1.2))



pdf(output)
 cutplot
 durplot
 n_durdata_plot
 d_durdata_plot
# maxconplot
# maxdilplot
# avgdilplot
 tcvplot
# sd for the geometric mean value (bar representation)
n_bplot
d_bplot
n_fplot
d_fplot

cat("writing plot in ", output, "\n");

### EOF
