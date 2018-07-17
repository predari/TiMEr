#!/usr/bin/Rscript --vanilla
# if you want see commands in output file
# options(echo=TRUE)

############ INPUT ARGS ############

args <- commandArgs(TRUE)
# print(args)

input <- args[1]
outputplot <- args[2]

if(length(args) != 2) { print("usage: input.data output.pdf"); q(); }

############ AUXILIARY ROUTINES ############

library(plyr)      # sudo apt-get install r-cran-plyr
library(ggplot2)   # sudo apt-get install r-cran-ggplot2

source('multiplot.R')

############ LOADING DATA ############

# load raw data
ALLDATA <- read.table(input)

############ PROCESSING DATA ############
# set column names

# COLLECTION GRAPH NVTXS NEDGES AVGDEG MINDEG MAXDEG METHOD NPARTS NTHREADS SEED UB EDGECUT NBISLANDS TIME CHECK
# method  collection   graph   topo  seeds  duration   maxCongestion    maxDilation    avgDilation  Cut   TCV
names(ALLDATA) <- c('FILE','COLLECTION','METHOD','TOPOLOGY','SEED','DURATION','MAXCON','MAXDIL','AVGDIL','CUT','TCV')

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
#METHODREF <- "SCOTCH"
#METHODREF <- "KAHIP"
cat("METHODREF = ", METHODREF, "\n")

############# IGNORE METHOD ############
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "GLOBAL_CLASSIC_100", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "GLOBAL_DIFF_100", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "GLOBAL_HYBRID_100", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_CLASSIC_100", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_CLASSIC_10", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_HYBRID_100", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_HYBRID_10", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_DIFF_1", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_DIFF_100", ]
## ALLDATA <- ALLDATA[ALLDATA$METHOD != "LOCAL_DIFF_10", ]

# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(COLLECTION, FILE, METHOD, TOPOLOGY),
                  transform,
                  DURATIONMEAN=mean(DURATION,na.rm=TRUE),
                  MAXCONMEAN=mean(MAXCON,na.rm=TRUE),
                  MAXDILMEAN=mean(MAXDIL,na.rm=TRUE),
                  AVGDILMEAN=mean(AVGDIL,na.rm=TRUE),
                  CUTMEAN=mean(CUT,na.rm=TRUE),
                  TCVMEAN=mean(TCV,na.rm=TRUE))

NORMDATA[c('FILE','METHOD','TOPOLOGY','CUTMEAN','DURATIONMEAN')]
NORMDATA <- ddply(NORMDATA, .(COLLECTION, FILE, TOPOLOGY),
                  transform,
                  CUTR=CUT/CUTMEAN[METHOD==METHODREF],
                  MAXCONR=MAXCON/MAXCONMEAN[METHOD==METHODREF],
                  MAXDILR=MAXDIL/MAXDILMEAN[METHOD==METHODREF],
                  AVGDILR=AVGDIL/AVGDILMEAN[METHOD==METHODREF],
                  TCVR=TCV/TCVMEAN[METHOD==METHODREF],
                  DURATIONR=DURATION/DURATIONMEAN[METHOD==METHODREF])
# dump
NORMDATA[c('FILE','COLLECTION','METHOD','TOPOLOGY','CUTR','DURATIONR')]

# remove useless columns
NORMDATA[c("CUT","DURATION","CUTMEAN","DURATIONMEAN")] <- list(NULL) 

# compute MEAN & SD & SE and summarise dataframe
MEANDATA <- ddply(NORMDATA,
                  .(COLLECTION, METHOD, TOPOLOGY), # compute separate mean for each COLLECTION x METHOD x NPARTS
                  summarise,
                  CUT=mean(CUTR,na.rm=TRUE),
                  CUTSD=sd(CUTR,na.rm=TRUE),
                  DURATION=mean(DURATIONR,na.rm=TRUE),
                  DURATIONSD=sd(DURATIONR,na.rm=TRUE),
                  MAXCON=mean(MAXCONR,na.rm=TRUE),
                  MAXCONSD=sd(MAXCONR,na.rm=TRUE),
                  MAXDIL=mean(MAXDILR,na.rm=TRUE),
                  MAXDILSD=sd(MAXDILR,na.rm=TRUE),       
                  AVGDIL=mean(AVGDILR,na.rm=TRUE),
                  AVGDILSD=sd(AVGDILR,na.rm=TRUE),       
                  TCV=mean(TCVR,na.rm=TRUE),
                  TCVSD=sd(TCVR,na.rm=TRUE)       
                  )


############ DUMP RESULTS ############
MEANDATA[c('COLLECTION','METHOD','TOPOLOGY','CUT','CUTSD','DURATION','DURATIONSD')] 

############ REORDER ALL METHODS ############

# reorder methods
# MEANDATA$METHOD <- factor(MEANDATA$METHOD, levels = c("SCOTCH","SCOTCHFM","KGGGP_G", "KGGGP_L","SCOTCHK","SCOTCHKK","SCOTCHK_C","SCOTCHK_L","METISKGGGP","METIS","KMETIS","PATOH","ZOLTAN"))

############ MY OWN COLORS ############

### configure my own palette based on brewer one

library(RColorBrewer)
PALSIZE <- max(length(unique(ALLDATA$METHOD)),3)
MYPAL <- brewer.pal(PALSIZE, "Blues")
MYPAL[1] <- "lightblue"
MYPAL[2] <- "blue"
MYPAL[3] <- "purple"
MYPAL[4] <- "yellow"
MYPAL[5] <- "magenta"
MYPAL[6] <- "red"
MYPAL[7] <- "green"
MYPAL[8] <- "pink"
MYPAL[9] <- "lightgreen"
# "darkgreen" 

############ BARPLOT ############

# display graphics with lines and errorbar


cutplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=CUT, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=CUT-CUTSD, ymax=CUT+CUTSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") +
        ## facet_grid(. ~ COLLECTION, scales="free") +
    facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") +    
    xlab("topologies") + ylab("relative cut") +
    coord_cartesian(ylim = c(0.5, 1.8)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  


   durplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=DURATION, fill=METHOD, colour=METHOD, group=METHOD)) +
       geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
       geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
       geom_errorbar(aes(ymin=DURATION-DURATIONSD, ymax=DURATION+DURATIONSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
       ##    facet_grid(. ~ COLLECTION, scales="free") +
       facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") +        
       scale_color_manual(values=MYPAL, name="Methods") + 
       xlab("topologies") + ylab("relative dur") +
                                        #    coord_cartesian(ylim = c(0.5, 1.3)) +
                                        #    scale_fill_brewer() +
       scale_fill_manual(values=MYPAL, name="Methods") +  
                                        #    scale_fill_manual(name="Methods",values=d_colors_methods) +
       theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
       scale_y_continuous(breaks=c(0.80,1.00,1.20))


   maxconplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=MAXCON, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=MAXCON-MAXCONSD, ymax=MAXCON+MAXCONSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
       facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") +        
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
    geom_errorbar(aes(ymin=MAXDIL-MAXDILSD, ymax=MAXDIL+MAXDILSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
       facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") +        
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
    geom_errorbar(aes(ymin=AVGDIL-AVGDILSD, ymax=AVGDIL+AVGDILSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
       facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") +        
    scale_color_manual(values=MYPAL, name="Methods") + 
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
    geom_errorbar(aes(ymin=TCV-TCVSD, ymax=TCV+TCVSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
       facet_wrap( ~ COLLECTION, ncol=2, scales="fixed") +        
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative tcv") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  





pdf(outputplot)

cutplot
durplot
maxconplot
maxdilplot
avgdilplot
tcvplot
cat("writing plot in ", outputplot, "\n");

### EOF
