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

############ PROCESSING DATA ############

# set column names
# COLLECTION GRAPH NVTXS NEDGES AVGDEG MINDEG MAXDEG METHOD NPARTS NTHREADS SEED UB EDGECUT NBISLANDS TIME CHECK
# method  collection   graph   topo  seeds  duration   maxCongestion    maxDilation    avgDilation  Cut   TCV
# alma data
names(ALLDATA) <- c('FILE','COLLECTION','METHOD','TOPOLOGY','SEED','DURATION','MAXCON','MAXDIL','AVGDIL','CUT','TCV')

# remove COLLECTION column
ALLDATA$FILE <- paste(ALLDATA$COLLECTION, ALLDATA$FILE, sep="_")
ALLDATA$COLLECTION <- list(NULL) 

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
METHODREF <- "SCOTCH"
#METHODREF <- "KAHIP"
cat("METHODREF = ", METHODREF, "\n")

ALLDATA[c('FILE','METHOD','TOPOLOGY','CUT','DURATION')]


# normalize data: TIME and EDGECUT relative to the reference method average
# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(FILE, METHOD, TOPOLOGY),
                  transform,
                  DURATIONMEAN=mean(DURATION,na.rm=TRUE),
                  MAXCONMEAN=mean(MAXCON,na.rm=TRUE),
                  MAXDILMEAN=mean(MAXDIL,na.rm=TRUE),
                  AVGDILMEAN=mean(AVGDIL,na.rm=TRUE),
                  CUTMEAN=mean(CUT,na.rm=TRUE),
                  TCVMEAN=mean(TCV,na.rm=TRUE))


NORMDATA[c('FILE','METHOD','TOPOLOGY','CUTMEAN','DURATIONMEAN')]
# normalize the results based on the reference method
NORMDATA <- ddply(NORMDATA, .(FILE, TOPOLOGY),
                  transform,
                  CUTR=CUT/CUTMEAN[METHOD==METHODREF],
                  MAXCONR=MAXCON/MAXCONMEAN[METHOD==METHODREF],
                  MAXDILR=MAXDIL/MAXDILMEAN[METHOD==METHODREF],
                  AVGDILR=AVGDIL/AVGDILMEAN[METHOD==METHODREF],
                  TCVR=TCV/TCVMEAN[METHOD==METHODREF],
                  DURATIONR=DURATION/DURATIONMEAN[METHOD==METHODREF])

NORMDATA[c('FILE','METHOD','TOPOLOGY','CUTR','DURATIONR')]

# remove useless columns
NORMDATA[c("CUT","DURATION","CUTMEAN","DURATIONMEAN")] <- list(NULL) 



# compute MEAN & SD & SE over all graphs and summarise data
MEANDATA <- ddply(NORMDATA,
                  .(METHOD, TOPOLOGY), # compute separate mean for each METHOD x NPARTS
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

# add TOPOLOGY column to use TOPO as a factor rather than numeric
# NORMDATA$NPARTSF <- factor(NORMDATA$NPARTS) 
# MEANDATA$NPARTSF <- factor(MEANDATA$NPARTS) 

############ REORDER ALL METHODS ############

# MEANDATA$METHOD <- factor(MEANDATA$METHOD, levels = c("SCOTCH","SCOTCHFM","KGGGP_G", "KGGGP_L","SCOTCHK","SCOTCHK_C","SCOTCHK_L","METIS","KMETIS","PATOH","ZOLTAN","KGGGPML","KGGGPML_RM","KGGGPML_SHEM","METISKGGGP","METISKGGGP_L"))

############ REMOVE A METHOD ############

## MEANDATA <- MEANDATA[MEANDATA$METHOD != 'SCOTCH', ]

############ DUMP RESULTS ############

# dump meandata
MEANDATA[c('METHOD','TOPOLOGY','CUT','CUTSD','DURATION','DURATIONSD')] 

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
   cutplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=CUT, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=CUT-CUTSD, ymax=CUT+CUTSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Methods") + 
    xlab("topologies") + ylab("relative cut") +
#    coord_cartesian(ylim = c(0.5, 1.3)) +
#    scale_fill_brewer() +
     scale_fill_manual(values=MYPAL, name="Methods") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
     theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  
# to remove diagonal lines in legend.key: + guides(fill = guide_legend(override.aes = list(colour = NULL)))


   durplot <- ggplot(MEANDATA, aes(x=TOPOLOGY, y=DURATION, fill=METHOD, colour=METHOD, group=METHOD)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    geom_errorbar(aes(ymin=DURATION-DURATIONSD, ymax=DURATION+DURATIONSD), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
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
MEANDATA[c("METHOD")] <- list(NULL)
sdfm <- MEANDATA[c('TOPOLOGY','CUTSD','DURATIONSD','MAXCONSD','AVGDILSD','MAXDILSD','TCVSD')]
# do not include avgdil
sdfm[c("AVGDILSD")] <- list(NULL)
# do not include duration either
sdfm[c("DURATIONSD")] <- list(NULL)

MEANDATA$DURATION <- MEANDATA$DURATION + 1
MEANDATA[c("CUTSD","DURATIONSD","MAXCONSD","AVGDILSD","MAXDILSD","TCVSD")] <- list(NULL)
# do not include avgdil
MEANDATA[c("AVGDIL")] <- list(NULL)
MEANDATA[c("DURATION")] <- list(NULL)
dfm <- melt(MEANDATA, id.vars = "TOPOLOGY")
cat("after melting dfm:\n")
#dfm[c("METHOD")] <- list(NULL) 
dfm
sdfm <- melt(sdfm, id.vars = "TOPOLOGY")
cat("sd data after melting sdfm:\n")
#dfm[c("METHOD")] <- list(NULL) 
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



splot <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    #geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    #geom_errorbar(aes(ymin=value-sdfm$value, ymax=value+sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    geom_point(position=position_dodge(width=0.9), size = 2) +
    geom_errorbar(aes(ymin=value-sdfm$value, ymax=value+sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.5, 1.5)) +
#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
#+        		     scale_y_continuous(breaks=c(0.80,1.00,1.20))  


bplot <- ggplot(dfm, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    # bars 
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_errorbar(aes(ymin=value-sdfm$value, ymax=value+sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    # scatter
    #geom_point(position=position_dodge(width=0.9), size = 2) +
    #geom_errorbar(aes(ymin=value-sdfm$value, ymax=value+sdfm$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +

    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
    scale_color_manual(values=MYPAL, name="Metrics") + 
    xlab("topologies") + ylab("relative results") +
    coord_cartesian(ylim = c(0.5, 1.5)) +
#    scale_fill_brewer() +
    scale_fill_manual(values=MYPAL, name="Metrics") +  
#    scale_fill_manual(name="Methods",values=d_colors_methods) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
#+        		     scale_y_continuous(breaks=c(0.80,1.00,1.20))  



pdf(output)
cutplot
durplot
maxconplot
maxdilplot
avgdilplot
tcvplot
splot
bplot
cat("writing plot in ", output, "\n");

### EOF
