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

############ LOADING DATA ############

# load raw data
ALLDATA <- read.table(input)

############ PROCESSING DATA ############

names(ALLDATA) <- c('FILE','METHOD','UBFACTOR','NPARTS','SEED','EC','TIME')

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
## METHODREF <- "SCOTCHFM"
cat("METHODREF = ", METHODREF, "\n")

# ubfactor
UBFACTOR <- 3
cat("UBFACTOR = ", UBFACTOR, "\n")

#dump
ALLDATA[c('FILE','NPARTS','EC','TIME')] 

# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(FILE, METHOD, NPARTS),
                  transform,
                  ECMEAN=mean(EC,na.rm=TRUE),
                  TIMEMEAN=mean(TIME,na.rm=TRUE))
#dump
NORMDATA[c('FILE','NPARTS','ECMEAN','TIMEMEAN')]


# remove useless columns
NORMDATA[c("EC","TIME")] <- list(NULL) 

 # compute MEAN and summarise dataframe so you have only one line for each (graph x topo)
 MEANDATA <- ddply(NORMDATA, .( NPARTS, FILE, NPARTS ),	
                   summarise,
                   ECMEAN=mean(ECMEAN,na.rm=TRUE),
		   TIMEMEAN=mean(TIMEMEAN,na.rm=TRUE))
#                   ECMEAN=mean(ECR,na.rm=TRUE),
#                   ECSD=sd(ECR,na.rm=TRUE),
#                   TIMEMEAN=mean(TIMER,na.rm=TRUE),
#                   TIMESD=sd(TIMER,na.rm=TRUE),
#                   UBMEAN=mean(UB,na.rm=TRUE),
#                   UBSD=sd(UB,na.rm=TRUE))
# #                  FAILSUM=sum(FAIL))

# add NPARTSF column to use NPARTS as a factor rather than numeric
NORMDATA$NPARTSF <- factor(NORMDATA$NPARTS) 
MEANDATA$NPARTSF <- factor(MEANDATA$NPARTS) 

# dump
MEANDATA[c('FILE','NPARTS','ECMEAN','TIMEMEAN')]

# dump only time for each (topo x graph)
MEANDATA[c('TIMEMEAN')]

#print(MEANDATA, row.names = FALSE)

MEANDATA256 <- MEANDATA[MEANDATA$NPARTS == '256', ]
MEANDATA256[c("NPARTS","ECMEAN","NPARTSF")] <- list(NULL) 
print(MEANDATA256, row.names = FALSE)

MEANDATA512 <- MEANDATA[MEANDATA$NPARTS == '512', ]
MEANDATA512
MEANDATA512[c("NPARTS","ECMEAN","NPARTSF")] <- list(NULL) 
print(MEANDATA512, row.names = FALSE)

mean(MEANDATA256[["TIMEMEAN"]])
exp(mean(log(MEANDATA256[["TIMEMEAN"]])))

mean(MEANDATA512[["TIMEMEAN"]])
exp(mean(log(MEANDATA512[["TIMEMEAN"]])))

library(RColorBrewer)
PALSIZE <- max(length(unique(ALLDATA$FILE)),13)

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



############ BARPLOT ############

## barplot <- TRUE

# display graphics with lines and errorbar
   edgecutplot <- ggplot(NORMDATA, aes(x=NPARTS, y=ECMEAN, fill=FILE, colour=FILE, group=FILE)) +
     geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    scale_color_manual(values=MYPAL, name="Files") + 
    xlab("nb parts") + ylab("relative edgecut") +
     scale_fill_manual(values=MYPAL, name="Files") +  
     theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) +        		
     scale_y_continuous(breaks=c(0.80,1.00,1.20))  
     
# to remove diagonal lines in legend.key: + guides(fill = guide_legend(override.aes = list(colour = NULL)))

timeplot <- ggplot(NORMDATA, aes(x=NPARTS, y=TIMEMEAN, fill=FILE, colour=FILE, group=FILE)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +        
   xlab("nb parts") + ylab("relative runtime") +
   coord_cartesian(ylim = c(0.5, 4.0)) +            
   scale_fill_manual(values=MYPAL, name="Files") +           
   theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) 

pdf(output)
edgecutplot
timeplot
cat("writing plot in ", output, "\n");


### EOF
