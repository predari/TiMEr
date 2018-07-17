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

names(ALLDATA) <- c('FILE','COLLECTION','METHOD','UBFACTOR','TOPO','SEED','TIME')

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
## METHODREF <- "SCOTCHFM"
cat("METHODREF = ", METHODREF, "\n")

# ubfactor
UBFACTOR <- 3
cat("UBFACTOR = ", UBFACTOR, "\n")

#dump
ALLDATA[c('FILE','TOPO','UBFACTOR','TIME')] 

# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(FILE, METHOD, TOPO),
                  transform,
                  TIMEMEAN=mean(TIME,na.rm=TRUE))
#dump
NORMDATA[c('FILE','TOPO','TIMEMEAN')]


# remove useless columns
NORMDATA[c("TIME")] <- list(NULL) 

 # compute MEAN and summarise dataframe so you have only one line for each (graph x topo)
 MEANDATA <- ddply(NORMDATA, .( TOPO, FILE),	
                   summarise,
		   TIMEMEAN=mean(TIMEMEAN,na.rm=TRUE))
#                   ECMEAN=mean(ECR,na.rm=TRUE),
#                   ECSD=sd(ECR,na.rm=TRUE),
#                   TIMEMEAN=mean(TIMER,na.rm=TRUE),
#                   TIMESD=sd(TIMER,na.rm=TRUE),
#                   UBMEAN=mean(UB,na.rm=TRUE),
#                   UBSD=sd(UB,na.rm=TRUE))
# #                  FAILSUM=sum(FAIL))

# add NPARTSF column to use NPARTS as a factor rather than numeric
#NORMDATA$NPARTSF <- factor(NORMDATA$NPARTS) 
#MEANDATA$NPARTSF <- factor(MEANDATA$NPARTS) 

# dump
MEANDATA[c('FILE','TOPO','TIMEMEAN')]

# dump only time for each (topo x graph)
MEANDATA[c('TIMEMEAN')]


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

timeplot <- ggplot(NORMDATA, aes(x=TOPO, y=TIMEMEAN, fill=FILE, colour=FILE, group=FILE)) +
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
