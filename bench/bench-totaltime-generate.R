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
# method  collection   graph   topo  seeds  duration   maxCongestion    maxDilation    avgDilation  Cut   Co
# alma data
names(ALLDATA) <- c('FILE','COLLECTION','METHOD','TOPOLOGY','SEED','Time','MAXCON','MAXDIL','AVGDIL','Cut','Co')

# remove COLLECTION column
ALLDATA$FILE <- paste(ALLDATA$COLLECTION, ALLDATA$FILE, sep="_")
ALLDATA$COLLECTION <- list(NULL) 

# select reference method
METHODREF <- as.vector(ALLDATA$METHOD)[1]
#METHODREF <- "SCOTCH"
#METHODREF <- "KAHIP"
#cat("METHODREF = ", METHODREF, "\n")

#ALLDATA[c('FILE','METHOD','TOPOLOGY','Cut','Time')]

ALLDATA[c("Cut","MAXDIL","AVGDIL","Co","MAXCON")] <- list(NULL) 



# normalize data: TIME and EDGECut relative to the reference method average
# compute mean over seeds
NORMDATA <- ddply(ALLDATA, .(FILE, METHOD, TOPOLOGY),
                  transform,
                  Time=mean(Time,na.rm=TRUE) #,
		  #minTime=Time,
		  #maxTime=Time
		  )

METHODREF <- as.vector(NORMDATA$METHOD)[1]
METHODS <- as.vector(NORMDATA$METHOD)[2]


# here I have the arithmetic mean, min and max values over the seeds for each graphxtopo for both methods
# dump data
#NORMDATA[c('FILE','METHOD','TOPOLOGY','Time')]

NORMDATA <- ddply(NORMDATA,
                  .(FILE, TOPOLOGY, SEED), # compute separate mean for each METHOD x NPARTS
                  transform,
		  sTime = Time[METHOD==METHODREF]
		  )

# dump data
#NORMDATA[c('FILE','METHOD','TOPOLOGY','Time','sTime')]


NORMDATA<- NORMDATA[NORMDATA$METHOD != METHODREF, ]


NORMDATA <- ddply(NORMDATA,
                  .( METHOD , FILE, TOPOLOGY, SEED), # compute separate mean for each METHOD x NPARTS
                  transform,
		  tTime = sTime + Time
		  )

NORMDATA[c("SEED","METHOD","COLLECTION")] <- list(NULL) # for scotch remove also topo and file ,"TOPOLOGY","FILE"
# remove linenumbers
# dump data 
#print(NORMDATA, row.names = FALSE)

NORMDATA[c("Time","sTime")] <- list(NULL) 
print(NORMDATA, row.names = FALSE)


