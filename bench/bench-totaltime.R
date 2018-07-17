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
names(ALLDATA) <- c('FILE','TOPOLOGY','IdentityTiMErraw','ScotchTiMErraw')


ALLDATA[c('FILE','TOPOLOGY','IdentityTiMErraw','ScotchTiMErraw')]


# # normalize data:
ALLDATA <- ddply(ALLDATA,
	 .(FILE, TOPOLOGY),
                  transform,
		  IdentityTiMErraw=IdentityTiMErraw/ScotchTiMErraw,
		  ScotchTiMErraw=ScotchTiMErraw/ScotchTiMErraw
		  )

ALLDATA[c('FILE','TOPOLOGY','IdentityTiMErraw','ScotchTiMErraw')]

# compute geom mean
MEANDATA <- ddply(ALLDATA,
                  .(TOPOLOGY), # compute separate mean for each METHOD x NPARTS
         summarise,	       

		  IdentityTiMEr=exp(mean(log(IdentityTiMErraw))),
		  IdentityTiMErSD=exp(sd(log(IdentityTiMErraw))),
		  ScotchTiMEr=exp(mean(log(ScotchTiMErraw))),
		  ScotchTiMErSD=exp(sd(log(ScotchTiMErraw)))
		  # IdentityTiMEr=mean(IdentityTiMErraw),
		  # IdentityTiMErSD=sd(IdentityTiMErraw),
		  # ScotchTiMEr=mean(ScotchTiMErraw),
		  # ScotchTiMErSD=sd(ScotchTiMErraw)
                  )


# remove useless columns
MEANDATA[c("IdentityTiMErraw","ScotchTiMErraw")] <- list(NULL) 




MEANDATA[c('TOPOLOGY', 'ScotchTiMEr','IdentityTiMEr',"ScotchTiMErSD","IdentityTiMErSD")]

library(reshape2)

 dd<- MEANDATA[c('TOPOLOGY','IdentityTiMEr','ScotchTiMEr')]
 sd<- MEANDATA[c('TOPOLOGY','IdentityTiMErSD','ScotchTiMErSD')]
 
 dd <- melt(dd, id.vars = "TOPOLOGY")
 sd <- melt(sd, id.vars = "TOPOLOGY")


dd
sd

 library(RColorBrewer)
 PALSIZE <- max(length(unique(ALLDATA$TOPOLOGY)),3)
 MYPAL <- brewer.pal(PALSIZE, "Set1")
 MYPAL  <- brewer.pal(PALSIZE, "Spectral")


 plot <- ggplot(dd, aes(x=TOPOLOGY, y=value, fill=variable, colour=variable, group=variable)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.9, size=0.2, color="black") +
    geom_hline(yintercept=1.0, color="red", size=0.2, linetype="dashed") +
    geom_errorbar(aes(ymin=value/sd$value, ymax=value*sd$value), position=position_dodge(width=0.9), colour="black", width=0.5, size=0.2) +
    scale_color_manual(values=MYPAL, name="Time") + 
    xlab("topologies") + ylab("relative runtime") +
    #coord_cartesian(ylim = c(0.0, 5)) +	  
    scale_fill_manual(values=MYPAL, name="Methods") +  
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
    #scale_y_continuous(breaks=c(0.5,1.00,1.50,2.00))
    # for scotch
    #scale_y_continuous(breaks=c(1.00,5.00,10.00,20.00,50.00))	

#summarize(MEANDATA, tmeanovertopo = mean(IdentityTiMEr, na.rm = T))


mean(MEANDATA[["IdentityTiMEr"]])



pdf(output)
plot

cat("writing plot in ", output, "\n");