rm(list=ls())

# Load data

#### Load packages
require(paleofire) # function: pretreatment()
require(caTools) # function: runmean()
require(mclust)
require(Matching)

setwd("/Users/wfinsing/Documents/GitHub/ChaR-Analysis copy")

## Load R source files
source("./R/pretreatment_full.r")
source("./R/CHARsm.r")
source("./R/CharThreshLocal.r")
source("./R/CharPeaks.r")


# 1. Load Charcoal data ####
CO <- read.csv("./Cores/CO/CO_charData.csv")

# 2. Pretreatment
CO.I <- pretreatment_full(params=CO[ ,1:5], serie=CO[ ,6], Int=T,
                                 first=-42, last=7444, yrInterp=15)
#plot(CO.I)

# 3. Smooth record
CO.sm <- CHARsm(x=CO.I, smoothing.yr=1000)
plot(CO.sm)

# 4. Peak detection (Local threshold)
CO.thresh <- CharThreshLocal(x=CO.sm, sm.meth="Robust Lowess", thresh.yr=1000,
                             thresh.values=c(0.95, 0.99,0.999,0.99), plot.thresh=T)

# 5. Identify peaks based on Cback and on local threshold
CO.peak <- CharPeaks(x=CO.thresh, minCountP=0.05)
plot(CO.peak)
plot.SNI(CO.peak)
