#' #' @title ChaR-Analysis
#' #' @description Diagnostic and anlystical tools for peak analysis in
#' #' sediment-charcoal records. This code is based on original Matlab code,
#' #' CharAnalysis reads in a data object or a filename and uses it to
#' #' parameterize the detection of charcoal peaks from background values.
#' #' 
#' #'  \item{ \code{charcoal} }{}
#' #'  \item{ \code{char.thresh} }{}}
#' #'  \item{ \code{pre.treatment} }{}
#' #'  \item{ \code{smoothing} }{}
#' #'  \item{ \code{peak.analysis} }{}
#' #'  \item{ \code{results} }{}
#' #' 
#' #' @section Note:
#' #' 
#' #' @examples \dontrun{
#' #' 
#' #' }
#' #' @references
#' #' CharAnalysis documentation: http://code.google.com/p/charanalysis/
#' #' @keywords IO connection
#' #' @export
#' 

#CharAnalysis <- function(site.name="CO", runname=NULL) {  

## Parameters that should go into main function, used to test here the code...
## SHOULD BE DELETED ONCE EVERYTHING IS READY
rm(list=ls())
setwd('/Users/wfinsing/Documents/GitHub/ChaR-Analysis/Cores')
site.name <- "CO"
runname <- NULL
runname <- "1"

#### Load packages
require(paleofire) # function: pretreatment()
require(caTools) # function: runmean()

# 0. Setup directories ####

#### Determine input directory
input.dir <- file.path("..", site.name)

#### Create output directory
if (is.null(runname)) {
#output.dir <- file.path(".", "Cores", site.name, "output")
output.dir <- file.path(".", site.name, "output")
} else {
  output.dir <- file.path(".", site.name, paste0("output", runname))
}

if (!dir.exists(output.dir)) {
  dir.create(output.dir)
} else {
  warning(paste("Note: The output directory", output.dir, "/n...already exists and was overwritten with new output!"))
}

# 1. Load Charcoal data ####
cat(' (1) Reading input files...')
Charcoal <- read.csv(file.path(".", site.name, paste0(site.name, "_charData.csv")))

# Load Parameters file
Params   <- read.csv(file.path(".", site.name, paste0(site.name, "_charParams.csv")),
                     header=T,
                     colClasses = c("NULL", "factor", "numeric", "NULL", "NULL"))

# Extract data from Parameters file
zones <- Params[1:9, ]
zones <- na.omit(zones)

zones <- zones[which(complete.cases(zones)), ]
zones <- zones[ ,2]

yr.interp     <- if(Params[10,2] == 0) {
  yr.interp = NULL
}

char.tr       <- Params[11,2]
char.sm.meth  <- Params[12,2]
char.sm.yr   <- Params[13,2]
cPeak         <- Params[14,2]
thresh.type   <- Params[15,2]
thresh.meth   <- Params[16,2]
thresh.values <- Params[17:20,2]
minCountP     <- Params[21,2]
peakFrequ     <- Params[22,2]

# Clean Environment
rm(Params)

# 2. Pretreatment ####
cat('\n (2) Pretreating charcoal data...')

# Get data and run pretreatment() function
char.series <- Charcoal[ ,6]
char.params <- Charcoal[ ,1:5]
Charcoal.I <- pretreatment(params = char.params, serie = char.series, Int = T,
                           first <- zones[1], last <- zones[length(zones)],
                           yrInterp = yr.interp)

# Clean Environment
rm(char.params, char.series)




# 3. Smooth Charcoal.I to estimate Low-frequency trends (i.e. Char.background) ####
cat('\n (3) Smoothing resampled CHAR to estimate low-frequency trends...')
cat('\n     and calculating peak CHAR...')

# Get data needed
yr.interp <- Charcoal.I$yrInterp
n.smooth <- round(char.sm.yr/yr.interp) # % number of data points over which to smooth the record.
span <- n.smooth/length(Charcoal.I$accI)

# Prepare empty data.frame
charAccIS <- data.frame(matrix(NA, nrow=length(Charcoal.I$cmI), ncol=6))

# # Lowess
#   charAccIS[ ,1] <- lowess(x=Charcoal.I$accI, f=span, iter = 3)$y

# Loess with default options
  in.loess <- data.frame(Charcoal.I$ybpI, Charcoal.I$accI)
  charAccIS[ ,1] <- loess(formula = Charcoal.I.accI ~ Charcoal.I.ybpI, data = in.loess,
                          span = span)$fitted
  
# Robust Loess
  charAccIS[ ,2] <- loess(formula = Charcoal.I.accI ~ Charcoal.I.ybpI, data = in.loess,
                     span = span, degree = 2, family="symmetric",
                     control=loess.control(iterations=4))$fitted

  rm(in.loess)
  
# Moving average
  charAccIS[ ,3] <- runmean(x=Charcoal.I$accI, k=n.smooth, alg="exact", endrule="mean",
                            align="center")

# Running median
  if (n.smooth %% 2 == 0) {
    s.smooth.rmed <- n.smooth-1
  } else {
    s.smooth.rmed <- n.smooth
  }
  r <- as.vector(runmed(x=Charcoal.I$accI, k=s.smooth.rmed, endrule="median"))
  charAccIS[ ,4] <- as.vector(runmed(x=Charcoal.I$accI, k=s.smooth.rmed, endrule="median"))
  
  rm(s.smooth.rmed)
  
# Running mode
  # To be done!
  

  
# Plot raw char and smoothed series
  
  # Set axis limits
  x.min <- floor(min(Charcoal.I$ybpI)/100) * 100
  x.max <- ceiling(max(Charcoal.I$ybpI)/100) * 100
  if(x.max - x.min > 1000) {
  x.by <- 1000
  } else {
    x.by <- round((x.max - x.min)/100) * 100
  }

  x.lim <- c(x.max, x.min)
  y.lim <- c(min(Charcoal.I$accI), 1.2*max(Charcoal.I$accI))
  
  # Plot
  pdf(paste0(output.dir, "/01_pretreatment.pdf"))
  par(mfrow=c(2,1), mar=c(0.5, 4, 0.5, 1), oma=c(5,1,1,1), cex=0.7)
  plot(Charcoal.I$ageTop, Charcoal.I$acc, type="h", col=grey(0.7), xlim=x.lim, ylim=y.lim,
       xlab="time (cal yr. BP)", ylab="CHAR (# cm^-2 yr^-1)", axes=F)
  par(new=T)
  plot(Charcoal.I$ybpI, Charcoal.I$accI, type="s", lwd=0.5, xlim=x.lim, ylim=y.lim,
       ylab="", axes=F)
  axis(1, at=seq(0, x.max, by=1000), labels=F)
  axis(2)
  legend("topleft", inset=c(0,0.2),
         legend=c("CHAR raw", "CHAR interpolated"),
         lwd=1.5, col=c("grey", "black"),
         border="NA", bty="n", horiz=F, cex=0.7)
  mtext(paste("a) ", site.name, ": C.raw and C.interpolated", sep=''), cex=0.8, line=-4)
  
  plot(Charcoal.I$ybpI, Charcoal.I$accI, type="s", xlim=x.lim, ylim=y.lim,
       xlab="", ylab="CHAR (# cm^-2 yr^-1)", lwd=1.5, axes=F)
  polygon(x=c(rev(Charcoal.I$ybpI), Charcoal.I$ybpI),
          y=c(rep(0, length(Charcoal.I$accI)), Charcoal.I$accI),
          col=gray(0.7), border=T, lwd=0.5)
  lines(Charcoal.I$ybpI, charAccIS[ ,1], type="l", col="red", lwd=1.5)
  lines(Charcoal.I$ybpI, charAccIS[ ,2], type="l", col="green", lwd=1.5)
  lines(Charcoal.I$ybpI, charAccIS[ ,3], type="l", col="blue", lwd=1.5)
  lines(Charcoal.I$ybpI, charAccIS[ ,4], type="l", col="orange", lwd=1.5)
  lines(Charcoal.I$ybpI, charAccIS[ ,5], type="l", col="mistyrose", lwd=1.5)
  legend("topleft", inset=c(0,0.2),
         legend=c("CHARinterpolated", "Lowess","Loess","Robust Loess","Moving average","Running median"),
         lwd=1.5,
         col=c("grey", "red","green","blue","orange","mistyrose"),
         border="NA", bty="n", horiz=F, cex=0.7)
  axis(1, at=seq(0, x.max, by=1000))
  axis(2)
  mtext(paste("b) CHARinterpolated and different options for a ", char.sm.yr,"yr CHARbackground",
              sep=''), cex=0.8, line=-4)
  dev.off()

# Calculate peak CHAR component by removing background CHAR
  if (min(Charcoal.I$acc == 0) && cPeak == 2) {
    warning('Cannot calculate C_peak when C_background values = 0; change parameters.')
  }
  
  if (cPeak == 1) {
    Charcoal.peak <- Charcoal.I$accI - charAccIS[ ,char.sm.meth]
  }
  if (cPeak == 2) {
    Charcoal.peak <- Charcoal.I$accI / charAccIS[ ,char.sm.meth]
  }

  cat('\n...done.')
  
  ## 4. Define possible threshold for peak identification
  cat("(4) Defining possible thresholds for peak identification...")
  
#}