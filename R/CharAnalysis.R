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
  require(mclust)
  
  ## Load R source files
  source("../R/pretreatment_edits.r")
  
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
  cat(' (1) Reading charcoal-data file...')
  Charcoal <- read.csv(file.path(".", site.name, paste0(site.name, "_charData.csv")))
  
  # Load Parameters file
  cat('     ...and reading parameters file...')
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
  Charcoal.I <- pretreatment_edits(params = char.params, serie = char.series, Int = T,
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
  if (n.smooth %% 2 == 0) { # if n.smooth is not an odd number
    s.smooth.rmed <- n.smooth-1
  } else {
    s.smooth.rmed <- n.smooth
  }
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
       xlab="", ylab="CHAR (# cm^-2 yr^-1)", lwd=1, axes=F)
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
  
  
  
  
  
  ## 4. Define possible threshold for peak identification ####
  cat("(4) Defining possible thresholds for peak identification...")
  
  if  (thresh.type == 2) {  # If threshold is defined locally...
    
    # [CharThresh] = CharThreshLocal(Charcoal,...
    #                                Smoothing, PeakAnalysis, site, Results);
    # based on PH 'CharThreshLocal.m'
    # Determines a threshold value for each interpolated sample, based on the
    # distribution of CHAR values within the selected window (yr) and either
    # a Gaussian mixture model or the assumption that the noise component 
    # of the peak charcoal record (C_peak) is normally distributed 
    # around 0 (if C_peak is defined by residuals)
    # or 1 (if C_peak is defined by ratios).
    
    # Create space for local variables
    threshYr <- char.sm.yr # [yr] Years over which to define threshold
    
    CharThresh.pos <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak),
                                        ncol=length(thresh.values))) # space for threshold values
    CharThresh.neg <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak),
                                        ncol=length(thresh.values))) # Space for negative thres.values
    muHat <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak), ncol=2)) # Space for mean of noise distribution
    sigmaHat <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak), ncol=2)) # Space for standard deviation of noise distribution
    propN <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak), ncol=2)) # Space for proportion of each Cluster-identified distribution
    CharThresh.SNI <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak), ncol=1)) # Space for SNI
    CharThresh.GOF <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak), ncol=1)) # Space for Goodness-of-fit
    
    # Set variables and space for plots
    j <- 1
    num.plots <- seq(from=round(1+(n.smooth/2)), to=length(Charcoal.peak), by=n.smooth)
    my.plots <- vector(length(num.plots), mode='list')
    
    # SELECT Charcoal.peak VALUES TO EVALUATE, BASED ON Smoothing.yr
    for (i in 1:length(Charcoal.peak)) {  #For each value in Charcoal.peak, find the threshold.
      cat(paste0("Calculating ", i, "th local threshold of ", length(Charcoal.peak)))
      
      if (i < round(0.5*(threshYr/yr.interp))+1) { # First 'threshYr' samples.
        #         X = Charcoal.peak(1:round(0.5*(threshYr/r))); % Pre June 2009.
        #         X = Charcoal.peak(1:round(threshYr/r)); % Modified, June 2009, PEH.
        X <- Charcoal.peak[1:round(0.5*(threshYr/yr.interp))+i] # Modified, % June 2009, PEH.
      }
      if (i > (length(Charcoal.peak)-round(0.5*(threshYr/yr.interp)))) {  # Last 'threshYr' samples.
        #             X = Charcoal.peak(length(Charcoal.peak)-...
        #                 round((threshYr/r)):end);   % Pre June 2009.
        X <- Charcoal.peak[(i-round(0.5*(threshYr/yr.interp))):length(Charcoal.peak)]   # Modified, June 2009, PEH.
        # As recommended by RK, this uses samples from 
        # a half-window before i, all the way to end of record.
      }
      if (i >= round(0.5*(threshYr/yr.interp))+1 && i <= (length(Charcoal.peak)-round(0.5*(threshYr/yr.interp)))) {
        # All samples between first and last 'thrshYr' samples.
        X <- Charcoal.peak[(i-round(0.5*(threshYr/yr.interp))) : (i+round(0.5*(threshYr/yr.interp)))]
      }
      
      
      ## ESTIMATE LOCAL NOISE DISTRIBUTION
      if (thresh.meth == 3) { # Estimate noise distribution with Guassian mixture model
        if (sum(X) == 0) {
          cat("NOTE: All C_peak values = 0; cannot fit noise distribution.")
          cat("\n      Mean and standard deviation forced to equal 0.") 
          cat("\n      Consider longer smoothing window or alternative for") 
          cat("\n      threshMethod parameter.")
          muHat[i, ] <- 0
          sigmaHat[i, ] <- 10^-100
          propN[i, ] <- 0
        } else {
          m <- densityMclust(data=X, G=2)
          # deviation of noise distribution using the Gaussian mixture model
          # implemented in the 'mclust' R package. Didn't check if it is exactly comparable to
          # CLUSTER GMM as from http://cobweb.ecn.purdue.edu/~bouman/software/cluster/
          
          # NOTRUN: plots density of two components with default plotting function: 
          #plot(m, what="density", data=X, breaks=50)
          #summary(m, parameters=T, classification=T)
          
          muHat[i, ] <- m$parameters$mean
          sigmaHat[i, ] <- cbind( sqrt(m$parameters$variance$sigmasq))
          propN[i, ] <- m$parameters$pro
          
          if (muHat[i,1] == muHat[i,2]) {
            warning('Poor fit of Gaussian mixture model')
          }
        }
      }
      
      ## Define local threshold, SNI, and GOF
      # Define range of threshold values, plus threshold value selected.
      CharThresh.pos[i, ] <- qnorm(p=thresh.values, mean=muHat[i,1], sd=sigmaHat[i,1])
      CharThresh.neg[i, ] <- qnorm(p=1-thresh.values, mean=muHat[i,1], sd=sigmaHat[i,1])
      sig_i <- X[which(X > qnorm(p=thresh.values[4], mean=muHat[i,1], sd=sigmaHat[i,1]))]
      noise_i <- X[which(X <= qnorm(p=thresh.values[4], mean=muHat[i,1], sd=sigmaHat[i,1]))]
      
      # SNI calculation based on Kelly et al. (2011).
      if (length(sig_i) > 0) {
        CharThresh.SNI[i, ] <- (1/length(sig_i)) * sum((sig_i - mean(noise_i)) / sd(noise_i)) *
          ((length(noise_i)-2)/length(noise_i))
      } else {
        CharThresh.SNI[i, ] <- 0
      }
      colnames(CharThresh.SNI) [1] <- "SNI"
      
      # Evaluate goodness-of-fit between modeled noise distribution and Cnoise samples
      # (Cnoise = CHAR samples less-than or equal to the threshold value) for this time window
      # NB: with ks.test(), 'y' can be a character string naming a continuous (cumulative) distribution function,
      # or such a function (e.g. "pnorm"). In this case, a one-sample test is carried out of the
      # null hypothesis that the distribution function which generated x is distribution y
      # with parameters specified by mean and sd...
      if (length(noise_i) > 3) {
        ksP <- ks.test(x=noise_i, y="pnorm", mean=muHat[i,1], sd=sigmaHat[i,1])$p.value
        CharThresh.GOF[i, ] = ksP
      }
      colnames(CharThresh.GOF) [1] <- "GOF"
      
      # PLOT SELECTED Charcoal.peak DISTRIBUTIONS
      
      
      if (any(i == num.plots)) {
        par(mfrow=c(1,1))
        h <- hist(x=X, breaks=50, plot=F)
        d <- hist(X, breaks=50, plot=F)$density
        pdf1 <- dnorm(x=d, mean=muHat[i, 1], sd=sigmaHat[i, 1])
        
        plot(h, freq=F, col="grey", border="grey", xlim=c(min(h$breaks), max(h$breaks)),
             ylab="Density", xlab='',
             main=paste(Charcoal.I$ybpI[i], "yr BP"))
        par(new=T)
        pdf1 <- curve(dnorm(x=x, mean=muHat[i, 1], sd=sigmaHat[i, 1]),
                      from=min(h$breaks), to=max(h$breaks),
                      ylim=c(0, max(d)), type="l", col="blue", lwd=1.5, axes=F, ylab='', xlab='')
        par(new=T)
        pdf2 <- curve(dnorm(x=x, mean=muHat[i, 2], sd=sigmaHat[i, 2]),
                      from=min(h$breaks), to=max(h$breaks),
                      ylim=c(0, max(d)), type="l", col="orange", lwd=1.5, axes=F, ylab='', xlab='')
        par(new=T)
        lines(x=c(qnorm(p=thresh.values[4], mean=muHat[i,1], sd=sigmaHat[i,1]),
                  qnorm(p=thresh.values[4], mean=muHat[i,1], sd=sigmaHat[i,1])), y=c(0, max(d)),
              type="l", col="red", lwd=1.5)
        mtext(text=paste0("SNIi= ", round(CharThresh.SNI[i, ], digits=2),
                          "\nKS p-val= ", round(CharThresh.GOF[i, ], digits=2)),
              side=3, las=0, line=-2)
        
        my.plots[[j]] <- recordPlot()
        j<- j+1
      }
      
      
    } # end loop for each Charcoal.peak
    
    # Print pdf with selected plots that were saved at the end of the loop above
    pdf(paste0(output.dir, '/02_threshold_determination.pdf'), onefile=TRUE, paper="a4")
    par(mfrow=(c(5,5)), mar = c(0.5,1,0.5,1), oma = c(1,1,0.5,1), cex=0.7)
    for (k in 1:length(num.plots)) {
      #k <- 1
      replayPlot(my.plots[[k]])
    }
    dev.off()
    
    ## Smooth thresholds with Loess smoother
    CharThresh.SNI$sSNI <- loess(formula = CharThresh.SNI[ ,1] ~ Charcoal.I$ybpI, span = span)$fitted
    # CharThresh.SNI (CharThresh.SNI < 0,1) = 0; # Matlab code...what does this mean?
    for (i in 1:length(thresh.values)) {
      CharThresh.pos[ ,i] <- loess(formula = CharThresh.pos[ ,i] ~ Charcoal.I$ybpI, span = span)$fitted
      CharThresh.neg[ ,i] <- loess(formula = CharThresh.neg[ ,i] ~ Charcoal.I$ybpI, span = span)$fitted
    }
    
  } # end part 4. Define possible threshold for peak identification 
  cat('...done.')
  
  
  
  
  
  # 5. Identify charcoal peaks based on possible thresholds ####
  cat('(5) Identifying peaks based on possible thresholds...')
  
  # based on 
  # function [Charcoal, CharThresh] = CharPeakID (Charcoal,Pretreatment,PeakAnalysis,CharThresh);
  #   Identifies charcoal samples that exceeds threshold value(s) determined 
  #   in CharThrshLocal or CharThreshGlobal, and screens these value
  #   according to the minimum-count criterion selected.
  
  ## PEAK IDENTIFICATION ALGORITHM
  # Create space for peaks, Charcoal.charPeaks
  if (thresh.type == 2) {
    # Create Charcoal.charPeaks matrix
    Charcoal.charPeaks <- as.data.frame(matrix(data=0,
                                               nrow=length(Charcoal.peak),
                                               ncol=length(CharThresh.pos)))
    thresholdValues <- CharThresh.pos
  }
  
  nThresholds <- length(thresholdValues)
  
  # Flag values exceeding thresholds
  for (i in 1:length(Charcoal.peak)) {  # For each value in Charcoal.peak
    for (j in 1:nThresholds) {          # For each threshold value
      if (Charcoal.peak[i] > thresholdValues[i, j]) {  # If Charcoal.peak exceeds threshold...
        Charcoal.charPeaks[i, j] <- 2                  # Charcoal.charPeaks = 2
      } else {
        Charcoal.charPeaks[i, j] <- 0                 # else Charcoal.charPeaks = 0
      }
    }
  }
  
  # Remove consecutive Charcoal.charPeaks
  for (i in 1:(length(Charcoal.peak)-1)) { # For each value in Charcoal.peak
    for (j in 1:nThresholds) {             # For each threshold value
      if (Charcoal.charPeaks[i, j] > 0
          && Charcoal.charPeaks[i+1, j] > 0) {  # if two consecutive values > 0 
        Charcoal.charPeaks[i, j] <- 1           # keep first as 2, mark subsequent as 1
      }
    }
  }
  
  for (i in 1:length(Charcoal.peak)) {
    for (j in 1:nThresholds) {
      if (Charcoal.charPeaks[i, j] < 2) {    # if value < 2
        Charcoal.charPeaks[i, j] <- 0        # mark sample as 0 (unflag Charcoal.charPeak)
      } else {
        Charcoal.charPeaks[i, j] <- 1        # else (if value=2) mark sample as 1 (flag as Charcoal.charPeak)
      }
    }
  }
  
  
  # Make a variable to hold the threshold value for each peak identified.
  # i.e. instead of 1 / 0 for peak / no peak, 1 is replaced with the 
  # threshold value used to identify the peak. This is purely for graphing 
  # purposes.
  Charcoal.charPeaksThresh <- Charcoal.charPeaks * 0
  for (i in 1:length(Charcoal.peak)) {
    for (j in 1:nThresholds) {
      Charcoal.charPeaksThresh[i, j] <- Charcoal.charPeaks[i,j] * thresholdValues[i,j]
    }
  }
  
  
  ## Minimum-count Analysis
  mcWindow <- round(150/yr.interp)*yr.interp # [yr] Years before and after a peak to look 
  # for the min. and max. value
  d <- as.data.frame(matrix(data=0, nrow=length(Charcoal.I$accI), ncol=nThresholds))
  CharThresh.minCountP <- as.data.frame(matrix(data=NA, nrow=length(Charcoal.I$acc), ncol=nThresholds))
  alphaPeak <- minCountP
  
  for (j in 1:nThresholds) {
    peakIndex <- which(Charcoal.charPeaks[ ,j] == 1) # Index to find peak samples
    
    if (length(peakIndex) > 1) {                     # Only proceed if there is > 1 peak
      
      for (i in 1:length(peakIndex)) {               # For each peak identified...
        peakYr <- Charcoal.I$ybpI[peakIndex[i]]      # Find age of peak and window around peak
        windowTime <- c(max(Charcoal.I$ybpI[which(Charcoal.I$ybpI <= peakYr+mcWindow)]),
                        min(Charcoal.I$ybpI[which(Charcoal.I$ybpI >= peakYr-mcWindow)]))
        windowTime_in <- c(which(Charcoal.I$ybpI == windowTime[1]), # Index to find range of window ages
                           which(Charcoal.I$ybpI == windowTime[2]))
        if (i == 1) {  # find the year of two adjacent Charcoal.charPeaks, unless first peak,
          #then use windowTime[2] as youngest
          windowPeak_in <- c(which(Charcoal.I$ybpI == Charcoal.I$ybpI[peakIndex[i+1]]),
                             which(Charcoal.I$ybpI == windowTime[2]))
        }
        if (i == length(peakIndex)) {  # if last peak, then use windowTime[1] as oldest age
          windowPeak_in <- c(which(Charcoal.I$ybpI == windowTime[1]),
                             which(Charcoal.I$ybpI == Charcoal.I$ybpI[peakIndex[i-1]]))
        }
        if (i > 1 && i < length(peakIndex)) {
          windowPeak_in <- c(which(Charcoal.I$ybpI == Charcoal.I$ybpI[peakIndex[i+1]]),
                             which(Charcoal.I$ybpI == Charcoal.I$ybpI[peakIndex[i-1]]))
        }
        if (windowTime_in[1] > windowPeak_in[1]) { # thus, if a peak falls within the time window defined by mcWindow
          windowTime_in[1] <- windowPeak_in[1] # replace the windowTime_in with the windowPeak_in
        }
        if (windowTime_in[2] < windowPeak_in[2]) { # thus, if a peak falls within the time window defined by mcWindow
          windowTime_in[2] <- windowPeak_in[2] # replace the windowTime_in with the windowPeak_in
        }
        
        # Final index value for search window: window (1) defines oldest sample,
        # window (2) defines youngest sample
        windowSearch <- c(windowTime_in[1], windowTime_in[2])
        
        # search for max and min Charcoal.charPeaks within this window.
        # [# cm^-3] Max charcoal concentration after peak.
        countMax <- round(max(Charcoal.I$countI[ windowSearch[2]:peakIndex[i] ]))
        # Index for location of max count.
        #countMaxIn <- windowSearch[2]-1 + which(round(Charcoal.I$countI[windowSearch[2]:peakIndex[i]]) == countMax)
        countMaxIn <- windowSearch[2]-1 + max(which(round(Charcoal.I$countI[windowSearch[2]:peakIndex[i]]) == countMax))
        # [# cm^-3] Min charcoal concentration before peak.
        countMin <- round(min(Charcoal.I$countI[peakIndex[i]:windowSearch[1] ]))
        # Index for location of Min count
        #countMinIn <- peakIndex[i]-1 + which(round(Charcoal.I$countI[peakIndex[i]:windowSearch[1]]) == countMin)
        countMinIn <- peakIndex[i]-1 + min(which(round(Charcoal.I$countI[peakIndex[i]:windowSearch[1]]) == countMin))
        
        volMax <- Charcoal.I$volI[countMaxIn]
        volMin <- Charcoal.I$volI[countMinIn]
        d[peakIndex[i], j] <- (abs(countMin-(countMin+countMax)*
                                     (volMin/(volMin+volMax)))-0.5)/(sqrt((countMin+countMax)*
                                                                            (volMin/(volMin+volMax))*(volMax/(volMin+volMax))))
        
        # Test statistic
        CharThresh.minCountP[peakIndex[i], j] <- 1-pt(d[peakIndex[i], j], df=Inf)
                            # Inverse of the Student's T cdf at 
                            # CharThresh.minCountP, with Inf degrees of freedom.
        # From Charster (Gavin 2005):
        # This is the expansion by Shuie and Bain (1982) of the equation by 
        # Detre and White (1970) for unequal 'frames' (here, sediment 
        # volumes). The significance of d is based on the t distribution 
        # with an infinite degrees of freedom, which is the same as the 
        # cumulative normal distribution.
      }
    }  
  }
  
  # Remove peaks that do not pass the minimum-count screening-peak test
  for (j in 1:nThresholds) {
    insig.peaks <- intersect(which(Charcoal.charPeaks[ ,j] > 0),
                             which(CharThresh.minCountP[ ,j] > alphaPeak)) # Index for
                  # Charcoal.charPeaks values that also have p-value > alphaPeak...thus insignificant
    Charcoal.charPeaks[insig.peaks, j] <- 0 # set insignificant peaks to 0
    Charcoal.charPeaksThresh[insig.peaks, j] <- 0
  }
  
  # Calculate sensitivity indices
  
  
  cat('      ...done.')
  
#}
