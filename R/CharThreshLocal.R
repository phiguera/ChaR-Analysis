## Define possible threshold for peak identification ####
# Based on PH 'CharThreshLocal.m'
# Determines a threshold value for each interpolated sample, based on the
# distribution of CHAR values within the selected window (yr) and either
# a Gaussian mixture model or the assumption that the noise component 
# of the peak charcoal record (C_peak) is normally distributed 
# around 0 (since C_peak is defined by residuals).

# NB: Threshold is defined locally...

# Reads output from CHARsm() function.



CharThreshLocal <- function(x, sm.meth=NULL, thresh.yr=NULL,
                            thresh.values=c(0.95, 0.99,0.999,0.99), plot.thresh=T) {
  
  if(is.null(sm.meth)) sm.meth <- 2
  if(is.null(thresh.yr)) thresh.yr <- x$smoothing.yr
  if(sm.meth == "Lowess") accIS.n <- 1
  if(sm.meth == "Robust Lowess") accIS.n <- 2
  if(sm.meth == "Loess") accIS.n <- 3
  if(sm.meth == "Moving Average") accIS.n <- 4
  if(sm.meth == "Moving Median") accIS.n <- 5
  if(sm.meth == "Moving Mode")  accIS.n <- 6
  
  accI <- x$accI
  accIS <- x[[10+accIS.n]]
  yr.interp <- x$yrInterp
  
  n.smooth <- round(x$smoothing.yr/yr.interp)
  # Proportion of datapoints used for Loess and Lowess smoothing functions:
  span <- n.smooth/length(accI) 
  Charcoal.I <- x$Charcoal.I
  
  # Calculate peak series
  Charcoal.peak <- as.matrix(accI - accIS)
  
  # Create space for local variables
  thresh.pos <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak),
                                      ncol=length(thresh.values))) # space for threshold values
  thresh.neg <- data.frame(matrix(data=NA, nrow=length(Charcoal.peak),
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
    #cat(paste0("Calculating ", i, "th local threshold of ", length(Charcoal.peak)))
    
    if (i < round(0.5*(thresh.yr/yr.interp))+1) { # First 'thresh.yr' samples.
      #         X = Charcoal.peak(1:round(0.5*(thresh.yr/r))); % Pre June 2009.
      #         X = Charcoal.peak(1:round(thresh.yr/r)); % Modified, June 2009, PEH.
      X <- Charcoal.peak[1:round(0.5*(thresh.yr/yr.interp))+i] # Modified, % June 2009, PEH.
    }
    if (i > (length(Charcoal.peak)-round(0.5*(thresh.yr/yr.interp)))) {  # Last 'thresh.yr' samples.
      #             X = Charcoal.peak(length(Charcoal.peak)-...
      #                 round((thresh.yr/r)):end);   % Pre June 2009.
      X <- Charcoal.peak[(i-round(0.5*(thresh.yr/yr.interp))):length(Charcoal.peak)]   # Modified, June 2009, PEH.
      # As recommended by RK, this uses samples from 
      # a half-window before i, all the way to end of record.
    }
    if (i >= round(0.5*(thresh.yr/yr.interp))+1 && i <= (length(Charcoal.peak)-round(0.5*(thresh.yr/yr.interp)))) {
      # All samples between first and last 'thrshYr' samples.
      X <- Charcoal.peak[(i-round(0.5*(thresh.yr/yr.interp))) : (i+round(0.5*(thresh.yr/yr.interp)))]
    }
    
    
    ## ESTIMATE LOCAL NOISE DISTRIBUTION
    # Estimate noise distribution with Guassian mixture model
    if (sum(X) == 0) {
      cat("NOTE: All C_peak values = 0; cannot fit noise distribution.")
      cat("\n      Mean and standard deviation forced to equal 0.") 
      cat("\n      Consider longer smoothing window or alternative for") 
      cat("\n      threshMethod parameter.")
      muHat[i, ] <- 0
      sigmaHat[i, ] <- 10^-100
      propN[i, ] <- 0
    } else {
      m <- densityMclust(data=X, G=2, verbose=F)
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
    
    rm(m)
    
    ## Define local threshold, SNI, and GOF
    # Define range of threshold values, plus threshold value selected.
    thresh.pos[i, ] <- qnorm(p=thresh.values, mean=muHat[i,1], sd=sigmaHat[i,1])
    thresh.neg[i, ] <- qnorm(p=1-thresh.values, mean=muHat[i,1], sd=sigmaHat[i,1])
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
    # if (length(noise_i) > 3) {
    #   ksP <- ks.test(x=noise_i, y="pnorm", mean=muHat[i,1], sd=sigmaHat[i,1])$p.value
    #   CharThresh.GOF[i, ] <- ksP
    # }
    # colnames(CharThresh.GOF) [1] <- "GOF"
    
    
    
    # PLOT SELECTED Charcoal.peak DISTRIBUTIONS
    if (any(i == num.plots)) {
      par(mfrow=c(1,1))
      h <- hist(x=X, breaks=50, plot=F)
      dens <- hist(X, breaks=50, plot=F)$density
      pdf1 <- dnorm(x=dens, mean=muHat[i, 1], sd=sigmaHat[i, 1])
      
      plot(h, freq=F, col="grey", border="grey", xlim=c(min(h$breaks), max(h$breaks)),
           ylab="Density", xlab='',
           main=paste(x$ybpI[i], "yr BP"))
      par(new=T)
      pdf1 <- curve(dnorm(x=x, mean=muHat[i, 1], sd=sigmaHat[i, 1]),
                    from=min(h$breaks), to=max(h$breaks),
                    ylim=c(0, max(dens)), type="l", col="blue", lwd=1.5, axes=F, ylab='', xlab='')
      par(new=T)
      pdf2 <- curve(dnorm(x=x, mean=muHat[i, 2], sd=sigmaHat[i, 2]),
                    from=min(h$breaks), to=max(h$breaks),
                    ylim=c(0, max(dens)), type="l", col="orange", lwd=1.5, axes=F, ylab='', xlab='')
      par(new=T)
      lines(x=c(qnorm(p=thresh.values[4], mean=muHat[i,1], sd=sigmaHat[i,1]),
                qnorm(p=thresh.values[4], mean=muHat[i,1], sd=sigmaHat[i,1])), y=c(0, max(dens)),
            type="l", col="red", lwd=1.5)
      mtext(text=paste0("SNIi= ", round(CharThresh.SNI[i, ], digits=2),
                        "\nKS p-val= ", round(CharThresh.GOF[i, ], digits=2)),
            side=3, las=0, line=-2)
      
      my.plots[[j]] <- recordPlot()
      j<- j+1
    }
    
  } # end loop for each Charcoal.peak
  
  rm(h, dens, pdf1, pdf2, noise_i, sig_i)
  
  
  
  # Print pdf with selected plots that were saved at the end of the loop above
  if (plot.thresh == T) {
    pdf('Threshold_determination.pdf', onefile=TRUE, paper="a4")
    par(mfrow=(c(5,5)), mar=c(0.5,1,0.5,1), oma=c(1,1,0.5,1), cex=0.7)
    for (k in 1:length(num.plots)) {
      replayPlot(my.plots[[k]])
    }
    dev.off()
  }
  rm(my.plots, num.plots)
  
  
  
  ## Smooth thresholds with Lowess smoother
  CharThresh.SNI$sSNI <- lowess(x=CharThresh.SNI[ ,1], f=span)$y
  # CharThresh.SNI (CharThresh.SNI < 0,1) = 0; # Matlab code...what does this mean?
  for (i in 1:length(thresh.values)) {
    thresh.pos[ ,i] <- lowess(x=thresh.pos[ ,i], f=span)$y
    thresh.neg[ ,i] <- lowess(x=thresh.neg[ ,i], f=span)$y
  }
  
  x <- unclass(x)
  output <- structure(list(cmI=x$cmI, ybpI=x$ybpI, countI=x$countI, volI=x$volI, conI=x$conI,
                          accI=x$accI, charBack=accIS, ageTop=x$ageTop, ageBot=x$ageBot, yrInterp=x$yrInterp,
                          acc=x$acc, Charcoal.peak=Charcoal.peak, SNI=CharThresh.SNI$SNI,
                          SNIsm=CharThresh.SNI$sSNI, thresh.pos=thresh.pos, thresh.neg=thresh.neg,
                          thresh.values=thresh.values))
  class(output) <- "CharThresh"
  return(output)
  
}
