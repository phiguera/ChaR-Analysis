# 3. Smooth Charcoal.I to estimate Low-frequency trends (i.e. Char.background) ####


CHARsm <- function(x, smoothing.yr=1000) {

  # Extract data from 'x'
  accI <- x$accI
    
  # Get parameters needed
  yr.interp <- x$yrInterp
  n.smooth <- round(smoothing.yr/yr.interp) # number of data points over which to smooth the record.
  # Value rounded when used in smooth functions. = 's' in Matlab version
  # Proportion of datapoints used for Loess and Lowess smoothing functions:
  span <- n.smooth/length(accI) 
  
  
  # Prepare empty data.frame
  charAccIS <- data.frame(matrix(NA, nrow=length(x$cmI), ncol=6))
  colnames(charAccIS) <- c("Lowess", "RobLowess", "Loess", "MovAver", "MovMed", "MovMode")
  
  # # Lowess (Method #1 in Matlab version)
  charAccIS[ ,1] <- lowess(x=x$accI, f=span, iter=0)$y
  
  
  # Robust Lowess (Method #2 in Matlab version)
  charAccIS[ ,2] <- lowess(x=x$accI, f=span, iter=4)$y
  
  
  # Loess with default options
  in.loess <- data.frame(x$ybpI, x$accI)
  charAccIS[ ,3] <- loess(formula = x.accI ~ x.ybpI, data = in.loess,
                          span = span)$fitted
  
  rm(in.loess)
  
  
  # Moving average (Method #3 in Matlab version)
  charAccIS[ ,4] <- runmean(x=x$accI, k=n.smooth, alg="exact", endrule="mean",
                            align="center")
  
  
  # Moving median (Method #4 in Matlab version; translated from Matlab version)
  for (i in 1:length(x$accI)) {
    if (i <= round(n.smooth/2)) { # if 1/2 n.smooth twords start
      CHARi_t <- x$accI[1:round(n.smooth)] # Charcoal.accI for year t
      charAccIS[i,5] <- median(CHARi_t)
    }
    if (i >= length(x$accI)-round(n.smooth)) { # if 1/2 s twords end
      CHARi_t <- x$accI[length(x$accI)-round(n.smooth/2):
                                   length(x$accI)]
      charAccIS[i,5] <- median(CHARi_t)
    }
    if (i > round(n.smooth/2) && i < length(x$accI)-round(n.smooth)) { # else, you're in the middle of the record
      CHARi_t <- x$accI[round(i-0.5*n.smooth):round(i+0.5*n.smooth)]
      charAccIS[i,5] <- median(CHARi_t)
    }
  }
  charAccIS[ ,5] <- lowess(x=charAccIS[ ,5], f=span, iter=0)$y
  # [# cm^-2 yr^-1] smoothed with lowess filter 
  
  
  # Running mode (Method #5 in Matlab version; translated from Matlab version)
  bin <- 100  # bins to divide x$accI into
  for (i in 1:length(x$accI)) {
    if (i <= round(n.smooth/2)) {                        # if 1/2 s twords start
      CHARi_t <- x$accI[1:round(n.smooth)]     # Charcoal.accI for year t
      mode_bin <- c(0, diff(range(CHARi_t))/bin, max(CHARi_t))
      n.x <- data.frame(n=hist(CHARi_t, bin, plot=F)$c, x=hist(CHARi_t, bin, plot=F)$mids)
      mode_in <- n.x$x[which(n.x$n == max(n.x$n))]
      charAccIS[i,6] <- median(mode_in)
    }
    if (i >= length(x$accI)-round(n.smooth)) {
      CHARi_t <- x$accI[length(x$accI)-round(n.smooth/2):
                                   length(x$accI)]
      mode_bin <- c(0, diff(range(CHARi_t))/bin, max(CHARi_t))
      n.x <- data.frame(n=hist(CHARi_t, bin, plot=F)$c, x=hist(CHARi_t, bin, plot=F)$mids)
      mode_in <- n.x$x[which(n.x$n == max(n.x$n))]
      charAccIS[i,6] <- median(mode_in)
    }
    if (i > round(n.smooth/2) && i < length(x$accI)-round(n.smooth)) {
      CHARi_t <- x$accI[round(i-0.5*n.smooth):round(i+0.5*n.smooth)]
      mode_bin <- c(0, diff(range(CHARi_t))/bin, max(CHARi_t))
      n.x <- data.frame(n=hist(CHARi_t, bin, plot=F)$c, x=hist(CHARi_t, bin, plot=F)$mids)
      mode_in <- n.x$x[which(n.x$n == max(n.x$n))]
      charAccIS[i,6] <- median(mode_in)
    }
  }
  charAccIS[ ,6] <- lowess(x=charAccIS[ ,6], f=span, iter=0)$y
  # [# cm^-2 yr^-1] smoothed with lowess filter
  
  # end running mode smoother
  
  # Clean Environment
  rm(bin, CHARi_t, mode_bin, n.x, mode_in)
  
  # Write output to Environment
  # output <- structure(list(Charcoal.I=x, smooth=charAccIS, smoothing.yr=smoothing.yr))
  # class(output) <- "CHARsm"
  # return(output)
  # 
  output <- structure(list(cmI=x$cmI, ybpI=x$ybpI, countI=x$countI, volI=x$volI, conI=x$conI,
                           accI=x$accI, ageTop=x$ageTop, ageBot=x$ageBot, yrInterp=x$yrInterp,
                           acc=x$acc, smLowess=charAccIS[ ,1], smRobLowess=charAccIS[ ,2], smLoess=charAccIS[ ,3],
                           smMovAver=charAccIS[ ,4], smMovMed=charAccIS[ ,5], smMovMode=charAccIS[ ,6],
                           smoothing.yr=smoothing.yr))
  class(output) <- "CHARsm"
  return(output)
}




# Plot

plot.CHARsm <- function(sm, xlim=NULL, ylim=NULL,xlab=NULL,ylab=NULL,...) {

  # Get data
  ageTop <- sm$ageTop
  acc <- sm$acc
  ageI <- sm$ybpI
  accI <- sm$accI
  sm.low <- sm$smLowess
  sm.roblow <- sm$smRobLowess
  sm.loe <- sm$smLoess
  sm.mav <- sm$smMovAver
  sm.med <- sm$smMovMed
  sm.mod <- sm$smMovMode

  # Set axis limits
  if(is.null(xlim)){
  x.min <- floor(min(ageI)/100) * 100
  x.max <- ceiling(max(ageI)/100) * 100
  if(x.max - x.min > 1000) {
    x.by <- 1000
  } else {
    x.by <- round((x.max - x.min)/100) * 100
  }
  x.lim <- c(x.max, x.min)
  }
  if(is.null(ylim)) {
  y.lim <- c(min(accI), 1.2*max(accI))
  }

  par(mfrow=c(2,1), mar=c(0.5, 4, 0.5, 1), oma=c(5,1,1,1), cex=0.7)
  plot(ageTop, acc, type="h", col=grey(0.7), xlim=x.lim, ylim=y.lim,
       xlab="time (cal yr. BP)", ylab="CHAR (# cm^-2 yr^-1)", axes=F)
  par(new=T)
  plot(ageI, accI, type="s", lwd=0.5, xlim=x.lim, ylim=y.lim, ylab="", axes=F)
  axis(1, at=seq(0, x.max, by=1000), labels=F)
  axis(2)
  legend("topleft", inset=c(0,0.2),
         legend=c("CHAR raw", "CHAR interpolated"),
         lwd=1.5, col=c("grey", "black"),
         border="NA", bty="n", horiz=F, cex=0.7)

  plot(ageI, accI, type="s", xlim=x.lim, ylim=y.lim,
       xlab="", ylab="CHAR (# cm^-2 yr^-1)", lwd=1, axes=F)
  polygon(x=c(rev(ageI), ageI),
          y=c(rep(0, length(accI)), accI), col="gray50", border=T, lwd=0.5)
  lines(ageI, sm.low, type="l", col="blue", lwd=1.5)
  lines(ageI, sm.roblow, type="l", col="forestgreen", lwd=1.5)
  lines(ageI, sm.loe, type="l", col="green", lwd=1)
  lines(ageI, sm.mav, type="l", col="red", lwd=1.5)
  lines(ageI, sm.med, type="l", col="cyan3", lwd=1.5)
  lines(ageI, sm.mod, type="l", col="orange", lwd=1.5)
  legend("topleft", inset=c(0,0.2),
         legend=c("CHARinterpolated", "Lowess","Robust Lowess","Loess","Moving average","Moving median","Moving Mode"),
         lwd=1.5,
         col=c("gray50","blue","forestgreen","green","red","cyan3","orange"),
         border="NA", bty="n", horiz=F, cex=0.7)
  axis(1, at=seq(0, x.max, by=1000))
  axis(2)
  mtext(paste("b) CHARinterpolated and different options for a ", sm$smoothing.yr,"yr CHARbackground",
              sep=''), cex=0.8, line=-2)
}