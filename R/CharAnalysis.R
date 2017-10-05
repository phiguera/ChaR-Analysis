#' @title ChaR-Analysis
#' @description Diagnostic and anlystical tools for peak analysis in 
#' sediment-charcoal records. This code is based on original Matlab code, 
#' CharAnalysis reads in a data object or a filename and uses it to 
#' parameterize the detection of charcoal peaks from background values.
#'
#'  \item{ \code{charcoal} }{}
#'  \item{ \code{char.thresh} }{}}
#'  \item{ \code{pre.treatment} }{}
#'  \item{ \code{smoothing} }{}
#'  \item{ \code{peak.analysis} }{}
#'  \item{ \code{results} }{}
#'
#' @section Note:
#'
#' @examples \dontrun{
#'
#' }
#' @references
#' CharAnalysis documentation: http://code.google.com/p/charanalysis/
#' @keywords IO connection
#' @export

CharAnalysis <- function(site.name="CO", runname=NULL)
{  

# Parameters that should go into main function, used to test here the code...
  # SHOULD BE DELETED ONCE EVERYTHING IS READY
#setwd('/Users/wfinsing/Documents/GitHub/ChaR-Analysis/Cores')
#site.name <- "CO"
#runname <- NULL
#runname <- "1"

#### Load packages
require(paleofire)

#### Determine input directory
input.dir <- file.path("..", site.name)

#### Create output directory
if (is.null(runname)) {
output.dir <- file.path(".", "Cores", site.name, "output")
output.dir <- file.path(".", site.name, "output")
} else {
  output.dir <- file.path(".", site.name, paste0("output", runname))
}

if (!dir.exists(output.dir)) {
  dir.create(output.dir)
} else {
  warning("Note: The output directory already exists and was overwritten with new output!")
}

# Load Charcoal data
cat(' (1) Reading input files...')
Charcoal <- read.csv(file.path(".", site.name, paste0(site.name, "_charData.csv")))
char.series <- Charcoal[ , 6]
char.params <- Charcoal[ , 1:5]

# Load Parameters file
Params   <- read.csv(file.path(".", site.name, paste0(site.name, "_charParams.csv")),
                     header=T,
                     colClasses = c("NULL", "factor", "numeric", "NULL", "NULL"))

# Extract data from Parameters file
zones <- Params[1:9, ]
zones <- na.omit(zones)

zones <- zones[which(complete.cases(zones)), ]
zones <- zones[ ,2]

yr.interp     <- if(Params[10, 2] == 0) {
  yr.interp = NULL
}
char.tr       <- Params[11, 2]
char.tr.meth  <- Params[12, 2]
char.smooth   <- Params[13, 2]
cPeak         <- Params[14, 2]
thresh.type   <- Params[15, 2]
thresh.meth   <- Params[16, 2]
thresh.values <- Params[17:20, 2]
minCountP     <- Params[21, 2]
peakFrequ     <- Params[22, 2]


# 2. Pretreatment
cat('\n (2) Pretreating charcoal data...')
Charcoal.I <- pretreatment(params = char.params, serie = char.series, Int = T,
                           first <- zones[1], last <- zones[length(zones)],
                           yrInterp = yr.interp)

# 3. Smooth Charcoal.I to estimate Low-frequency trends (i.e. Char.background)
cat('\n (3) Smoothing resampled CHAR to estimate low-frequency trends...')
cat('\n     and calculating peak CHAR')


}