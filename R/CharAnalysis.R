#' @title ChaR-Analysis
#' @description Diagnostic and anlystical tools for peak analysis in 
#' sediment-charcoal records. This code is based on original Matlab code, 
#' CharAnalysis reads in a data object or a filename and uses it to 
#' parameterize the detection of charcoal peaks from background values.
#'
#' @importFrom RJSONIO fromJSON
#' @param x Optional parameter for a \code{site}, \code{dataset}, or \code{dataset_list}.
#' @author Philip Higuera \email{phiguera@@uidaho.edu}
#' @return This command returns either object of class \code{"try-error"}' (see \code{\link{try}}) definined by the error returned from the Neotoma API call, or an object of class \code{charanalysis}, containing a set of objects:
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

CharAnalysis <- function(x){                      }

require(paleofire)


# 1. Load data
load("data/CO_charData.RData")
