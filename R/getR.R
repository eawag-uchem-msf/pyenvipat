#' Interpolation of MS measurement resolution
#' 
#' Given a set of MS measurement resolutions (R) as a function of measurement mass (m/z), getR interpolates R for any given molecular mass(es)
#' calculated by check_chemform using smooth.spline.
#' 
#' @param checked Dataframe produced by check_chemform.
#' @param resmass Dataframe with two columns, resolution and mass; such as the list entries in resolution_list.
#' @param nknots Integer number of knots to use for the smoothing spline. Default = 13. See also smooth.spline.
#' @param spar Smoothing parameter, (0,1]. See also smooth.spline.
#' @param plotit Plot results, TRUE/FALSE?
#' 
#' @note
#' check_chemform gives molecular masses (m/z) for z=+/-1 only. 
#' If z>1 or z<-1 is required, molecular mass entries in argument checked have to be divided accordingly to be consistent.
#' 
#' @return Vector with resolutions.
#' 
#' @author Martin Loos, Christian Gerber
#' 
#' @seealso smooth.spline, check_chemform, resolution_list
#' 
#' @examples
#' data(resolution_list)
#' resmass <- resolution_list[[4]]
#' data(isotopes)
#' data(chemforms)
#' checked <- check_chemform(isotopes, chemforms)
#' resolution <- getR(checked, resmass, nknots=13, spar=0.1, plotit=TRUE)
#' 
#' # same for z=-2:
#' checked <- check_chemform(isotopes, chemforms)
#' checked[,3] <- (checked[,3]/abs(-2))
#' resolution <- getR(checked, resmass, nknots=13, spar=0.1, plotit=TRUE)
#' 
#' @export
getR <-
function(checked, resmass, nknots = 13, spar = 0.1, plotit = TRUE){

    ############################################################################
    # (1) issue warnings #######################################################
	if(length(resmass[, 1]) < 10) stop("stop: not enough data points in resmass\n")
    if(any(checked[, 3] < min(resmass[, 1])) || any(checked[, 3] > max(resmass[, 1]))) stop("stop: some mean_mass out of range of resmass\n")
    if(nknots < 3 || nknots > length(resmass[,1])) stop("stop: invalid nknots\n")
    if(spar <= 0 || spar >= 1) stop("stop: invalid spar; spar=(0,1]")
    if(!is.logical(plotit)) stop("stop: plotit invalid. TRUE, FALSE.\n")
    if(any(checked[, 1] == TRUE)) stop("stop: in checked[,1] ... invalid chemical formula!")
    ############################################################################
    # (2) fit gam model ########################################################
    model <- smooth.spline(resmass[, 1], resmass[, 2], cv = TRUE, all.knots = FALSE, nknots = nknots, spar = spar)
    if(plotit){
		plot(model, pch=19, type = "l", lwd = 2, col = "red", xlab = "mass", ylab = "Resolution")
		points(resmass[, 1], resmass[, 2], pch = 19, cex = 0.7);
    }
    ############################################################################
    # (3) predict ##############################################################
    resolution <- c()
    for(i in 1:length(checked[, 1])){
		it <- predict(model,checked[i, 3])
		if(plotit) points(checked[i, 3], it$y, pch = 19, cex = 0.7, col = "green")
		resolution <- c(resolution,it$y)  
    }
    return(resolution)
    ############################################################################

}
