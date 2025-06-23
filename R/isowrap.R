#' Combined (batch) calculation of isotope pattern, envelope and centroids/intensoids/valleys on interpolated resolutions
#' 
#' Wrapper combining the functions getR, isopattern, envelope and vdetect.
#' Uses chemical formulas from check_chemform as argument.
#' 
#' @param isotopes Dataframe listing all relevant isotopes, such as isotopes.
#' @param checked Output dataframe from check_chemform with correct chemical formulas.
#' @param resmass For resolution interpolation: dataframe with two columns, resolution and mass; see getR. 
#'        Otherwise, set to FALSE and use argument resolution to utilize a single resolution.
#' @param resolution Single resolution value. Only used if argument resmass is set to FALSE.
#' @param nknots Number of knots, see getR. Ignored if argument resmass set to FALSE.
#' @param spar Smoothing parameter, see getR. Ignored if argument resmass set to FALSE.
#' @param threshold Abundance below which isotope peaks are omitted, see isopattern.
#' @param charge z in m/z, see isopattern.
#' @param emass Electrone mass. Only relevant if charge is not set to FALSE, see isopattern.
#' @param algo Which algorithm to use? Type 1 or 2. See details section in isopattern.
#' @param ppm Set stick discretization, see details section of envelope.
#' @param dmz Set stick discretization, see details section of envelope.
#' @param frac Set stick discretization, see details section of envelope.
#' @param env Peak shape function; either "Gaussian" or "CauchyLorentz".
#' @param detect Detection from "centroid", "intensoid" or "valley".
#' @param plotit Should results be plotted, TRUE/FALSE?
#' @param verbose Verbose, TRUE/FALSE?
#' 
#' @return List with length equal to length of list \code{profiles}, with equal names of list entries.
# Each entry in that list contains the centroids, intensoids or valley of the envelope in two columns:
#' \item{m/z}{m/z}
#' \item{abundance}{area(centroid) or abundance (intensoid, valley)}
#' 
#' @author Martin Loos, Christian Gerber
#' 
#' @examples
#' data(isotopes)
#' data(resolution_list)
#' data(chemforms)
#' chemforms <- chemforms[1:10]
#' 
#' checked <- check_chemform(
#'     isotopes,
#'     chemforms
#'   )
#' 
#' resmass <- resolution_list[[1]]
#' 
#' centro <- isowrap(
#'   isotopes,
#'   checked,
#'   resmass=resolution_list[[4]],
#'   resolution=FALSE,
#'   nknots=4,
#'   spar=0.2,
#'   threshold=0.1,
#'   charge=1,
#'   emass=0.00054858,
#'   algo=2,
#'   ppm=FALSE,
#'   dmz="get",
#'   frac=1/4,
#'   env="Gaussian",
#'   detect="centroid",
#'   plotit=TRUE
#' )
#' 
#' @seealso vdetect
#' 
#' @export
isowrap <-
function(
  isotopes,
  checked,
  resmass,
  resolution = FALSE,
  nknots = 6,
  spar = 0.2,
  threshold = 0.1,
  charge = 1,
  emass = 0.00054857990924,
  algo = 2,
  ppm = FALSE,
  dmz = "get",   # retrieve dm from R=m/dm
  frac = 1/4,
  env = "Gaussian",
  detect = "centroid",
  plotit = FALSE,
  verbose = TRUE
){
    
 
    ############################################################################
    # (1) issue warnings #######################################################
    #if(resmass==FALSE & resolution==FALSE){stop("WARNING: either resmass or resolution; one set to FALSE")}
    #if(resmass!=FALSE & resolution!=FALSE){stop("WARNING: either resmass or resolution; one set to FALSE")}
    if(any(checked[, 1])) stop("WARNING: checked with incorrect chemical detected!")
    ############################################################################
    # (2) interpolate resolutions ##############################################
    if(length(resmass) > 1){
		resolution <- getR(
			checked,
			resmass = resmass,
			nknots = nknots,
			spar = spar,
			plotit = plotit
		)
    }
    ############################################################################
    # (3) isotope pattern ######################################################
    pattern <- isopattern(
		isotopes,
		checked[, 2],
		threshold = threshold,
		charge = charge,
		emass = emass,
		plotit = plotit,
		algo = algo,
		verbose = verbose
    )
    ############################################################################
    # (4) profiles #############################################################
    profiles <- envelope(
        pattern,
        ppm = ppm,
        dmz = dmz, 
        frac = frac,
        env = env,
		resolution = resolution,
		plotit = plotit,
		verbose = verbose
    )
    ############################################################################
    # (5) centroidization ######################################################
    centro <- vdetect(
		profiles,
		detect = detect,
		plotit = plotit,
		verbose = verbose
    )
    ############################################################################
    return(centro)
    ############################################################################

}
