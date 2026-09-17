#' enviPat: Calculation of isotope patterns, stick profiles (envelopes) and centroids/intensoids for mass spectrometry
#' 
#' Fast and memory-efficient calculation of isotope patterns (fine structures) for up to very large molecules, based on three different algorithms.
#' Subsequent convolution of isotope patterns with a peak shape function to theoretical envelopes (profiles).
#' Based on envelopes, valley detection and centroidization/intensoid calculation.
#' Allows for batch processing of chemical formulas and interpolation of measurement resolutions.
#' Includes a wrapper combining all of the above functionalities.
#' 
#' Furthermore, includes:
#' \itemize{
#'   \item a check for consistency of chemical formulas
#'   \item a check for molecules with overlapping isotope patterns
#'   \item a list of all stable isotopes
#'   \item a list of different resolution data sets for Thermo Orbitrap and QExactive high-resolution mass spectrometers
#'   \item a list of adducts formed during electorspray ionization (ESI)
#' }
#' 
#' A web-based GUI for enviPat is freely available under \url{https://www.envipat.eawag.ch}.
#' 
#' @name enviPat-package
#' @aliases enviPat
#' 
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://www.envipat.eawag.ch/}
#'   \item \url{https://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00941}
#'   \item Report bugs at \url{https://github.com/blosloos/enviPat/issues}
#' }
NULL

.onAttach <- function(lib, pkg)
{
	packageStartupMessage(paste0("\n \n Welcome to enviPat version 2.8 \n Check https://www.envipat.eawag.ch for an interactive online version\n"));

}
