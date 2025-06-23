#' Resolution list
#'
#' List of different resolutions R=f(m/z) for various (high-resolution) mass spectrometers.
#' For each of the instruments, different resolution settings are available.    
#' Here, R is defined as R=(m/z)/(dm/z), with dm/z = peak width at half maximum.
#' Serves as input to \link[enviPat]{getR} to interpolate R from given molecular masses.
#'
#' @format A list where each element is a data frame containing two columns:
#' \describe{
#'   \item{mass}{The m/z values}
#'   \item{resolution}{The corresponding resolution values at each m/z}
#' }
#'
#' @examples
#' data(resolution_list)
#' resolution_list[[1]]  # Access the first resolution dataset
#' @keywords datasets
"resolution_list"
