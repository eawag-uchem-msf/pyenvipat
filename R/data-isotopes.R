#' Isotope data
#'
#' A dataset containing stable isotope information for various elements, including their
#' masses and natural abundances.
#'
#' @format A data frame containing isotope information with the following variables:
#' \describe{
#'   \item{element}{Element symbol}
#'   \item{isotope}{Isotope name or designation}
#'   \item{mass}{Exact mass of the isotope}
#'   \item{abundance}{Natural abundance of the isotope (fraction)}
#'   \item{\code{ratioC}}{Maximum number of atoms of an element for one C-atom in a molecule, based on 99.99 \% of case molecules.}
#' }
#'
#' @details
#' The \code{ratioC}-value stems from a database survey conducted by Kind&Fiehn (2007); to disable, set value to 0.
#' The list serves as input into several package nontarget-functions. The first column of the data frame also
#' contains names of specific isotopes used for labeled compounds.
#' 
#' @source \url{https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl}
#'
#' @references Kind, T. and Fiehn, O., 2007. Seven golden rules for heuristic filtering  of molecular formulas obtained by
#' accurate mass spectrometry. BMC Bioinformatics, 8:105.
#'
#' @examples
#' data(isotopes)
#' isotopes
#' @keywords datasets
"isotopes"
