#' Check for overlapping molecules.
#' 
#' Check for molecules overlapping in m/z, based on isotope fine structures from isopattern
#' or on centroids/intensoids from envelope.
#' 
#' @param pattern Output from isopattern or from envelope.
#' @param dmz m/z window. In combination with ppm=TRUE set as ppm or with ppm=FALSE set as absolute m/z.
#' @param ppm Should m/z window be set in ppm (TRUE) or absolute m/z (FALSE)?
#' 
#' @details
#' Overlaps in m/z among molecules are screened for within the m/z tolerance defined by the arguments dmz and ppm.
#' 
#' @return Dataframe with 4 columns, with number of rows equal to the length of argument pattern
#' \item{compound}{Chemical formula of the compound}
#' \item{warning}{Overlap detected?}
#' \item{to?}{If overlap: with wich other compound(s)? Refers to row number, recycled for peak_number.}
#' \item{peak_number}{If overlap: with which peak(s) of the other compound(s)? Refers to peak number.}
#' 
#' @author Martin Loos, Christian Gerber
#' 
#' @seealso isopattern, envelope
#' 
#' @examples
#' data(isotopes)
#' data(chemforms)
#' pattern <- isopattern(
#'   isotopes,
#'   chemforms,
#'   threshold=0.1,
#'   plotit=TRUE,
#'   charge=FALSE,
#'   emass=0.00054858,
#'   algo=1
#' )
#' check_several(pattern, dmz=0.001, ppm=FALSE)
#' 
#' @export
check_several <-
function(pattern,dmz,ppm=TRUE){

    ############################################################################
    # (1) issue warnings #######################################################
    if(length(pattern[[1]]) < 2 || !all(colnames(pattern[[1]])[1:2] == c("m/z", "abundance"))) stop("WARNING: pattern has invalid entries\n")
    if(ppm == TRUE & dmz < 0) stop("\n WARNING: ppm=TRUE -> dmz must be >0\n")
    if(ppm != TRUE & ppm != FALSE) stop("WARNING: ppm TRUE or FALSE")
    if(length(pattern) < 2) stop("WARNING: nothing to compare")
    ############################################################################
    # (2) check ################################################################    
    # (2.1) disassemble ######################################################## 
    getit1 <- getit2 <- getit3 <- c()
    for(i in 1:length(pattern)){
        getit1 <- c(getit1, rep(i,length(pattern[[i]][, 1])))
        getit2 <- c(getit2, pattern[[i]][, 1]);
        getit3 <- c(getit3, seq(1, length(pattern[[i]][, 1]), 1))
    }
    getit1 <- getit1[order(getit2)]
    getit3 <- getit3[order(getit2)]    
    getit2 <- getit2[order(getit2)]
    # (2.2) ####################################################################
    get1 <- get2 <- get3 <- rep("", length(pattern))
    for(i in 2:length(getit1)){
      if(!ppm){
        if(getit1[i - 1] != getit1[i]){
          if((getit2[i - 1] + dmz) >= getit2[i]){
             get1[getit1[i]] <- "TRUE"
             get1[getit1[i - 1]] <- "TRUE"
             get2[getit1[i]] <- paste0(get2[getit1[i]], getit1[i - 1], "/")
             get2[getit1[i-1]] <- paste0(get2[getit1[i - 1]], getit1[i], "/")
             get3[getit1[i]] <- paste0(get3[getit1[i]], getit3[i - 1], "/")
             get3[getit1[i-1]] <- paste0(get3[getit1[i - 1]], getit3[i], "/")        
          }
        }
      }else{ # ppm==TRUE
        if(getit1[i - 1] != getit1[i]){
          if((getit2[i - 1] + (getit2[i - 1] * dmz / 1e6)) >= getit2[i]){
             get1[getit1[i]] <- "TRUE"
             get1[getit1[i - 1]] <- "TRUE"
             get2[getit1[i]] <- paste0(get2[getit1[i]], getit1[i-1], "/")
             get2[getit1[i-1]] <- paste0(get2[getit1[i - 1]], getit1[i], "/")
             get3[getit1[i]] <- paste0(get3[getit1[i]], getit3[i - 1], "/")
             get3[getit1[i-1]] <- paste0(get3[getit1[i - 1]], getit3[i], "/",)        
          }
        }
      }
    }
    if(any(get1 != "")) message("Overlaps detected!")
    checked <- data.frame(names(pattern), get1, get2, get3)
    names(checked) <- c("compound", "warning", "to?", "peak_number")
    return(checked)
    ############################################################################

}
