#' Multiply a chemical formula
#' 
#' Multiply all atom numbers in a chemical formula by a factor
#' 
#' @param formula_in Chemical formula to be multiplied, vector of character strings
#' @param fact Factor to multiply with
#' 
#' @details
#' Useful for adduct calculations, check adducts.
#' Chemical formulas must conform to what is described in check_chemform.
#' 
#' @return Multiplied chemical formula, character string
#' 
#' @author Martin Loos
#' 
#' @seealso adducts
#' 
#' @examples
#' formula_in <- "C10[13]C2H10Cl10"
#' multiform(formula_in, 3)
#' 
#' @export
multiform <-
function(formula_in, fact){

	formula_out <- c()
	for(i in 1:length(formula_in)){
  
		##############################################################################
		# on first chemical formula
		formula_in[i] <- gsub("D(?![a-z])","[2]H", formula_in[i], perl = T)
		len <- nchar(formula_in[i])
		element1 <- number1 <- c()
		##############################################################################
		# on formula_in (deduct)
		j <- 1
		while( j <= len){
	
			if(substr(formula_in[i], j, j) == "["){
				b <- j
				while(!any(substr(formula_in[i], j, j) == "]")) j <- j + 1
				k <- j
				while(!any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
				m <- j - 1
				element1 <- c(element1, substr(formula_in[i], b, m))
			}
	  
			if(!any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))){
				k <- j
				while(!any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
				m <- j <- j - 1
				element1 <- c(element1, substr(formula_in[i], k, m))
			}
		  
			if(any(substr(formula_in[i],j,j) == c("0","1","2","3","4","5","6","7","8","9"))){
				k <- j
				while(any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
				m <- j <- j - 1
				number1 <- c(number1, as.numeric(substr(formula_in[i], k, m)))
			}
	  
			j <- j + 1
		}
		##############################################################################
		# multiply
		number1 <- fact * number1
		formula_fin <- ""
		for(i in 1:length(element1)) formula_fin <- paste0(formula_fin, element1[i], number1[i])
		formula_out <- c(formula_out, formula_fin)
	
	}
	
	return(formula_out)

}

