

sigmoid = function(x) 1 / (1 + exp(-x))

alrinv <- function(y) {
  x <- cbind(exp(y), 1)
  x / rowSums(x)
}

getBetaFromCodacore = function(SBP_codacore, coeffs_codacore, p){
  if(!is.matrix(SBP_codacore)) SBP_codacore = matrix(SBP_codacore)
  if(ncol(SBP_codacore) != length(coeffs_codacore)){
    stop("getBetaFromCodacore: SBP and coefficients don't match")
  }
  kplus = apply(SBP_codacore, 2, function(col) sum(col == 1))
  kminus = apply(SBP_codacore, 2, function(col) sum(col == -1))
  reciprocals = matrix(
    0, nrow = nrow(SBP_codacore), ncol = ncol(SBP_codacore))
  for(i in 1:ncol(SBP_codacore)){
    reciprocals[SBP_codacore[, i] == 1, i] = 1 / kplus[i]
    reciprocals[SBP_codacore[, i] == -1, i] = -1 / kminus[i]
  }
  ReciprocalstimesCoeffs = matrix(
    NA, nrow = nrow(SBP_codacore), ncol = ncol(SBP_codacore))
  for(i in 1:ncol(ReciprocalstimesCoeffs)){
    ReciprocalstimesCoeffs[, i] = reciprocals[, i] * coeffs_codacore[i]
  }
  beta = rowSums(ReciprocalstimesCoeffs)
  return(beta)
}