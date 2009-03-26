### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### luminophore concentrations from luminophore numbers and slice limits
### ============================================================================


numtoconc <- function(values, limits){
  numtoconc <- NULL
  for(slice in 1:length(values)){
    numtoconc <- c(numtoconc, values[slice] / (limits[slice+1] - limits[slice]))
  }
  return(numtoconc)
}