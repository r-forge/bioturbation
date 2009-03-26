### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### slice midpoints from slice limits
### ============================================================================

midpoints <- function(limits){
  midpoints <- NULL
  for(point in 1:(length(limits)-1)){
    midpoints <- c(midpoints,(limits[point]+limits[point+1])/2)
  }
  return(midpoints)
}

