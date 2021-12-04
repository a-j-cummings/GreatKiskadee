# performs a single update via the Metropolis algorithm
metropolis_update <- function(x0, log_density, w){
  x1 <- x0 + rnorm(1,0,w)
  ld_x0 <- log_density(x0)
  ld_x1 <- log_density(x1)
  log_r <- ld_x1 - ld_x0 # log metropolis ratio
  if (runif(1) < exp(log_r)){
    return(x1)
  } else {
    return(x0)
  }
}
