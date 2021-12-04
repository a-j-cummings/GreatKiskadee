# Slice sampling
single_slice_sample <- function(x0, log_density, w){
  ld_x0 <- log_density(x0)
  # take slice
  y <- log(runif(1, 0, exp(ld_x0)))
  # set window
  x_l <- x0 - runif(1)*w
  x_u <- x_l + w
  # stepping out
  while (y < log_density(x_l)){
    x_l <- x_l - w
  }
  while (y < log_density(x_u)){
    x_u <- x_u + w
  }
  # sample/stepping in
  x1 <- runif(1, x_l, x_u)
  ld_x1 <- log_density(x1)
  while (ld_x1 < y){
    if (x1 < x0){
      x_l <- x1
    } else {
      x_u <- x1
    }
    x1 <- runif(1, x_l, x_u)
    ld_x1 <- log_density(x1)
  }
  return(x1)
}
