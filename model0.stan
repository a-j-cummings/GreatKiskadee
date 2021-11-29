// data block
data {
  int<lower=0> n_classes;
  int<lower=0> n_copies;
  int<lower=0> n_ind;
  int<lower=0> n_loc;
  int<lower=0, upper=n_ind> ind[n_ind];
  int<lower=0, upper=n_copies> copy[n_ind*n_copies];
  int<lower=0, upper=1> X[n_ind, n_loc, n_copies];
  real<lower=0> alpha_q[n_classes];
  real<lower=0> alpha_p[2];
}

// parameters block 
parameters {
   real<lower=0, upper=1> Z[n_ind, n_loc, n_copies, n_classes];
   matrix[n_ind, n_classes] Q;
   real<lower=0> P[n_classes, n_loc, 2];
}


// model block
model {
  // prior model
  for (i in 1:n_ind){
    Q[i,] ~ dirichlet(alpha_q);
    for (m in 1:n_loc){
      for (c in 1:n_copies){
        Z[i,m,c,] ~ multinomial(Q);
      }
    }
  }
  for (j in 1:n_classes){
    for (m in 1:n_loc){
      P[j,m,] ~ dirichlet(alpha_p)
    }
  }
  // data model
}

