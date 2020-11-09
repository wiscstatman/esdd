dpmm_beta <-
function(x0, alpha, beta, m0, burn_in = 200, sample_size = 10, thinning = 5){
#  source("Initial_beta.R")
#  source("Update_beta.R")
  
  ## MCMC for random clustering
  m = nrow(x0)
  n = ncol(x0)
  KK = rep(0,sample_size)
  NN = matrix(0, sample_size, 2*m)
  CC = matrix(0,sample_size, m)
  ### Initialization
  cl = Initial_beta(x0, m0)
  ### burn-in stage
  for (i in 1:burn_in) {
    cl = Update_beta(cl, x0, alpha, beta, m0)
  }
    
  ### Sample clustering assignments
  for (i in 1:sample_size) {
    tmp = 0
    while (tmp < thinning) {
      cl = Update_beta(cl, x0, alpha, beta, m0)
      tmp = tmp + 1
    }
    KK[i] = cl$K
    NN[i, ] = cl$N
    CC[i, ] = cl$C
  }
  
  cl_sample = list(KK = KK, NN = NN, CC = CC)
  return(cl_sample)
}
