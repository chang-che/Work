Vit_hmm <- function(obs.vec, pi.mar, tran.m, mean.vec, sd.vec){
  N <- nrow(tran.m)
  Ti <- length(obs.vec)
  density.m <- matrix(rep(as.double(0), Ti*N), nrow = Ti, ncol = N)
  for (i in 1:Ti){
    density.m[i,] <- mapply(
    dnorm, 'x' = obs.vec[i], 'mean' = mean.vec, 'sd' = sd.vec,'log' = T, SIMPLIFY = T)
  }
  tran.m.log <- log(tran.m)
  delta_log <- matrix(rep(NA, N*Ti), ncol = N)
  delta_log[1,] <- log(pi.mar) + density.m[1,]
  #keep track of states
  bt <- rep(NA, Ti)
  delta_ind <- matrix(rep(NA, (Ti-1)*N), ncol = N)
  for (i in 2:Ti){
    delta_max <- sapply(1:N, function(x) tran.m[,x]+delta_log[i-1, x]+density.m[i,x])
    delta_ind[i-1, ] <- apply(delta_max, 1, which.max)
    delta_log[i,] <- apply(delta_max,1,max)
  }
  bt[Ti] <- which.max(delta_log[Ti, ])
  for (i in Ti-1:1) {
    bt[i] <- delta_ind[i,bt[i+1]]
  }
  return(bt)
}