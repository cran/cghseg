getmBIC <- function(K,lv,mu,CGHo){  
  M = length(names(mu))
  N = M*mu[[1]]$end[dim(mu[[1]])[1]]
  Ent = sum(unlist(lapply(mu,FUN = function(x){log(x$end-x$begin+1)})))
 # mBIC = 0.5*(N-K+1)*(2*lv/N + 1 + log(2*pi)-log(N))-0.5*Ent-(K-M)*log(N)

  if (CGHo["calling"]==FALSE){
    mBIC = ((N-K+1)/2)*(lv*(2/N)+1+log(2*pi))-0.5*Ent -(K-M)*log(N)+lgamma((N-K+1)/2)-((N-K+1)/2)*log(N)
#    mBIC = 0.5*(N-K+1)*(2*lv/N + 1 + log(2*pi)-log(N))-0.5*Ent-(K-M)*log(N)
  } else {
    P  = CGHo["nblevels"]
    Np = lapply(mu,FUN=function(x){
      nk = x$end-x$begin+1
      rep(x$levels,nk)
    })
    Np   = tabulate(unlist(Np))
    Ent  = sum(log(Np))
    mBIC = ((N-P+1)/2)*(lv*(2/N)+1+log(2*pi))-0.5*Ent-(K-M)*log(N)+lgamma((N-P+1)/2)-((N-P+1)/2)*log(N)
#    mBIC = 0.5*(N-P+1)*(2*lv/N + 1 + log(2*pi)-log(N))-0.5*Ent-(K-M)*log(N)
  }
  return(mBIC)
}
