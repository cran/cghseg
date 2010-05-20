unisegclust <- function(Y,CGHo,Kmax){

  P            = CGHo["nblevels"]
  loglik       = c(-Inf,Kmax)
  param        = list()  
  n.com        = length(Y)
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  lmax         = floor(CGHo["lmax"]*length(x))  
  out.hybrid   = hybrid(x,P,Kmax,CGHo["lmin"],lmax,vh = TRUE,CGHo["fast"])
  param        = out.hybrid$param
  loglik       = c(out.hybrid$Linc)
  
  if (CGHo["select"]=="none"){
    Kselect = Kmax
  } else {    
    Kseq    = c(P:Kmax)
    Kselect = PKselection(Y,Kseq,param,loglik,CGHo)
  }
  phi        = param[[Kselect]]$phi
  rupt       = param[[Kselect]]$rupt
  logdensity = t(apply(rupt,1,FUN=function(y) logdens(   x[ y[1]:y[2] ] ,P,phi)))
  tau        = Estep(logdensity,phi)
  pop        = apply(tau,1,which.max)
  cluster    = c()
  t.est      = bpwmissing.calls(t.est=rupt[,2],present.data,n.com,Kselect)
  if (Kselect>1){
    rupt[,1]   = c(1,t.est[1:(Kselect-1)]+1)
  } else {
    rupt[,1] = 1
  }
  rupt[,2]   = t.est
  cluster    = c()
  for (k in (1:Kselect)){
    cluster[rupt[k,1]:rupt[k,2]] = rep(pop[k],rupt[k,2]-rupt[k,1]+1)
  }
  mu = data.frame(begin = rupt[,1],
    end   = rupt[,2],
    mean  = phi[pop],
    levels= pop) 
  invisible(list(mu=mu,loglik=loglik))      

}

bpwmissing.calls <- function(t.est,present.data,n.com,Kselect){
  for (h in 1:Kselect){
    if (length(which(t.est[h]==0))!=0){
      t.est[h][-which(t.est[h]==0)] = present.data[t.est[h]]
    } 
    else {t.est[h] = present.data[t.est[h]]}
  }
  t.est[Kselect] = n.com
  invisible(t.est)
}



