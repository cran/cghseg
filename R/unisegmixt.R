unisegmixt <- function(Y,CGHo,Kmax,phi){
  
  P            = CGHo["nblevels"]
  n.com        = length(Y)
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  #lmax         = floor(CGHo["lmax"]*length(x))
  #out          = segmixt(x,P,Kmax,phi,CGHo["lmin"],lmax)
  out          = segmixtGR(x,P,Kmax,phi)
  loglik       = -out$J.est
  t.est        = bpwmissing(out$t.est,present.data,n.com)
  
  if (CGHo["select"]=="none"){
    Kselect = Kmax
  } else if (CGHo["select"]!="none"){
    Kseq    = c(1:Kmax)
    Kselect = Kselection(Y,out$J.est,Kseq,CGHo)
  }
    
  th        = t.est[Kselect,1:Kselect]
  rupt      = matrix(ncol=2,c(c(1,th[1:Kselect-1]+1),th))    
  mu        = data.frame(begin = rupt[,1],
    end   = rupt[,2],
    mean  = apply(rupt,1,FUN=function(z) mean(Y[z[1]:z[2]], na.rm=T)))  

  invisible(list(mu=mu,loglik=loglik,t.est=t.est))
  
}

bpwmissing <- function(t.est,present.data,n.com){
  for (h in 1:ncol(t.est)){
    if (length(which(t.est[,h]==0))!=0){
      t.est[,h][-which(t.est[,h]==0)] = present.data[t.est[,h]]
    } 
    else {t.est[,h] = present.data[t.est[,h]]}
  }
  diag(t.est) = n.com
  invisible(t.est)
}
