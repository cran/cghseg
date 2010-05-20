unisegmean <- function(Y,CGHo,Kmax){
  
  n.com        = length(Y)
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  #lmax         = floor(CGHo["lmax"]*length(x))    
  #out         = segmean(x,Kmax,CGHo["lmin"],lmax,vh = TRUE)  
  #out          = segmeanGR(x,Kmax)
   out          = segmeanCO(x,Kmax)
  loglik       = -(length(x)/2)*(log(2*pi*out$J.est/(length(x)))+1)
  #t.est        = bpwmissing(out$t.est,present.data,n.com)
  
  if (CGHo["select"]=="none"){
    Kselect = Kmax
  } else {    
    Kseq    = c(1:Kmax)
    Kselect = Kselection(x,out$t.est,out$J.est,Kseq,CGHo)
  }
  
  t.est        = bpwmissing(out$t.est,present.data,n.com)
  th      = t.est[Kselect,1:Kselect]
  rupt    = matrix(ncol=2,c(c(1,th[1:Kselect-1]+1),th))    
  mu      = data.frame(begin = rupt[,1],
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
