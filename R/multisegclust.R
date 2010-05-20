setMethod(f = "multisegclust",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){

            P            = CGHo["nblevels"]
            tol          = 1e-6
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            command      = parse(text = "mu = multisegclust.output(.Object,mu,out.EM$phi,out.EM$tau) \n invisible(list(mu = mu, loglik = loglik,nbiter = iter))")
            nbdata       = lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) 
            nbdata       = sum(unlist(nbdata))    
            iter         = 0
            eps          = Inf
            
            mu        = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            out.DP2EM = DP2EM(.Object,mu)
            phi       = EMinit(out.DP2EM$signal,out.DP2EM$rupt,P,vh=TRUE)
            out.EM    = EMalgo(out.DP2EM$signal, phi, out.DP2EM$rupt, P, vh = TRUE)
            n.com     = length(.Object@Y[[1]])
            mu.test   = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau)
            nk        = unlist(lapply(mu.test,function(x){x$end-x$begin+1}))
            muk       = unlist(lapply(mu.test,function(x){x$mean}))
            param     = list(tm1 = rep(muk,nk),
                               t = rep(muk,nk),
                             tp1 = rep(muk,nk))


            ## second iteration to initialize the epsilon algorithm
            ## initialize param$t
            mu                  = multisegmixt(.Object,CGHo,uniKmax,multiKmax,out.EM$phi)$mu
            out.DP2EM           = DP2EM(.Object,mu)
            out.EM              = EMalgo(out.DP2EM$signal, out.EM$phi, out.DP2EM$rupt, P, vh = TRUE)
            mu.test             = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau) 
            pred                = lapply(names(.Object@Y),FUN = function(m){
              nk  = mu.test[[m]]$end-mu.test[[m]]$begin+1
              muk = mu.test[[m]]$mean
              invisible(rep(muk,nk))
            })
            names(pred)   = names(.Object@Y)
            pred          = unlist(pred,use.names=TRUE)
            param$t       = pred
            param.dot.tm2 = param$t


            
            while ( (eps > tol) & (iter < CGHo@itermax) ){    

              iter      = iter+1              
              mu        = multisegmixt(.Object,CGHo,uniKmax,multiKmax,out.EM$phi)$mu
              out.DP2EM = DP2EM(.Object,mu)
              out.EM    = EMalgo(out.DP2EM$signal,out.EM$phi, out.DP2EM$rupt, P, vh = TRUE)
              mu.test   = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau) 
              pred      = lapply(names(.Object@Y),FUN = function(m){
                nk  = mu.test[[m]]$end-mu.test[[m]]$begin+1
                muk = mu.test[[m]]$mean
                invisible(rep(muk,nk))
              })
              names(pred) = names(.Object@Y)
              pred        = unlist(pred,use.names=TRUE)
              param$tp1   = pred           

              param.dot.tm1 = param$t + invnorm( invnorm(param$tm1-param$t) + invnorm(param$tp1-param$t) )            
              param$tm1     = param$t
              param$t       = param$tp1         
              eps           = sum( (param.dot.tm1-param.dot.tm2)^2 )
              if (is.nan(eps)) {eps=0} 
              param.dot.tm2 = param.dot.tm1
              
            } # end while
            
######   output   #####################################################################
            
            loglik       = lvmixt.MSC(.Object,mu,out.EM$phi)
            select(CGHo) =  select.tmp
            eval(command)
            
          } #end function
          )

######   auxiliary functions for ILSclust      ########################################

setMethod(f = "multisegclust.output",signature = "CGHdata",
          definition = function(.Object,mu,phi,tau){
            
            M               = length(names(.Object@Y))
            levels          = apply(tau,1,which.max)
            mutmp           = cbind(mean = phi[levels],levels=levels)
            nk              = lapply(mu,FUN=function(x){length(x$mean)})
            end             = cumsum(nk)
            begin           = c(1,end[1:(length(end)-1)]+1)
            rupt            = cbind(begin,end)
            rownames(rupt)  = names(.Object@Y)
            mutmp           = apply(rupt,1,FUN = function(y){
              tmp           = data.frame(matrix(mutmp[y[1]:y[2],],ncol=2))
              colnames(tmp) = c("mean","levels")
              invisible(tmp)
            })
            mutmp = lapply(1:M,FUN = function(m){cbind(mu[[m]][,-3],mutmp[[m]])})
            names(mutmp) = names(.Object@Y)  
            invisible(mutmp)
          }
          )

setMethod(f = "lvmixt.MSC",signature = "CGHdata",
          definition = function(.Object,mu,phi){
            P     = length(phi)/3
            chrom = names(.Object@Y)
            lv = sum(unlist( lapply(names(.Object@Y),FUN = function(m){
              Ym    = .Object@Y[[m]]
              ruptm = mu[[m]][-3]
              lvinc.MSC(Ym,phi,ruptm,P)
            })))
            invisible(lv)
          }
          )


lvinc.MSC  <- function (x,phi,rupt,P){
  logdensity  = t(apply(rupt, 1,FUN = function(y){
    xk  = x[y[1]:y[2]]
    xk  = xk[!is.na(xk)]
    invisible(logdens(xk,P, phi))
  }))
  K       = nrow(logdensity)
  P       = ncol(logdensity)
  tau     = sapply(1:P,FUN = function(p){logdensity[,p]+log(phi[p+2*P])})
  tau     = matrix(tau,ncol=P)
  tau_max = apply(tau,1,max)
  tau     = exp(tau-matrix(rep(tau_max,P),ncol=P))
  lvinc   = sum(log( apply(tau,1,sum)) + tau_max)
  return(lvinc)
}


invnorm <- function(x){x/sum(x^2)}
