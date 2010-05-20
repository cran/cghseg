setMethod(f = "ILS",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){

            tol          = 1e-6
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            command      = parse(text = "invisible(list(mu = mu, theta = B,loglik = loglik,nbiter = iter))")
    
            nbdata   = lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) 
            nbdata   = sum(unlist(nbdata))    
            M        = length(names(.Object@Y))
            n.com    = length(.Object@Y[[1]])
            eps      = Inf
            iter     = 0
            
            ## first iteration to initialize the epsilon algorithm
            mu       = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            nk       = unlist(lapply(mu,function(x){x$end-x$begin+1}))
            muk      = unlist(lapply(mu,function(x){x$mean}))
            B        = list(waveffect = rep(0,n.com), GCeffect = rep(0,n.com))
            param    = list(tm1 = rep(muk,nk),
                              t = rep(muk,nk),
                            tp1 = rep(muk,nk))
         
            ## second iteration to initialize the epsilon algorithm
            B                   = getbias(.Object,CGHo,mu,B)		
            removebias(.Object) = B$waveffect+B$GCeffect
            mu                  = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            revertbias(.Object) = B$waveffect+B$GCeffect            
            pred  = lapply(names(.Object@Y),FUN = function(m){
              nk  = mu[[m]]$end-mu[[m]]$begin+1
              muk = mu[[m]]$mean
              invisible(B$waveffect+B$GCeffect+ rep(muk,nk))
            })
            names(pred)   = names(.Object@Y)
            pred          = unlist(pred,use.names=TRUE)
            param$t       = pred          
            param.dot.tm2 = param$t
            
            while ( (eps > tol) & (iter < CGHo@itermax)){

              iter     = iter+1

              B                   = getbias(.Object,CGHo,mu,B)		
              removebias(.Object) = B$waveffect+B$GCeffect
              mu                  = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
              revertbias(.Object) = B$waveffect+B$GCeffect
              
              pred  = lapply(names(.Object@Y),FUN = function(m){
                nk  = mu[[m]]$end-mu[[m]]$begin+1
                muk = mu[[m]]$mean
                invisible(B$waveffect+B$GCeffect+ rep(muk,nk))
              })
              names(pred) = names(.Object@Y)
              pred        = unlist(pred,use.names=TRUE)
              param$tp1   = pred
                            
              param.dot.tm1 = param$t + invnorm( invnorm(param$tm1-param$t) + invnorm(param$tp1-param$t) )

              param$tm1     = param$t
              param$t       = param$tp1         
              eps           = sum( (param.dot.tm1-param.dot.tm2)^2 )
              param.dot.tm2 = param.dot.tm1
             
            } # end while
            
            loglik       = ILS.loglik(.Object,mu,B)
            select(CGHo) = select.tmp
            eval(command)
            
          })


setMethod(f = "ILS.loglik",signature = "CGHdata",
          definition = function(.Object,mu,bias){
            n   = lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) 
            n   = sum(unlist(n))    
            RSS = lapply(names(.Object@Y),FUN = function(m){
              nk      = mu[[m]]$end -  mu[[m]]$begin + 1
              rss     = sum( (.Object@Y[[m]] - rep(mu[[m]]$mean,nk) - bias$waveffect - bias$GCeffect)^2, na.rm = TRUE)
            })
            RSS    = 0.5* sum(unlist(RSS))
            loglik = -  (n/2)*(log(2*pi*RSS/n)+1)
            invisible(loglik)
          })

invnorm <- function(x){x/sum(x^2)}
