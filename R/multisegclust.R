setMethod(f = "multisegclust",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){			  
            
            P            = CGHo["nblevels"]
            tol          = 1e-2
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            command      = parse(text = "mu = multisegclust.output(.Object,mu,out.EM$phi,out.EM$tau) \n invisible(list(mu = mu, loglik = loglik,nbiter = iter))")
            nbdata       = Reduce("sum",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}))
            iter         = 0
            eps          = Inf		
            
            if (CGHo@nbprocs>1){
				if (Sys.info()["sysname"] == "Windows"){
					## Initial data sends, will be reused but not resend
					## Data are emulated to belong to .GlobalEnv
					## since worker function will also belong to .GlobalEnv
					assign("Y.ref", .Object@Y, envir = .GlobalEnv)
					clusterExport(CGHo@cluster, "Y.ref")
					assign("uniKmax.ref", uniKmax, envir = .GlobalEnv)
					clusterExport(CGHo@cluster, "uniKmax.ref")
					assign("CGHo.ref", CGHo, envir = .GlobalEnv)
					clusterExport(CGHo@cluster, "CGHo.ref")
				}
            }
            
	    mu        = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            out.DP2EM = DP2EM(.Object,mu)
            phi       = compactEMinit(out.DP2EM$xk,out.DP2EM$x2k,out.DP2EM$nk,P,CGHo@nbprocs,vh=TRUE)
            out.EM    = compactEMalgo(out.DP2EM$xk,out.DP2EM$x2k,phi,out.DP2EM$nk,P,vh=TRUE)            
            n.com     = length(.Object@Y[[1]])
            mu.test   = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau)
	    mu.tmp    = mu.test	
            
            while ( (eps > tol) & (iter < CGHo@itermax) ){   
              iter       = iter+1 
              mu         = multisegmixt(.Object,CGHo,uniKmax,multiKmax,out.EM$phi)$mu
              out.DP2EM  = DP2EM(.Object,mu)
	      out.EM     = compactEMalgo(out.DP2EM$xk,out.DP2EM$x2k,phi,out.DP2EM$nk,P,vh=TRUE)			  
              mu.test    = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau) 
              eps        = max(sapply(names(.Object@Y),FUN=function(m,x,y){xk = rep(x[[m]]$mean,x[[m]]$end-x[[m]]$begin+1); yk =rep(y[[m]]$mean,y[[m]]$end-y[[m]]$begin+1) ; return(max(abs((xk-yk)/xk)))},mu.tmp,mu.test))
	      mu.tmp     = mu.test
            } # end while
			
            
######   output   #####################################################################
            out.DP2EM    = DP2EM(.Object,mu)
            loglik       = quicklvinc(out.DP2EM$xk,out.DP2EM$x2k,out.EM$phi,out.DP2EM$nk,P,vh=TRUE)$lvinc
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
