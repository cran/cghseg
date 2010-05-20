setMethod(f = "multisegmean",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){
           
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            multiKselect = multiKmax
            M            = length(names(.Object@Y))
            Kseq         = c(M:multiKmax)

######   Individual segmentations for patients   
            
            Res = lapply(names(.Object@Y), FUN = function(m){
              n                           = length(which(!is.na(.Object@Y[[m]])))
              Kmax                        = uniKmax[[m]]
              out                         = unisegmean(.Object@Y[[m]],CGHo,Kmax)
              J.est                       = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
              invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
            }) 
            names(Res) = names(.Object@Y)

######   Segment Repartition segments across patients 
            
            J.est              = lapply(Res,FUN = function(x){x$J.est})
            nbdata             = lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) 
            nbdata             = sum(unlist(nbdata))    
            out.ibp            = segibp(data.frame(J.est = unlist(J.est)),unlist(uniKmax),multiKmax)
            multiloglik        = -(nbdata/2)*(log(2*pi*out.ibp[[1]]/nbdata)+1)
            seg.rep            = out.ibp[[2]]
            row.names(seg.rep) = names(.Object@Y)


######   Model Selection  
                        
            if (select.tmp=="none"){
              multiKselect = multiKmax
              dimll        = length(multiloglik)
            } else if (select.tmp=="mBIC"){    
              mBIC = c()
              i    = 0
              for (K in Kseq){
                i       = i+1
                mu      = multisegout(.Object,seg.rep,Res,K)
                mBIC[i] = getmBIC(K,multiloglik[i],mu,CGHo)        
              }
              multiKselect = Kseq[which.max(mBIC)]
              dimll        = multiKselect
            }
            
            
######   Outputs   
            
            mu           = multisegout(.Object,seg.rep,Res,multiKselect)
            select(CGHo) = select.tmp
            invisible(list(mu=mu,loglik=multiloglik[dimll],nbiter=0))
          } 
          )
