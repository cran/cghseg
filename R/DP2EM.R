setMethod(f = "DP2EM",signature = "CGHdata",
          definition = function(.Object,mu){
            
            # warning : missing values are replaced by the mean for the position in each group
            n.com = length(.Object@Y[[1]])
            M     = length(names(.Object@Y))  
            nbseg = lapply(mu,FUN = function(x){dim(x)[1]})
            phase = rep(n.com,M)
            phase = rep(c(0,cumsum(phase))[1:M],unlist(nbseg))
            
            correct = mapply(.Object@Y,FUN=function(x){x})  
            correct = apply(correct,1,FUN = function(x){mean(x,na.rm=TRUE)})
            
            signal= lapply(.Object@Y, FUN = function(x){
              x[which(is.na(x))] = correct[which(is.na(x))]
              invisible(x)
            })
            signal        = unlist(signal)
            
            end = lapply(mu,FUN = function(x){x$end})
            end = unlist(end)
            begin = lapply(mu,FUN = function(x){x$begin})
            begin = unlist(begin)
            
            rupt = cbind(begin+phase,end+phase)
            invisible(list(rupt=rupt,signal=signal))
          }
          )


