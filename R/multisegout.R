setMethod(f = "multisegout",signature = "CGHdata",
          definition = function(.Object,seg.rep,Res,Kselect){
            
            M  = length(names(.Object@Y))
            
            out = lapply(names(.Object@Y),FUN = function(m){
              i        = which(row.names(seg.rep) == m)
              k        = seg.rep[i,Kselect-M+1]
              rupt     = matrix(Inf,ncol = 2 , nrow= k)
              rupt[,2] = Res[[i]]$t.est[k,1:k]    
              if (k==1){
                rupt[1,1] = 1
              } else {
                rupt[,1] = c(1,rupt[1:(k-1),2]+1)
              }
              mu       = data.frame(begin = rupt[,1],
                end   = rupt[,2],
                mean  = apply(rupt,1,FUN=function(z) mean(.Object@Y[[m]][z[1]:z[2]], na.rm=T)))
              invisible(mu)
            })
            names(out) = names(.Object@Y)
            invisible(out)
          })
