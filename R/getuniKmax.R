setMethod(f = "getuniKmax",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax=NULL){
            
            if (is.null(uniKmax)){              
              if (CGHo["select"] == "none"){
                cat("[getuniKmax] if no selection is performed while uniKmax is not specified \n")
                cat("[getuniKmax] uniKmax is initialized by default\n")
              }
              uniKmax = lapply(names(.Object@Y),FUN = function(ell){
                floor(sum(!is.na(.Object@Y[[ell]]))*CGHo["alpha"])
              })
              names(uniKmax) = names(.Object@Y)                                           
            } else {              
              if ( sum(names(.Object@Y) != names(uniKmax))>0){                
                cat("[getuniKmax] uniKmax should be a list with the same names as CGHd \n")
                cat("[getuniKmax] the names of uniKmax should be patients names \n")      
                stop()
              }              
              lapply(names(.Object@Y), FUN = function(m){
                if (CGHo["calling"]){
                  if ((CGHo["nblevels"]>uniKmax[[m]]) ){
                    cat("[getuniKmax] Error in uniKmax \n")
                    cat("[getuniKmax] Error for profile: ", m, "\n")
                    cat("[getuniKmax] The number of clusters must be lower than the number of segments in this profile","\n")
                    cat("[getuniKmax] Check CGHo[\"nblevels\"] and Kmax for this profile","\n")
                    stop()
                  }
                }
                n.com  = length(.Object@Y[[1]])
#                if (uniKmax[[m]]> floor(n.com/CGHo["lmin"])){
#                  cat("[getuniKmax] Error in uniKmax \n")
#                  cat("[getuniKmax] Error for profile: ", m, "\n")
#                  cat("[getuniKmax] The max number of segments for this profile must be lower than [length(x)/lmin]","\n")
#                  cat("[getuniKmax] Check CGHo[\"lmin\"] and uniKmax for this profile","\n")                
#                  stop()
#                }                
#                if (uniKmax[[m]]*CGHo["lmax"]<1){
#                  cat("[getuniKmax] Error in uniKmax \n")
#                  cat("[getuniKmax] Error for profile: ", m, "\n")
#                  cat("[getuniKmax] The max number of segments for this profile must be greater than [length(x)/lmax]","\n")
#                  cat("[getuniKmax] Check CGHo[\"lmax\"] and uniKmax for this patient","\n")
#                  stop()          
#                }
              })
            }
            return(uniKmax)  
          })



























