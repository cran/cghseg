setMethod(f = "getmultiKmax",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){
            
            M       = length(names(.Object@Y))
            if (is.null(uniKmax)){
              uniKmax = getuniKmax(.Object,CGHo,uniKmax)
            }
            
            if (is.null(multiKmax)){
              if (CGHo["select"] == "none"){
                message("[check getmultiKmax] if no selection is performed while multiKmax is not specified \n")
                message("[check getmultiKmax] multiKmax is initialized by default\n")
                message("**********************************************************\n")
              }
              multiKmax    = floor(sum(unlist(uniKmax))) * CGHo["beta"]
            } else {
              if (multiKmax < M){
                stop("[getmultiKmax]  multiKmax should be greater that the number of patients (",M,") \n")
              } else {
                if (multiKmax > sum(unlist(uniKmax))){
                stop("[getmultiKmax] multiKmax should be lower that the sum of the number allowed segments in every profile (",sum(unlist(uniKmax)),") \n")
                }
              }
            }
            return(multiKmax)
          }
          )
