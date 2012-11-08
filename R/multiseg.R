setMethod(f = "multiseg",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){
            
            uniKmax    = getuniKmax(.Object,CGHo,uniKmax)
            multiKmax  = getmultiKmax(.Object,CGHo,uniKmax,multiKmax)  
            CGHr       = new("CGHresults",CGHd=.Object,CGHo=CGHo)
            from(CGHr) = "multiseg"
            Res        = list()
            
            if ((CGHo["GCnorm"]=="linear") & (is.null(.Object@GCcontent))){
              cat("[check multiseg] GCnorm can not be performed \n")
              cat("[check multiseg] check that GCcontent is a record in the data\n")
              stop()
            }
            
            if (CGHo@nbprocs>1){	
				if (Sys.info()["sysname"] == "Windows"){	
              		CGHo@cluster <- makeCluster(getOption("cl.cores", CGHo@nbprocs))
				}
            }
            
            if (CGHo["select"] != "none"){
              if ( (CGHo["calling"]==FALSE) & (CGHo["wavenorm"]=="none")  & (CGHo@GCnorm=="none")){
                cat("[multiseg] multisegmean running \n")
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
                Res = multisegmean(.Object,CGHo,uniKmax,multiKmax)
              } else {
                Kh        = golden.search(.Object,CGHo,uniKmax,multiKmax)
                multiKmax = Kh
                eval(fun2run(CGHo))
              }              
            } else {
              eval(fun2run(CGHo))
            }
            
            if (CGHo@nbprocs>1){
				if (Sys.info()["sysname"] == "Windows"){
              	stopCluster(CGHo@cluster)
				}
            }
            cat("\n")
            
            if ( (CGHo["wavenorm"]!="none") | (CGHo["GCnorm"]!="none")  ){
              theta(CGHr) =  Res$theta
            }
            
            if (!is.null(.Object["genomic.position"])){
              x      = .Object["genomic.position"]
              Res$mu = lapply(Res$mu,FUN = function(y){
                y$begin = x[y$begin]
                y$end   = x[y$end]
                return(y)
              })
            }
            
            mu(CGHr)     = Res$mu 
            loglik(CGHr) = list(multiloglik = Res$loglik)
            nbiter(CGHr) = Res$nbiter
            return(CGHr)
            
          } # end function
          )
