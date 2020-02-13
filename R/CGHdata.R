setClassUnion("numOrNULL", c("numeric","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))


setClass(Class = "CGHdata",representation(Y = "list",GCcontent="numOrNULL",genomic.position="numOrNULL",probeID="factorOrNULL"))

setMethod(
          f          = "initialize",
          signature  = "CGHdata",
          definition =          
          function(.Object,Y) {
            
            Y = data.frame(Y)            
            patient.status = ""

            if ("GCcontent" %in% colnames(Y)){
              j = which(colnames(Y)=="GCcontent")
              .Object@GCcontent = Y[,j]
              Y = Y[,-j]
            } else {
              .Object@GCcontent=NULL
            }
            if ("genomic.position" %in% colnames(Y)){
              if (is.unsorted(Y$genomic.position)){
                stop("[raw data: error] sort genomic.positions \n")
              } else {              
                j = which(colnames(Y)=="genomic.position")
                .Object@genomic.position = Y[,j]
                Y = Y[,-j]
              }
            } else {
              .Object@genomic.position = NULL
            }
            
            if ("probeID" %in% colnames(Y)){
              if (sum(duplicated(Y$probeID)>0)){
                stop("[raw data: error] remove duplicated names in probeID \n")
                } else {                
                j = which(colnames(Y)=="probeID")              
                .Object@probeID = Y[,j]
                Y = Y[,-j]
              }
            } else {
              .Object@probeID = NULL
            }
            if (dim(Y)[2]==1){
              patient.status = "U"
            } else {
              patient.status = "M"
            }
                        
            ########################################
            #                                      #
            # check data format on the dataframe   #
            #                                      #
            ########################################
            
            message("[raw data: check] minimum number of positions per signal \n")
            n.com    = dim(Y)[1]
            min.ncom = 10
            
            if ( n.com<min.ncom){
             stop("[raw data: check] Problem with dataset \n It must contain more than ", min.ncom , " points \n")
            }            
            
            if (patient.status=="M"){
              message("[raw data: check] number of records per position \n")
              j        = which(apply((apply(Y,1,is.na)),2,sum)==dim(Y)[2])   
              if (length(j) > 0){
                message("[raw data: check] Problem with non observed positions \n")
                message("Position(s) number", j ," is/are missing in every patient(s) \n")
              }
            } else {
              message("[raw data: check] number of records per position \n")
              if (sum(is.na(Y))>0){
                message("[raw data: check] Problem with non observed positions \n")
                message("Position(s) number", which(is.na(Y))," is/are missing \n")
              }
            }
            message("[raw data: format] changing dataframe to compact format (list) \n")
            
            Res = lapply(names(Y), FUN = function(m){
              Y[,names(Y)==m]
            })
            names(Res) = names(Y)
            .Object@Y  = Res            
            return(.Object)
          }
          )
          
          
            


###############################################################
#                                                             #  
# GENERIC DEFINITION FOR CLASS CGHDATA                        #
#                                                             #
###############################################################


setGeneric( name = "getuniKmax"            ,def = function(.Object,CGHo,uniKmax=NULL){standardGeneric("getuniKmax")})
setGeneric( name = "getmultiKmax"          ,def = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){standardGeneric("getmultiKmax")})
setGeneric( name = "multisegout"           ,def = function(.Object,seg.rep,Res,Kselect){standardGeneric("multisegout")})
setGeneric( name = "multisegmean"          ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("multisegmean")})
setGeneric( name = "multisegmixt"          ,def = function(.Object,CGHo,uniKmax,multiKmax,phi){standardGeneric("multisegmixt")})
setGeneric( name = "ILS"                   ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("ILS")})
setGeneric( name = "ILSclust.output"       ,def = function(.Object,mu,phi,tau){standardGeneric("ILSclust.output")})
setGeneric( name = "ILSclust"              ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("ILSclust")})
setGeneric( name = "multisegclust.output"  ,def = function(.Object,mu,phi,tau){standardGeneric("multisegclust.output")})
setGeneric( name = "multisegclust"         ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("multisegclust")})
setGeneric( name = "multiseg"              ,def = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){standardGeneric("multiseg")})
setGeneric( name = "golden.search"         ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("golden.search")})
setGeneric( name = "uniseg"                ,def = function(.Object,CGHo,uniKmax=NULL){standardGeneric("uniseg")})
setGeneric( name = "correctdata"           ,def = function(.Object,bias){standardGeneric("correctdata")})
setGeneric( name = "revertdata"            ,def = function(.Object,bias){standardGeneric("revertdata")})
setGeneric( name = "getbias"               ,def = function(.Object,CGHo,mu,bias,phi=c(),tau=c()){standardGeneric("getbias")})
setGeneric( name = "DP2EM"                 ,def = function(.Object,mu,theta=NULL){standardGeneric("DP2EM")})




setGeneric("removebias<-",function(object,value){standardGeneric("removebias<-")})
setReplaceMethod(
                 f="removebias",
                 signature="CGHdata",
                 definition=function(object,value){
                   object@Y = lapply(object@Y,FUN=function(y){y-value})
                   return (object)
                 }
                 )
setGeneric("revertbias<-",function(object,value){standardGeneric("revertbias<-")})
setReplaceMethod(
                 f="revertbias",
                 signature="CGHdata",
                 definition=function(object,value){
                   object@Y = lapply(object@Y,FUN=function(y){y+value})
                   return (object)
                 }
                 )
setMethod(
          f = "[",
          signature = "CGHdata",
          definition = function(x,i,j,drop){   
            if (i=="Y")                  {return(x@Y)}                else {}
            if (i=="GCcontent")          {return(x@GCcontent)}        else {}
            if (i=="genomic.position")   {return(x@genomic.position)} else {}
            if (i=="probeID")            {return(x@probeID)}          else {}            
          }          
          )
setMethod(
          f = "summary",
          signature = "CGHdata",
          definition = function(object){
            message("****** Summary of CGHd object ******\n")
            message("[CGHd summary] Patients ID \n");
            message(names(object@Y),fill=TRUE,"\n")
            message("\n")
            message("[CGHd summary] number of points\n")
            message(length(object@Y[[1]]))
            message("[CGHd summary] probeID records\n")
            message(!is.null(object@probeID))
            message("[CGHd summary] genomic position \n")
            message(!is.null(object@genomic.position))
            message("[CGHd summary] GC content records\n")
            message(!is.null(object@GCcontent))
          }          
          )
setMethod(
          f = "show",
          signature = "CGHdata",
          definition = function(object){
            message("****** CGHdata show ******\n")
            message("[CGHd show] Data are in the list format [[patient]]\n")
            message("[CGHd show] Data sample: \n")
            message("Y[[",names(object@Y)[1],"]]\n",sep="")
            message(object@Y[[1]][1:5])
            if (length(names(object@Y))>1){
              message("Y[[",names(object@Y)[2],"]] \n",sep="")
              message(object@Y[[2]][1:5])
            }
            message("[CGHd show] probeID sample: \n")
            message(object@probeID[1:5])
            message("[CGHd show] genomic positions sample: \n")
            message(object@genomic.position[1:5])
            message("[CGHd show] GC content sample: \n")
            message(object@GCcontent[1:5])
          }          
          )

setMethod(
          f = "print",
          signature = "CGHdata",
          definition = function(x){
            message("****** CGHdata print ******\n")
            message("[CGHd print] Data are in the list format [[patient]]\n")
            message("[CGHd print] Data sample: \n")
            message("Y[[",names(x@Y)[1],"]]\n",sep="")            
            message(x@Y[[1]][1:5])
            if (length(names(x@Y))>1){
              message("Y[[",names(x@Y)[2],"]] \n",sep="")            
              message(x@Y[[2]][1:5])              
            }
            message("[CGHd print] probeID sample: \n")
            message(x@GCcontent[1:5])
            message("[CGHd print] genomic positions sample: \n")
            message(x@genomic.position[1:5])
            message("[CGHd print] GC content sample: \n")
            message(x@GCcontent[1:5])
          }          
          )

