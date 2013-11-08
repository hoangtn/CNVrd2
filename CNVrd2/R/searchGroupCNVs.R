##This function find suitable n components##################
##BIC is used##############################################

setMethod("searchGroupCNVs", "clusteringCNVs",
          function(Object, nGroupMin = 1, nGroupMax = 6, 
                   leftLimit = NULL, rightLimit = NULL,...){
              nGroupMax <- nGroupMax
              ObjectSearch <- Object
              verbose = Object@verbose
              if (!is.null(rightLimit))
                  ObjectSearch@x <- ObjectSearch@x[ObjectSearch@x <= rightLimit]
              if (!is.null(leftLimit))
                  ObjectSearch@x <- ObjectSearch@x[ObjectSearch@x >= leftLimit]
              testBIC <-  c()
              for (i in nGroupMin:nGroupMax){
                  ObjectSearch@k <- i
                  xLabel <- emnormalCNV(Object = ObjectSearch)
                  testBIC[i] <- xLabel$bic
                  }
            if (nGroupMin > 1){
                testBIC[1:(nGroupMin -1)] <- min(testBIC[nGroupMin:nGroupMax])
                }
              components <- paste("Using BIC: the best is ",
                                  which(testBIC == max(testBIC)), "components")
              nComponents <- which(testBIC == max(testBIC))
                        
    return(list(Groups = components, nComponents = nComponents))
  })


