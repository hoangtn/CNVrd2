setMethod("emnormalCNV", "clusteringCNVs",
          function(Object, init.method = "kmeans", ...){
            k <- Object@k
            EV <- Object@EV
            x <- Object@x
            p <- Object@p
            m <- Object@m
            sigma <- Object@sigma
            small <- Object@small
            nMax <- Object@nMax
            eee <- Object@eee
            nmaxInit = Object@nmaxInit
            nChangeVariance = Object@nChangeVariance
            verbose = Object@verbose
            groupDistance = Object@groupDistance
            stopEV = 0

#####Verify the distance between groups##############################
            if (!is.null(groupDistance))
              if ((max(x) - min(x))/k < groupDistance)
                stop("\nPlease change the distance between groups: should be smaller than ", groupDistance, "\n")
#################Log likelihood function#############################
            loglk <- function(x, p, m, sigma){
              return(sum(vapply(x, function(x){
                log(sum(p*dnorm(x, mean = m, sd = sigma)))}, FUN.VALUE = 0)))}
########################Obtaining initial values####################
###################################################################
            initialValues <- function(x, nCenters, nmaxInit, init.method = init.method){
                p.init <- mu.init <- sigma.init <- c()
###########Use kmeans to obtain initial values: run nmaxInit times

                if (nCenters > 1){
                    mTemp <- pTemp <- sTemp <- matrix(0, ncol = nCenters,
                                                      nrow = nmaxInit)
                    mSum <- rep(NA, nCenters)
                                           
###########Use means of these blocks to input initial values for kmeans
                    for (jj in 1:nmaxInit){
                       tempAll <- kmeans(x, centers = nCenters)
                        mTemp[jj, ] <- tempAll$centers
                        
                        breakInit <- TRUE
                        if (!is.null(groupDistance)){

                          for (iiM in 2:nCenters){
                            mTempTest <- sort(mTemp[jj, ])
                                                       
                            if (abs(mTempTest[iiM] - mTempTest[iiM -1]) < groupDistance){
                              breakInit <- FALSE
                              break
                              }}
                           }
                        
                        if (breakInit == TRUE){

                          pTemp[jj, ] <- tempAll$size/length(x)
                          sTemp[jj, ] <- sqrt(tempAll$withinss/(tempAll$size))

                          mSum[jj] <- tempAll$tot.withinss

                     }
                    }
                    if (all(is.na(mSum)))
                      stop("\n\nPlease lessen groupDistance or increase nmaxInit (group numbers = ", nCenters, ") to obtain initial values\n\n")
                      
                    bInDex <- which(mSum == min(mSum, na.rm = TRUE))[1]
                    
                    
                    tempAll <- data.frame(
                        mTemp[bInDex, ],
                        pTemp[bInDex, ],
                        sTemp[bInDex, ])
                    tempAll <- tempAll[order(tempAll[, 1]),]
                    mu.init <- tempAll[, 1]
                    p.init <- tempAll[, 2]
                    sigma.init <- tempAll[, 3]


                    }
                else {
                    mu.init <- mean(x)
                    sigma.init <- sd(x)
                    p.init <- 1
                    }
                if (verbose){
                  message("===================================")
                  message("==========Initial values============")
                  message("p: ", p.init, "\nm: ", mu.init, "\nsigma: ", sigma.init)
                  message("===================================\n")
                }
                list(p.init = p.init,
                     mu.init = mu.init,
                     sigma.init = sigma.init)
                }
#########################################################################
#########################################################################
            inits <- initialValues(x = x, nCenters = k, nmaxInit = nmaxInit, init.method = init.method)
            warn <- options(warn = -1)

            if (is.null(p)) p = inits$p.init
            if (is.null(m)) m = inits$mu.init
            if (is.null(sigma))
              sigma <- ifelse(EV, rep(sd(x), k),
                              inits$sigma.init)
            eStop <- 2
            LLK <- c()
            count <- 0
            ndata <- length(x)
            z <- data.frame(matrix(0, nrow = ndata, ncol = k))

##################Main program############################################
            if (k == 1){
                llk1 <- sum(log(dnorm(x, mean = m, sd = sigma)))
                z <- matrix(1, nrow = length(x), ncol = 1)
                count <- 1
                LLK <- llk1
                }
            else {
##############################################################################
#############################################################################
#######Loops if n components > 1#############################################
                repeat {
                    llk0 <- loglk(x, p, m, sigma)
                    m0 <- m
                    if (verbose){
                        message("Loop: ", count, "\n  p:     ", p, "\n  m:     ", m, "\n  sigma: ", sigma)
                        }
###########################################################################
##########################################################################
                    ####E steps
                    z <- t(sapply(x, function(x){p*dnorm(x, mean = m, sd = sigma)}))
                    for (iZ in 1:dim(z)[1]){
                        if (sum(z[iZ, ])!=0)
                            z[iZ, ] <- z[iZ, ]/sum(z[iZ, ])
                        }

                    ####M steps
                    ##################################
                    colSumsZ <- colSums(z)
                    colSumsZ[which(colSumsZ == 0)] <- eee
                    m <- t(x)%*%z/colSumsZ
                    p <- colSumsZ/ndata
                    zk <- sapply(x, function(x) return(x - m))^2
                    ###Equal variances (EV)
                    if (EV == TRUE){
                        sigma <- sqrt(sum(z*t(zk))/ndata)
                        sigma <- rep(sigma, k)
                        }
                    ###Unequal variances (UEV)
                    ###If there is any sigma == 0 after nChangeVariance, UEV will be changed to EV
                    if (EV == FALSE){
                        sigma <- sqrt(diag(zk%*%z)/colSumsZ)
                          if (length(which(sigma <= eee)) > 0){
                              sigma[which(sigma <= eee)] <- 1.5*eee
                              }
                        testSigma <- any(sigma <= 2*eee)
                        if (testSigma == TRUE)
                            stopEV <- stopEV + 1
                        else
                            stopEV <- stopEV
                        if (stopEV >= nChangeVariance)
                            EV <- TRUE
                        }

                     llk1 <- loglk(x, p, m, sigma)
                     count <- count + 1
                    LLK[count] <- llk1
                    testLoglk <- abs(1 - llk1/llk0)
                    testRun <- (testLoglk > small) & (count <= nMax)

                    if (testRun == FALSE){
                        break}
                    }
                }
            if (EV == FALSE){
                bic <- 2*llk1 - (3*k - 1)*log(ndata)
                ev <- "Unequal variances"
                }
            else {
                bic <- 2*llk1 - (2*k)*log(ndata)
                ev <- "Equal variances"
                }
              message("=====================================")
              message(count - 1, " iterations")
              message(k, " components with ", ev)
              message("m: ", m)
              message("p: ", p)
              message("sigma: ", sigma)
              message("====================================")

            if (count >= nMax)
                warning(count -1, "iterations are maximum")

            

            return(list(loglk = llk1, p = p, m = m,
                        sigma = sigma, count =count, z = z,
                        bic = bic))

            })
