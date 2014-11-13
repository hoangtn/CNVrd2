groupBayesianCNVs <- function(xData, nGroups,
                              lambda0,
                              sd0,
                              alpha0,
                              distanceBetweenGroups ,
                              inits = NULL,

                              precisionOfGroupMeans = 3000,
                                                            
                              sdOftau = NULL,
                              n.adapt = 100, nUpdate = 1000, n.iter = 20000, 
				thin = 5, n.chains = 1,
				heidel.diag = FALSE,
                              leftLimit = NULL, rightLimit = NULL) {

    #f(x|...) ~ sum[pj*fj(x|mj, sj)]
    ####################################Check data##################################
    ############Set limits at both sides of segmentation scores################
            if (is.null(leftLimit))
                leftLimit <- min(xData)
            if (is.null(rightLimit))
                rightLimit <- max(xData)


###########Keep only values being inside the limits########################
            x <- xData
            xData <- x[x>=leftLimit & x<= rightLimit]

    ################################################################################
    xData <- sort(xData)
    N <- length(xData) #length of data vector
    nGroups <- nGroups #number of groups

    lambda0 <- lambda0 # initial means
    nLambda <- length(lambda0)
    
    sd0 <- sd0 #initial standard deviations
    tau0 <- 1/sd0^2 #precisions
    
    alpha0 <- alpha0 #prior information for the proportions of groups

    distanceBetweenGroups <- distanceBetweenGroups #prior distances between groups
    
    precisionOfGroupMeans <- precisionOfGroupMeans #prior precisions of groups

    ##tau ~ dgamma(alpha, beta)
    ##Make alphas' and betas' values from lamda0

    alphaOfTau <- betaOfTau <- c()

    if (is.null(sdOftau))
        sdOftau <- rep(1, nLambda) #Set standard deviation for tau


    for (ii in 1:nLambda){

        modeOfTau = tau0[ii] #Set tau0 as mode of tau

        betaOfTau[ii] = ( modeOfTau + sqrt( modeOfTau^2 + 4*sdOftau[ii]^2 ) ) / ( 2 * sdOftau[ii]^2 ) #Calculate shape parameter of tau
        alphaOfTau[ii] = 1 + modeOfTau*betaOfTau[ii] #Calculate rate parameter of tau
           }

       
    
    if (is.null(inits)){
        inits = list(lambda = c(lambda0, rep(NA, nGroups - nLambda)))
    }

        ##################################################################
        #########################################Make a data file
        dataJags = list(
            xData = xData, #data
            nGroups = nGroups, #number of groups
            N = N,
            alpha0 = alpha0,

            T = c(1, rep(NA, N-2), nGroups),

             
    #tau ~ dgamma(alpha, beta)
            alphaSD = alphaOfTau, #shape parameters for tau
            betaSD = betaOfTau, #rate parameters for tau

            lambda0 = lambda0, #
            nLambda = nLambda,
            distanceBetweenGroups = distanceBetweenGroups,
            precisionOfGroupMeans = precisionOfGroupMeans
            )
        


   jagsFile = "
        
        model {
            for( i in 1:N ) {
                xData[i] ~ dnorm(mu[i], tau[i])
                mu[i] <- lambda[T[i]]
                tau[i] <- theta[T[i]]
                T[i] ~ dcat(pi[])
                }
            pi[1:nGroups] ~ ddirch(alpha0[])
            for (hL in 1:nLambda){
                lambda[hL] ~ dnorm(lambda0[hL], precisionOfGroupMeans)
                }
            for (j in (nLambda+1):nGroups) {
                hh[j] ~ dnorm(distanceBetweenGroups, precisionOfGroupMeans)
                lambda[j] <- lambda[j -1] + hh[j]
                }
            for (gg in 1:(nLambda)){
                theta[gg] ~ dgamma(alphaSD[gg], betaSD[gg])
                }
            for (ii in (nLambda+1):nGroups){
                theta[ii] ~ dgamma(alphaSD[nLambda], betaSD[nLambda])
                }

            for ( kk in 1:nGroups ) {
                sigma[kk] <- sqrt(1/theta[kk])
                }
              }
         "

        mixture <- jags.model(textConnection(jagsFile), inits = inits,
               data = dataJags, n.chains = n.chains, n.adapt = n.adapt)

              


        update(mixture, nUpdate)

        xem <- coda.samples(mixture, c('lambda', 'pi', 'sigma'), n.iter = n.iter, thin = thin)
        
	hTest <- NULL
        if (heidel.diag == TRUE){
		hTest <- heidel.diag(xem)
		if (any(is.na(hTest[[1]][, 5])))

		warning("Not convergent\n", hTest)

		        if (any(is.na(hTest[[1]][, 5]))){
            message("===============================")
            message("==========Not convergent=======")
            message("===============================")
	}}



        allResults <- t(apply(xem[[1]], 2, function(x)
                    bb <- c(mean(x), quantile(x, c(0.25, 0.75)))))

       
            

        ##Obtain means, proportions and standard deviations
        m1 <- allResults[1:nGroups, 1]
        p1 <- allResults[(nGroups +1):(2*nGroups), 1]
        s1 <- allResults[(2*nGroups +1):(3*nGroups), 1]


    ###########################################################################
    ###################Determine group for data points


    zz <- matrix(0, nrow = length(xData), ncol = nGroups)
    for (ii in 1:length(xData)){
        for (jj in 1:nGroups){
            zz[ii, jj] <- p1[jj]*dnorm(xData[ii], mean = m1[jj], sd = s1[jj])/sum(p1*dnorm(xData[ii], mean = m1, sd = s1))
            }}

    zz0 <- apply(zz, 1, function(x)
                 which(x == max(x)))

    par(mfrow = c(2, 1))
    plot(xData, zz0, xlab = 'Segmentation score', ylab = 'Group', main = '')
    hist(xData, 100, xlab = 'Segmentation score', main = '')

    xydata <- data.frame(names(xData),  zz0, zz, xData)
    rownames(xydata) <- names(xData)
    colnames(xydata) <- c("Name", "Classification", paste("Group", 1:nGroups, sep = ""), "score")

############Force outliers into the largest groups or smallest group###
            x2 <- x[-pmatch(names(xData), names(x))]
            fx1 <- data.frame(names(x2), ifelse(x2 < leftLimit, 1, nGroups),
                              ifelse(x2< leftLimit, 1, 0),
                              data.frame(matrix(0, nrow = length(x2), ncol = nGroups -2)),
                              ifelse(x2 > rightLimit, 1, 0), x2)
            colnames(fx1) <- c("Name", "Classification", paste("Group", 1:nGroups, sep = ""), "score")
            xydata <- rbind(xydata, fx1)

            xydata <- xydata[order(rownames(xydata)),]





    ####################Return##############################################3

        returnValues <- list(mcmcChains = xem, m1 = m1, p1 = p1, s1 = s1, allGroups = xydata, hTest = hTest)

        return(returnValues)
    }

