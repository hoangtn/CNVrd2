setMethod("groupCNVs", "clusteringCNVs",
          function(Object, leftLimit = NULL, rightLimit = NULL, autoDetermineGroup = FALSE,
                   nGroupMin = 1, nGroupMax = 6, ...){
            k <- Object@k
            EV <- Object@EV
            x <- Object@x
            p <- Object@p
            m <- Object@m
            sigma <- Object@sigma
            small <- Object@small
            nMax <- Object@nMax
            eee <- Object@eee
            verbose = Object@verbose
            if (is.null(names(x)))
                names(x) <- paste('Sample', 1:length(x), sep = '')
############Set limits at both sides of segmentation scores################
            if (is.null(leftLimit))
                leftLimit <- min(x)
            if (is.null(rightLimit))
                rightLimit <- max(x)


###########Keep only values being inside the limits########################
            x1 <- x[x>=leftLimit & x<= rightLimit]
            Object@x <- x1
##########################################################################
	    if (autoDetermineGroup){
                k <- searchGroupCNVs(Object = Object, nGroupMin = nGroupMin, nGroupMax = nGroupMax)$nComponents
                Object@k <- k
		cat(k, " was chosen\n")

            }
      

            
############Run EM algorithm##############################################
            resultsEM <- emnormalCNV(Object)
            z <- resultsEM$z
            getMax <- function(x)
                ifelse(length(which(x == max(x))) > 1, 2, which(x == max(x)))
############Make data frame for clustering results#######################
            x1Label <- data.frame(names(x1), apply(z, 1, getMax), z, x1)
            colnames(x1Label) <- c("Name", "Classification", paste("Group", 1:k, sep = ""), "score")
############Force outliers into the largest groups or smallest group###
            x2 <- x[-pmatch(x1, x)]
            fx1 <- data.frame(names(x2), ifelse(x2 < leftLimit, 1, k),
                              ifelse(x2< leftLimit, 1, 0),
                              data.frame(matrix(0, nr = length(x2), ncol = k -2)),
                              ifelse(x2 > rightLimit, 1, 0), x2)
            colnames(fx1) <- c("Name", "Classification", paste("Group", 1:k, sep = ""), "score")
            xLabel <- rbind(x1Label, fx1)
#############################################Plots##################
            par(mfrow = c(2, 1))
            plot(xLabel[, k + 3], xLabel[, 2], xlab = 'Segmentation score', ylab = 'Group',
                 main = '', axes = FALSE)
            axis(2, 0:k)
            axis(1)
            abline(v = c(leftLimit, rightLimit), col = 'green')
            text(leftLimit, k, paste('leftLimit = ', round(leftLimit, 2), sep = ''), cex=0.7, srt= 90, col = 'red', pos = 2)
            text(rightLimit, k, paste('rightLimit = ', round(rightLimit, 2), sep = ''), cex=0.7, srt= 90, col = 'red', pos = 2)
            hist(xLabel[, k +3], main = '',
                 xlab = 'Segmentation score', 100)
  ########################################################################
            return(list(allGroups = xLabel, means = resultsEM$m, sigma = resultsEM$sigma, loglk = resultsEM$loglk, p = resultsEM$p))

})
