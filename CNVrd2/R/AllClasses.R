setClass("CNVrd2",
         representation(windows = "numeric", chr = "character",
                        st = "numeric", en = "numeric",
                        dirBamFile = "character", dirCoordinate = "character",
                        genes = "numeric", geneNames = "character"),
         prototype(windows = NULL, chr = NA_character_,
                   st = NULL, en = NULL,
                   dirBamFile = NA_character_, 
                   dirCoordinate = NA_character_,
                   genes = NULL,
                   geneNames = NA_character_))

setClassUnion("numericOrNULL", c("numeric","NULL"))

setClass("clusteringCNVs", representation(x = "numeric", k = "numeric", p = "numericOrNULL",
                                     m = "numericOrNULL", sigma= "numericOrNULL", small = "numeric",
                                     nMax = "numeric", EV = "logical", eee = "numeric",
                                     nmaxInit = "numeric", nChangeVariance = "numeric",
                                           verbose = "logical", groupDistance = "numericOrNULL"),
         prototype(x = NULL, k = 3, p = NULL, m = NULL, sigma = NULL, small = 1e-4,
                   nMax = 50, EV = FALSE, eee = .Machine$double.eps,
                   nmaxInit = 100, nChangeVariance = 3, verbose = FALSE, groupDistance = 0.25))
