# Apdapted from celda package 
.sampleNamedLl <- function(llProbs) {
		# names(llProbs) should only be integers
		if (is.null(names(llProbs))) {
				stop("Log-probability should be named vector.")
		}
    probsSub <- exp(llProbs - max(llProbs))
    probsNorm <- probsSub / sum(probsSub)
    probsSelect <- sample(
        #x = names(probsNorm),
				x = names(llProbs),
        size = 1L,
        replace = TRUE,
        prob = probsNorm
        )  # return the cluster name (in character type)
    probsSelect = as.integer(probsSelect) # change character type into integer type 
    return(probsSelect)  
}





Initialization = function(counts=sim$counts, SampleLabel=sim$sampleLabel, K =20, seed=123 ) {

    set.seed(seed)

    s.int = as.integer(SampleLabel)
    Sample = length(unique(SampleLabel))
    nC.byS = countwFixedL(labelVec=s.int, Length=Sample)

    nC = ncol(counts)  # Number of cells in the dataset 
    nG = nrow(counts)  # Number of genes in the dataset 

    #Initialize topic label for each cell 
    k.byC = sample(x=K, size=nC, replace=TRUE)
    #
    nC.KbyS = sapply(1:Sample, FUN=function(s) {
			     nC.sampleS = countwFixedL(labelVec=k.byC[s.int==s], Length=K) 
			     names(nC.sampleS) = as.character(1:K) 
			     #nC.sampleS = nC.sampleS[order(nC.sampleS, decreasing=TRUE)]
			     return(nC.sampleS)
				 }, simplify=FALSE) 
    #
    N.GbyK = countNFixedCol(M=counts, labelVec=k.byC, nG=nG, nCol=nC, RnCol=K) 
    # 
    #nC.byK = countwFixedL(labelVec=k.byC, Length=K) 


    iter.updates = list("k.byC"=k.byC, 'nC.KbyS'=nC.KbyS,'N.GbyK'=N.GbyK, 's.int'=s.int) 


    return(iter.updates)
}



iterationDP = function(counts = sim$counts, iter.updates, alpha=2, eta=1,K=20) {

    nC = ncol(counts) 
    nG = nrow(counts) 
    #Iterate through each cell 
    for (Cindex in 1:nC) { 

    N.byG.Cindex = counts[, Cindex]
    k.Cindex = iter.updates$k.byC[Cindex]
    s.Cindex = iter.updates$s.int[Cindex] 

    N.GbyK = iter.updates$N.GbyK
    N.GbyK[, k.Cindex] = N.GbyK[, k.Cindex] - N.byG.Cindex

    #Calculate log-probability (posterior-Dir term) of the current cell in each topic/cluster 
    Nplus.GbyK = N.GbyK + N.byG.Cindex
    logp.multi = colSums(lgamma(Nplus.GbyK + eta)) - lgamma(colSums(Nplus.GbyK) + nG*eta) + 
	   lgamma(colSums(N.GbyK) + nG*eta) - colSums(lgamma(N.GbyK + eta))  
    names(logp.multi) = as.character(1:K) 


    nC.byK.s = iter.updates$nC.KbyS[[s.Cindex]] 
    logp.multi = logp.multi[ order(nC.byK.s, decreasing=TRUE) ]
    nC.byT.s = nC.byK.s[ order(nC.byK.s, decreasing=TRUE) ]

    #Calculate log-probabiliry of the topic/cluster proportion (stick-breaking process) 
    b_SBP = reverseCumSum(nC.byT.s)
    logp.pi = SBP2Prop(a=nC.byT.s, b=b_SBP,alpha=alpha) 


    #Sample new cluster label for current cell
    new.k.Cindex = .sampleNamedLl( logp.pi + logp.multi)
    #Update
    N.GbyK[, new.k.Cindex] = N.GbyK[, new.k.Cindex] + N.byG.Cindex
    #Update
    nC.byK.s[ k.Cindex] = nC.byK.s[ k.Cindex] - 1
    nC.byK.s[ new.k.Cindex] = nC.byK.s[ new.k.Cindex] + 1

    iter.updates$k.byC[Cindex] = new.k.Cindex
    iter.updates$nC.KbyS[[s.Cindex]] = nC.byK.s
    iter.updates$N.GbyK = N.GbyK

    }

    return(iter.updates) 

}


#' @export
celdaC_dp = function(counts, SampleLabel, alpha =2, eta=1, K=20, seed=12345, iter = 30, stopIteri = 3){
    iter.updates = Initialization(counts = counts, SampleLabel = SampleLabel, K=K, seed=seed) 

		temp_z = iter.updates$k.byC
    temp_K = length(unique(temp_z))
		numIterWithoutImprovement = 0

    for (i in 1:iter) { 
        cat(format(Sys.time(),"%a %b %d %X %Y"), " running iteration: ",  i, "\n") 
        iter.updates = iterationDP(counts = counts, iter.updates=iter.updates, K=K, alpha=alpha, eta=eta) 

				#if (all.equal(temp_z, iter.updates$k.byC) == TRUE) {
        ## labels could shuffle after being updated
        match_tb = table(temp_z, iter.updates$k.byC)
        if (length(which(match_tb>0)) == temp_K) {
					numIterWithoutImprovement = numIterWithoutImprovement + 1
				} else {
					temp_z = iter.updates$k.byC
          temp_K = length(unique(temp_z))
					numIterWithoutImprovement = 0
				}

				if (numIterWithoutImprovement == stopIteri) {
					break
				}
    }


    return(list("resList" = iter.updates, "params" = list("alpha"=alpha, "eta"=eta, "K"=K, "iter"=iter))) 
}

