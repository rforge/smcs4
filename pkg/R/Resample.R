#smc_resample_systematic <- function(logWeights,...) {
#    sw <- sum(weights)
#    u1 <- runif(1,0,sw/N)
#    u <- c(0,seq(sw/N,sw*(N-1)/N,length=N-1)) + u1
#    times <- hist(u,breaks=c(0,cumsum(weights)),plot=FALSE)$counts
#    ids <- rep(1:N,times=times)
#    return(ids)
#}

smc_resample_multinomial <- function(logWeights,...) {
    N <- length(logWeights)
    weights <- exp(logWeights)
	  times <- rowSums(rmultinom(N,size=1,prob=weights))
    ids <- rep(1:N,times=times)
    return(ids)
}

smc_resample_residual <- function(logWeights,...) {
    N <- length(logWeights)
	  weights <- exp(logWeights)
	  Nw <- floor(N*weights/sum(weights))
	  p <- weights - Nw/N
	  times <- Nw + rowSums(rmultinom(N,size=1,prob=p))
    ids <- rep(1:N,times=times)
    return(ids)
}

setMethod("Resample",signature(object="ParticleMatrix"),
	function(object,type=c("systematic","residual","multinomial","stratified"),...) {
		type <- match.arg(type)
		logW <- logWeights(object,...)
		N <- object@N
		switch(type,
			systematic = {
				ids <- .Call("resample_systematic",logWeights = logW,PACKAGE="SMCS4")
			},
			residual = {
				ids <- smc_resample_residual(logW,...)
			},
			multinomial = {
				ids <- smc_resample_multinomial(logW,...)
			},
			stratified = {
				ids <- .Call("resample_stratified",logWeights = logW,PACKAGE="SMCS4")
			}
		)
		if(object@p_margin == 2) object@particles <- object@particles[,ids,drop=FALSE] else object@particles <- object@particles[ids,,drop=FALSE]
		logWeights(object) <- rep(log(1/N),N)
		object@unifWeights <- TRUE
		object
	}
)