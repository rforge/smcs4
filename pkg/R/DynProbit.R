#require(SMCS4)
#require(mvtnorm)
#require(msm)

# Dynamic Probit Regression
setClass("DynProbit",
	contains="ParticleMatrix",
	representation(
	    Q = "matrix",
	    Ptt = "matrix",
        yp = "vector",
        S = "numeric",      
        nx = "integer"
	)
)

DynProbitInit <- function(nx,mu0,Sigma0,Q,N,resampleC = .5) {
    new("DynProbit",
        particles = t(rmvnorm(N,mean=mu0,sigma=Sigma0)),
        logWeights = rep(log(1/N),N),
        unifWeights=TRUE,
        p_margin = as.integer(2),
        yp = rep(0,N),
        S = 0,
        Q = Q,
        Ptt = Sigma0,
        nx = as.integer(nx),
        resampleC = resampleC,
        N= as.integer(N)
    )
}

DynProbit <- function(formula,mu0,Sigma0,Q,N = 2000,data,resampleC = .5,output=c("mean","filter","smoothing")) {
    output <- match.arg(output)
    if(!missing(data)) {
        mf <- model.frame(formula,data=data)
    } else {
        mf <- model.frame(formula)
    }
    y <- model.response(mf)
    x <- as.matrix(model.matrix(terms(mf),mf))
    if(is.factor(y)) {
        y <- as.numeric(y != levels(y)[1])
    } else {
        if(!all(y %in% c(0,1))) stop("response should be a factor or binary variable")
    }
    nx <- ncol(x)
    
    if(length(mu0)!=nx) stop("mu0 should have length",nx)
    if(!is.matrix(Sigma0) | ncol(Sigma0) !=nx) stop("Sigma0 should be a square matrix with",nx,"rows")
    if(!is.matrix(Q) | ncol(Q) !=nx) stop("Q should be a square matrix with",nx,"rows")
    
    pf <- DynProbitInit(nx=nx,mu0=mu0,Sigma0=Sigma0,Q=Q,N=N,resampleC = resampleC)
    
    T <- length(y)
    
    if(output == "smoothing") {
        particles <- array(NA,dim=c(T,nx,N))
    } else {
        pmean <- matrix(NA,nrow=T,ncol=nx)
    }
    for(t in 1:T) {
        pf <- ParticleMove(pf,x=x[t,],y=y[t])
        pf <- UpdateWeights(pf,x=x[t,],y=y[t])
        if(output == "smoothing") {
            particles[t,,] <- particles(pf)
        }
        if(output == "mean") pmean[t,] <- mean(pf)
        if(pf@resampleC == 1 || ESS(pf) < resampleC*pf@N) {
            if(output == "smoothing") {
                idx <- .Call("resample_systematic",logWeights=logWeights(pf),PACKAGE="SMCS4")
                particles <- particles[,,idx]
                particles(pf) <- particles(pf)[,idx]
            } else {
                pf <- Resample(pf)
            }
        }
    }
    if(output == "mean") return(pmean)
    if(output == "filter") return(pf)
    if(output == "smoothing") return(list(particles=particles,weights=getNormWeights(pf)))
}

setMethod("ParticleMove",signature(object="DynProbit"),
	function(object,x,y,...) {
        object@Ptt <- object@Ptt + object@Q
        # assumes y is binary
        object@S <- as.numeric(t(x)%*%object@Ptt%*%x + 1)
        object@yp <- as.vector(t(x)%*%object@particles)
        # sample ypart
        if(y==0) ypart <- rtnorm(getN(object),mean=object@yp,sd=sqrt(object@S),lower=-Inf,upper=0) else ypart <- rtnorm(getN(object),mean=object@yp,sd=sqrt(object@S),lower=0,upper=Inf)
        object@particles <- object@particles + object@Ptt%*%x%*%(1/object@S)%*%(ypart-object@yp)
        object@Ptt <- (diag(object@nx)-(object@Ptt%*%x%*%(1/object@S))%*%t(x))%*%object@Ptt
		object
	}
)

setMethod("SmcIterate",signature(object="DynProbit"),
	function(object,x,y,...) {
		object <- ParticleMove(object,x,y,...)
		object <- UpdateWeights(object,x,y,...)
		ess <- ESS(object,...)
		if(ess < (object@resampleC*object@N)) object <- Resample(object,...)
		object
	}
)

setMethod("UpdateWeights",signature(object="DynProbit"),
    function(object,x,y,...) {
        weights <- logWeights(object) + pnorm(-object@yp/sqrt(object@S),lower.tail=!as.logical(y[1]),log.p=TRUE)
        #weights <- pnorm(-object@yp/sqrt(object@S),lower.tail=as.logical(y[1]),log.p=TRUE) 
        logWeights(object) <- weights
        object@unifWeights <- FALSE
        object
    }
)

setMethod("logLik",signature(object="DynProbit"),
    function(object,x,y,...) {
        log(sum(getNormWeights(object)*pnorm(-object@yp/sqrt(object@S),lower.tail=as.logical(y[1]))))
    }
)
