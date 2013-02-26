setClass("ParticleBase",
	representation(
	    particles = "ANY",
		logWeights = "vector", # log of particle weights
		unifWeights = "logical", # are current weights uniform?
		p_move = "function",
		mcmc_move = "function",
		lW_update = "function",
		logLik = "function",
		resampleC = "numeric",
		N = "integer"
	)
)

setMethod("particles",signature(object="ParticleBase"),
    function(object,...) {
        return(object@particles)
    }
)

setReplaceMethod(
    f = "particles",
    signature = "ParticleBase",
    definition = function(object,value) {
        object@particles <- value
        return(object)
    }
)


setMethod("ParticleMove",signature(object="ParticleBase"),
	function(object,...) {
		object@particles <- object@p_move(particles(object,...),...)
		object
	}
)

setMethod("McmcMove",signature(object="ParticleBase"),
	function(object,...) {
		object@particles <- object@mcmc_move(particles(object,...),...)
		object
	}
)

setMethod("SmcIterate",signature(object="ParticleBase"),
	function(object,...) {
		object <- ParticleMove(object,...)
		object <- UpdateWeights(object,...)
		ess <- ESS(object,...)
		if(ess < (object@resampleC*object@N)) object <- Resample(object,...)
		object
	}
)

setMethod("UpdateWeights",signature(object="ParticleBase"),
	function(object,data,...) {
		logWeights <- object@lW_update(getParticles(object,...),logWeights(object,...),...)
		# normalize to sensible values
		max <- max(-.Machine$double.xmax,logWeights)
		object@logWeights <- logWeights - max
		
		#dMaxWeight = -std::numeric_limits<double>::infinity();
    # for(int i = 0; i < N; i++)
    #   dMaxWeight = max(dMaxWeight, pParticles[i].GetLogWeight());
    # for(int i = 0; i < N; i++)
    #   pParticles[i].SetLogWeight(pParticles[i].GetLogWeight() - (dMaxWeight));  
		
		object@unifWeights <- FALSE
		object
	}
)


setMethod("ESS",signature(object="ParticleBase"),
	function(object,...) {
		#w <- exp(object@logWeights)
		sumw <- sum(getWeights(object,...))
		sumsq <- sum(exp(2*logWeights(object,...)))
		exp(-log(sumsq) + 2*log(sumw))
	}
)

setMethod("logWeights",signature(object="ParticleBase"),
  function(object,...) {
    return(object@logWeights) 
  }
)

setReplaceMethod(
    f = "logWeights",
    signature = "ParticleBase",
    definition = function(object,value) {
        # normalize to sensible values
		max <- max(-.Machine$double.xmax,value)
		object@logWeights <- value - max
		if(length(unique((value))) == 1) object@unifWeights <- TRUE
        return(object)
    }
)




setMethod("getWeights",signature(object="ParticleBase"),
  function(object,...) {
    if(!object@unifWeights) return(exp(logWeights(object,...))) else return(rep(1,object@N))
  }
)

setMethod("getNormWeights",signature(object="ParticleBase"),
  function(object,...) {
    if(!object@unifWeights) {
        weights <- getWeights(object,...)
        weights <- weights/sum(weights)
        return(weights)
    } else {
        return(rep(1/object@N,object@N))
    }
  }
)

setMethod("getN",signature(object="ParticleBase"),
    function(object,...) {
        return(object@N)
    }
)

setMethod("logLik",signature(object="ParticleBase"),
    function(object,...) {
        return(object@logLik(particles(object,...),...))
    }
)

