\name{ParticleBase-class}

\docType{class}

\alias{ParticleBase-class}

\alias{particles}
\alias{McmcMove}
\alias{ParticleMove}
\alias{Resample}
\alias{SmcIterate}
\alias{UpdateWeights}
\alias{doMcmcMove}
\alias{doParticleMove}
\alias{doResample}
\alias{doSmcIterate}
\alias{doUpdateWeights}
\alias{getLogWeights}
\alias{getN}
\alias{getNormWeights}
\alias{getParticles}
\alias{getWeights}
\alias{setLogWeights<-}
\alias{setParticles<-}

\alias{particles,ParticleBase-method}
\alias{particles<-,ParticleBase-method}
\alias{getParticles,ParticleBase-method}
\alias{setParticles,ParticleBase-method}
\alias{ParticleMove,ParticleBase-method}
\alias{doParticleMove,ParticleBase-method}
\alias{SmcIterate,ParticleBase-method}
\alias{doSmcIterate,ParticleBase-method}
\alias{UpdateWeights,ParticleBase-method}
\alias{doUpdateWeights,ParticleBase-method}
\alias{ESS}
\alias{ESS,ParticleBase-method}
\alias{getLogWeights,ParticleBase-method}
\alias{setLogWeights<-,ParticleBase-method}
\alias{getWeights,ParticleBase-method}
\alias{getNormWeights,ParticleBase-method}

\alias{McmcMove,ParticleBase-method}
\alias{doMcmcMove,ParticleBase-method}
\alias{getN,ParticleBase-method}
\alias{setLogWeights,ParticleBase-method}
\alias{setParticles<-,ParticleBase-method} 

\title{Class "ParticleBase"}

\description{Basic class to implement a Sequential Monte Carlo sampler.}

\section{Slots}{

	\describe{
        particles = "ANY",
		logWeights = "vector", # log of particle weights
		unifWeights = "logical", # are current weights uniform?
		p_move = "function",
		lW_update = "function",
		logLik = "function",
		resampleC = "numeric",
		N = "integer"

    \item{\code{particles}:}{Set of particles. Format depends on implementation.}
    
	\item{\code{logWeights}:}{Vector containing the log (unnormalised) particle weights.}

	\item{\code{unifWeights}}{Logical indicating whether the logWeights are uniform.}

	\item{\code{p_move}:}{Function to move the particles to a new position.}
	
	\item{\code{mcmc_move}:}{Function to perform a Monte Carlo Markov Chain move.}

	\item{\code{lW_update}:}{Function to do update the \code{logWeights}.}
	
	\item{\code{logLik}:}{Function to compute the \code{logLik}.}

	\item{\code{resampleC}:}{Numeric value (between 0 and 1) indicating when to perform resampling.}
	
	\item{\code{N}:}{Integer indicating the total number of particles.}

	}
}

\section{Methods}{
\describe{
    \item{particles}{\code{signature(object = "ParticleBase")}: return particles}
    \item{particles<-}{\code{signature(object = "ParticleBase")}: set particles}
    \item{ParticleMove}{\code{signature(object = "ParticleBase")}: move particles}
    #\item{doParticleMove}{\code{signature(object = "ParticleBase")}: move particles}
    \item{SmcIterate}{\code{signature(object = "ParticleBase")}: perform a full SMC iteration}
    #\item{doSmcIterate}{\code{signature(object = "ParticleBase")}: move particles}
    
    \item{UpdateWeights}{\code{signature(object = "ParticleBase")}: update weights}
    #\item{doUpdateWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{ESS}{\code{signature(object = "ParticleBase")}: Effective Sample Size}
    \item{logWeights}{\code{signature(object = "ParticleBase")}: get the log of the particles weights}
    \item{logWeights<-}{\code{signature(object = "ParticleBase")}: set the log of the particles weights}
    #\item{getLogWeights}{\code{signature(object = "ParticleBase")}: move particles}
    #\item{setLogWeights}{\code{signature(object = "ParticleBase")}: move particles}
    \item{getWeights}{\code{signature(object = "ParticleBase")}: get the exponentiated (unnormalized) weights}
    \item{getNormWeights}{\code{signature(object = "ParticleBase")}: get the normalized weights}
  }
}

\author{Maarten Speekenbrink}

\keyword{classes}
