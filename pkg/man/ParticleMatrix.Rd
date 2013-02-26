\name{ParticleMatrix-class}

\docType{class}

\alias{ParticleMatrix-class}
\alias{ParticleCov}
\alias{ParticleMean}
\alias{ParticleVar} 

\alias{particles,ParticleMatrix-method}
\alias{particles<-,ParticleMatrix-method}
\alias{ParticleMove,ParticleMatrix-method}
\alias{SmcIterate,ParticleMatrix-method}
\alias{UpdateWeights,ParticleMatrix-method}
\alias{ESS,ParticleMatrix-method}
\alias{logWeights,ParticleMatrix-method}
\alias{logWeights<-,ParticleMatrix-method}
\alias{getWeights,ParticleMatrix-method}
\alias{getNormWeights,ParticleMatrix-method}
\alias{ParticleMean,ParticleMatrix-method}
\alias{ParticleCov,ParticleMatrix-method}
\alias{ParticleVar,ParticleMatrix-method}
\alias{mean,ParticleMatrix-method}
\alias{cov,ParticleMatrix-method}
\alias{var,ParticleMatrix-method}
\alias{Resample,ParticleMatrix-method}


\title{Class "ParticleMatrix"}

\description{Basic class to implement a Sequential Monte Carlo sampler.}

\section{Slots}{

  Code: N lW_update logLik logWeights mcmc_move p_margin p_move
        particles resampleC unifWeights
  Inherited: particles logWeights unifWeights p_move mcmc_move
        lW_update logLik resampleC N
  Docs: Particles code{p_margin}: lW_update logWeights mcmc_move p_move
        resampleC unifWeights


	\describe{
    \item{\code{particles}:}{Matrix with particle values}
    
    \item{\code{p_margin}:}{Integer indicating which margin in the matrix represents particle ids}

	\item{\code{logWeights}:}{Vector containing the log (unnormalised) particle weights.}

	\item{\code{unifWeights}}{Logical indicating whether the logWeights are uniform.}

	\item{\code{p_move}:}{Function to move the particles to a new position.}
	
	\item{\code{mcmc_move}:}{Function to perform a Monte Carlo Markov Chain move.}

	\item{\code{lW_update}:}{Function to do update the \code{logWeights}.}

	\item{\code{resampleC}:}{Numeric value (between 0 and 1) indicating when to perform resampling.}

	}
}

\section{Methods}{
\describe{
    \item{particles}{\code{signature(object = "ParticleMatrix")}: return particles}
    \item{particles<-}{\code{signature(object = "ParticleMatrix")}: set particles}
    \item{ParticleMove}{\code{signature(object = "ParticleMatrix")}: move particles}
    #\item{doParticleMove}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{SmcIterate}{\code{signature(object = "ParticleMatrix")}: do a full SMC iteration}
    #\item{doSmcIterate}{\code{signature(object = "ParticleMatrix")}: move particles}
    
    \item{UpdateWeights}{\code{signature(object = "ParticleMatrix")}: calculate updated weights}
    #\item{doUpdateWeights}{\code{signature(object = "ParticleMatrix")}: move particles}
    \item{ESS}{\code{signature(object = "ParticleMatrix")}: effective sample size}
    \item{logWeights}{\code{signature(object = "ParticleMatrix")}: return log weights}
    \item{logWeights<-}{\code{signature(object = "ParticleMatrix")}: set log weights}
    \item{getWeights}{\code{signature(object = "ParticleMatrix")}: return the unnormalized weights}
    \item{getNormWeights}{\code{signature(object = "ParticleMatrix")}: return the self-normalized weights}
    \item{mean}{\code{signature(object = "ParticleMatrix")}: return the estimated mean of each particle dimension}
    \item{var}{\code{signature(object = "ParticleMatrix")}: return the estimated variance of each particle dimension}
    \item{cov}{\code{signature(object = "ParticleMatrix")}: return the estimated covariance matrix of the particle dimensions}

  }

}

\author{Maarten Speekenbrink}

\keyword{classes}
