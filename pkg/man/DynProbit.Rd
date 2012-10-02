\name{DynProbit}
\alias{DynProbit-class}
\alias{DynProbit}
\alias{DynProbitInit}

\alias{particles,DynProbit-method}
\alias{getParticles,DynProbit-method}
\alias{setParticles,DynProbit-method}
\alias{ParticleMove,DynProbit-method}
\alias{doParticleMove,DynProbit-method}
\alias{SmcIterate,DynProbit-method}
\alias{doSmcIterate,DynProbit-method}
\alias{UpdateWeights,DynProbit-method}
\alias{doUpdateWeights,DynProbit-method}
\alias{ESS,DynProbit-method}
\alias{getLogWeights,DynProbit-method}
\alias{getWeights,DynProbit-method}
\alias{getNormWeights,DynProbit-method}

\title{Dynamic Probit regression}

\description{A class to implement a Sequential Monte Carlo sampler for a dynamic Probit regression model.}

\usage{
    DynProbit(formula,mu0,Sigma0,Q,N = 2000,data,resampleC = .5,
        output=c("mean","filter","smoothing"))
}

\arguments{
	
	\item{formula}{A formula with the response and covariates in the model. See Details.}
	
	\item{mu0}{Mean vector of the prior distribution over the Probit regression parameters.}
	
	\item{Sigma0}{Covariance matrix of the prior distribution over the Probit regression parameters.}
	
	\item{Q}{Covariance matrix of the update distribution over the Probit regression parameters.}
	
	\item{N}{Total number of particles.}
	
	\item{data}{An optional \code{data.frame} to interpret the variables in
		the formula.}
	
	\item{resampleC}{value between 0 and 1, indicating when to resample}
	
	\item{output}{specifying what is returned. See details.}
		
}

\details{
    
    The function DynProbit initializes a particle filter for a dynamic Probit regression model, and runs it over the 
    data specified in the formula. The response can be either a factor, in which the first level is taken as ``failure''
    and the remaining ones as ``success'', or a binary vector.
    
    The prior distriution over the parameters at $t=1$ is a (multivariate) Normal distribution with mean \code{mu0} and covariance \code{Sigma0}.
    
    If \code{output = means}, a matrix with the posterior means at each time point is returned.
    
    If \code{output = filter}, an particle filter object of class \code{DynProbit} is returned. The filter is in the state at t=T.
    
    If \code{output = smoother}, a list is returned with the particle paths of the N particles and the weights at t=T. Note that this
    is a naive implementation of a particle smoother, and the results are not expected to be good for initial time points. This
    option is to be used mainly for didacting reasons.
}

\author{Maarten Speekenbrink}

\examples{
    set.seed(120912)
    nt <- 100
    x <- cbind(rnorm(nt,sd=1),rnorm(nt,sd=1))
    theta <- cbind(
        rnorm(1,mean=0,sd=sqrt(.2))+c(0,cumsum(rnorm(nt-1,sd=sqrt(.01)))),
        rnorm(1,mean=0,sd=sqrt(.2))+c(0,cumsum(rnorm(nt-1,sd=sqrt(.01)))),
        rnorm(1,mean=0,sd=sqrt(.2))+c(0,cumsum(rnorm(nt-1,sd=sqrt(.01))))
    )
    y <- rbinom(nt,prob=1-pnorm(rowSums(cbind(1,x)*theta)),size=1)
    
    means <- DynProbit(y~x[,1] + x[,2],mu0=c(0,0,0),Sigma0=.2*diag(3),Q=.01*diag(3),N=2000)
    
    \notrun{
    layout(matrix(1:3,ncol=3))
    plot(c(1,100),c(min(c(means[,1],theta[,1])),max(c(means[,1],theta[,1]))),type="n")
    lines(1:100,theta[,1],lty=2)
    lines(1:100,means[,1],type="l",lty=1)
    plot(c(1,100),c(min(c(means[,2],theta[,2])),max(c(means[,2],theta[,2]))),type="n")
    lines(1:100,theta[,2],lty=2)
    lines(1:100,means[,2],type="l",lty=1)
    plot(c(1,100),c(min(c(means[,3],theta[,3])),max(c(means[,3],theta[,3]))),type="n")
    lines(1:100,theta[,3],lty=2)
    lines(1:100,means[,3],type="l",lty=1)
    }


}
