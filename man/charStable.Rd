\name{charStable}
\alias{charStable}
\title{The characteristic function of a stable distribution}
\description{
 It computes the theoretical characteristic function of a stable distribution for two different parametrizations. It is used in the vignette to illustrate the estimation of the parameters using GMM.}
\usage{
charStable(theta, tau, pm = 0)
}
\arguments{
 \item{theta}{Vector of parameters of the stable distribution. See details.}
\item{tau}{A vector of numbers at which the function is evaluated.}
\item{pm}{The type of parametization. It takes the values 0 or 1.}
}

\value{
It returns a vector of complex numbers with the dimension equals to \code{length(tau)}.
}

\details{
The function returns the vector \eqn{\Psi(\theta,\tau,pm)} defined as \eqn{E(e^{ix\tau}}, where \eqn{\tau} is a vector of real numbers, \eqn{i} is the imaginary number, \eqn{x} is a stable random variable with parameters \eqn{\theta} = \eqn{(\alpha,\beta,\gamma,\delta)} and \code{pm} is the type of parametrization. The vector of parameters are the characteristic exponent, the skewness, the scale and the location parameters, respectively. The restrictions on the parameters are: \eqn{\alpha \in (0,2]}, \eqn{\beta\in [-1,1]} and \eqn{\gamma>0}. For mode details see Nolan(2009).
}


\references{
Nolan J. P. (2020), Univariate Stable Distributions - Models for Heavy Tailed Data.
  \emph{Springer Series in Operations Research and Financial Engineering}.
  URL \url{https://edspace.american.edu/jpnolan/stable/}.
}

\examples{

# GMM is like GLS for linear models without endogeneity problems

pm <- 0
theta <- c(1.5,.5,1,0) 
tau <- seq(-3, 3, length.out = 20)
char_fct <- charStable(theta, tau, pm)

}

