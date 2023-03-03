\name{KTest}
\alias{KTest}
\alias{print.gmmTests}

\title{Compute the K statistics of Kleibergen}
\description{The test is proposed by Kleibergen (2005). It is robust to weak identification.
}
\usage{
KTest(obj, theta0 = NULL, alphaK = 0.04, alphaJ = 0.01)
\method{print}{gmmTests}(x, digits = 5, ...)
}
\arguments{
 \item{obj}{Object of class "gmm" returned by \code{\link{gmm}}}
 \item{theta0}{The null hypothesis being tested. See details.}
 \item{alphaK, alphaJ}{The size of the J and K tests when combining the two. The overall size is alphaK+alphaJ.}
 \item{x}{An object of class \code{gmmTests} returned by \code{KTest}}
 \item{digits}{The number of digits to be printed}
 \item{...}{Other arguments when \code{print} is applied to another class object}
}

\details{
The function produces the J-test and K-statistics which are robust to weak identification. The test is either \eqn{H0:\theta=theta_0}, in which case theta0 must be provided, or \eqn{\beta=\beta_0}, where \eqn{\theta=(\alpha', \beta')'}, and \eqn{\alpha} is assumed to be identified. In the latter case, theta0 is NULL and obj is a restricted estimation in which \eqn{\beta} is fixed to \eqn{\beta_0}. See \code{\link{gmm}} and the option "eqConst" for more details.  
}

\value{
Tests and p-values
}

\references{

 Keibergen, F. (2005),
  Testing Parameters in GMM without assuming that they are identified.
  \emph{Econometrica}, \bold{73},
  1103-1123,
 
}

\examples{

library(mvtnorm)
sig <- matrix(c(1,.5,.5,1),2,2)
n <- 400
e <- rmvnorm(n,sigma=sig)
x4 <- rnorm(n)
w <- exp(-x4^2) + e[,1]
y <- 0.1*w + e[,2]
h <- cbind(x4, x4^2, x4^3, x4^6)
g3 <- y~w
res <- gmm(g3,h)

# Testing the whole vector:

KTest(res,theta0=c(0,.1))

# Testing a subset of the vector (See \code{\link{gmm}})

res2 <- gmm(g3, h, eqConst=matrix(c(2,.1),1,2))
res2
KTest(res2)

}
