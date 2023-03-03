\name{estfun}
\alias{estfun.gmmFct}
\alias{estfun.gmm}
\alias{estfun.gel}
\alias{estfun.tsls}
\alias{model.matrix.tsls}
\title{Extracts the empirical moment function}
\description{
It extracts the matrix of empirical moments so that it can be used by the \code{\link{kernHAC}} function. 
}
\usage{
\method{estfun}{gmmFct}(x, y = NULL, theta = NULL, ...)
\method{estfun}{gmm}(x, ...)
\method{estfun}{gel}(x, ...)
\method{estfun}{tsls}(x, ...)
\method{model.matrix}{tsls}(object, ...)
}
\arguments{
\item{x}{A function of the form \eqn{g(\theta,y)} or a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,y_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n} or an object of class \code{gmm}. See \code{\link{gmm}} for more details. For \code{\link{tsls}}, it is an object of class \code{tsls}.}
\item{object}{An object of class \code{tsls}.}
\item{y}{The matrix or vector of data from which the function \eqn{g(\theta,y)} is computed if \code{g} is a function.}
\item{theta}{Vector of parameters if \code{g} is a function.}
\item{...}{Other arguments when \code{estfun} is applied to another class object}
}

\details{
For \code{estfun.gmmFct}, it returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,y_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. It is only used by \code{gmm} to obtain the estimates.

For \code{estfun.gmm}, it returns the matrix of first order conditions of \eqn{\min_\theta \bar{g}'W\bar{g}/2}, which is a \eqn{n \times k} matrix with the \eqn{t^{th}} row being \eqn{g(\theta, y_t)W G}, where \eqn{G} is \eqn{d\bar{g}/d\theta}. It allows to compute the sandwich covariance matrix using \code{\link{kernHAC}} or \code{\link{vcovHAC}} when \eqn{W} is not the optimal matrix.

The method if not yet available for \code{gel} objects.

For tsls, model.matrix and estfun are used by \code{vcov()} to compute different covariance matrices using the \code{\link{sandwich}} package. See \code{\link{vcov.tsls}}. \code{model.matrix} returns the fitted values frin the first stage regression and \code{esfun} the residuals.
}

\value{
A \eqn{n \times q} matrix (see details). 
}

\references{
   Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \doi{10.18637/jss.v016.i09}.

}

\examples{

n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]
H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H
res <- gmm(g, x,weightsMatrix = diag(5))

gt <- res$gt
G <- res$G

foc <- gt%*%G
foc2 <- estfun(res)

foc[1:5,]
foc2[1:5,]

}




