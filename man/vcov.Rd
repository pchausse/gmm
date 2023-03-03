\name{vcov}
\alias{vcov.gmm}
\alias{vcov.gel}
\alias{vcov.tsls}
\alias{vcov.ategel}
\title{Variance-covariance matrix of GMM or GEL}
\description{
 It extracts the matrix of variances and covariances from \code{gmm} or \code{gel} objects.
}
\usage{
\method{vcov}{gmm}(object, ...)
\method{vcov}{gel}(object, lambda = FALSE, ...)
\method{vcov}{tsls}(object, type=c("Classical","HC0","HC1","HAC"),
                    hacProp = list(), ...)
\method{vcov}{ategel}(object, lambda = FALSE, robToMiss = TRUE, ...)
}

\arguments{
 \item{object}{An object of class \code{gmm} or \code{gmm} returned by the function \code{\link{gmm}} or \code{\link{gel}}}
\item{lambda}{If set to TRUE, the covariance matrix of the Lagrange multipliers is produced.}
\item{type}{Type of covariance matrix for the meat}
\item{hacProp}{A list of arguments to pass to \code{\link{kernHAC}}}
\item{robToMiss}{If \code{TRUE}, it computes the robust to
 misspecification covariance matrix}
\item{...}{Other arguments when \code{vcov} is applied to another class object}
}

\details{
For tsls(), if vcov is set to a different value thand "Classical", a sandwich covariance matrix is computed.

}

\value{
A matrix of variances and covariances
}


\examples{

# GMM #
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n = n,list(order = c(2,0,1), ar = phi, ma = thet, sd = sd)), ncol = 1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H

res <- gmm(g, x)
vcov(res)

## GEL ##

t0 <- c(0,.5,.5)
res <- gel(g, x, t0)
vcov(res)
vcov(res, lambda = TRUE)

}

