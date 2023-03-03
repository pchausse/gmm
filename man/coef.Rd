\name{coef}
\alias{coef.gel}
\alias{coef.gmm}
\title{Coefficients of GEL or GMM}
\description{
 It extracts the coefficients from \code{gel} or \code{gmm} objects.
}
\usage{
\method{coef}{gmm}(object, ...)
\method{coef}{gel}(object, lambda = FALSE, ...)
}
\arguments{
 \item{object}{An object of class \code{gel} or \code{gmm} returned by the function \code{\link{gel}} or \code{\link{gmm}}}
\item{lambda}{If set to TRUE, the lagrange multipliers are extracted instead of the vector of coefficients}
\item{...}{Other arguments when \code{coef} is applied to an other class object}
}

\value{
Vector of coefficients
}


\examples{

#################
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
t0 <- c(0,.5,.5)

res <- gel(g, x, t0)

coef(res)
coef(res, lambda = TRUE)
###################
res <- gmm(g, x)
coef(res)

}

