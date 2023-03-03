\name{marginal}
\alias{marginal}
\alias{marginal.ategel}
\title{Marginal effects Summary}
\description{
It produces the summary table of marginal effects for GLM estimation with
GEL. Only implemented for ATEgel.
}
\usage{
\method{marginal}{ategel}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{ategel} returned by the function
    \code{\link{ATEgel}}}
  \item{...}{Other arguments for other methods}
}

\value{
It returns a matrix with the marginal effects, the standard errors based on the Delta
method when the link is nonlinear, the t-ratios, and the pvalues.} 


\references{
  Owen, A.B. (2001),
  Empirical Likelihood.
  \emph{Monographs on Statistics and Applied Probability 92, Chapman and
  Hall/CRC}
}


\examples{
## We create some artificial data with unbalanced groups and binary outcome
genDat <- function(n)
    {
        eta=c(-1, .5, -.25, -.1)
        Z <- matrix(rnorm(n*4),ncol=4)
        b <- c(27.4, 13.7, 13.7, 13.7)
        bZ <- c(Z\%*\%b)
        Y1 <- as.numeric(rnorm(n, mean=210+bZ)>220)
        Y0 <- as.numeric(rnorm(n, mean=200-.5*bZ)>220)
        etaZ <- c(Z\%*\%eta)
        pZ <- exp(etaZ)/(1+exp(etaZ))
        T <- rbinom(n, 1, pZ)
        Y <- T*Y1+(1-T)*Y0
        X1 <- exp(Z[,1]/2)
        X2 <- Z[,2]/(1+exp(Z[,1]))
        X3 <- (Z[,1]*Z[,3]/25+0.6)^3
        X4 <- (Z[,2]+Z[,4]+20)^2
        data.frame(Y=Y, cbind(X1,X2,X3,X4), T=T)
    }

dat <- genDat(200)
res <- ATEgel(Y~T, ~X1+X2+X3+X4, data=dat, type="ET", family="logit")
summary(res)

marginal(res)

}

