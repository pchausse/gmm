\name{confint}
\alias{confint.gel}
\alias{confint.ategel}
\alias{confint.gmm}
\alias{print.confint}
\title{Confidence intervals for GMM or GEL}
\description{
It produces confidence intervals for the coefficients from \code{gel} or \code{gmm} estimation.
}
\usage{
\method{confint}{gel}(object, parm, level = 0.95, lambda = FALSE,
                        type = c("Wald", "invLR", "invLM", "invJ"),
                        fact = 3, corr = NULL, ...)
\method{confint}{gmm}(object, parm, level = 0.95, ...)
\method{confint}{ategel}(object, parm, level = 0.95, lambda = FALSE,
                            type = c("Wald", "invLR", "invLM", "invJ"), fact = 3,
                            corr = NULL, robToMiss=TRUE, ...)
\method{print}{confint}(x, digits = 5, ...)
}
\arguments{
 \item{object}{An object of class \code{gel} or \code{gmm} returned by the function \code{\link{gel}} or \code{\link{gmm}}}
\item{parm}{A specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector
          of names.  If missing, all parameters are considered.}
\item{level}{The confidence level}
\item{lambda}{If set to TRUE, the confidence intervals for the Lagrange
  multipliers are produced.}
\item{type}{'Wald' is the usual symetric confidence interval. The thee
  others are based on the inversion of the LR, LM, and J tests.}
\item{fact}{This parameter control the span of search for the inversion
  of the test. By default we search within plus or minus 3 times the
  standard error of the coefficient estimate.}
\item{corr}{This numeric scalar is meant to apply a correction to the
  critical value, such as a Bartlett correction. This value depends on
  the model (See Owen; 2001)}
\item{x}{An object of class \code{confint} produced by \code{confint.gel} and \code{confint.gmm}}
\item{digits}{The number of digits to be printed}
\item{robToMiss}{If \code{TRUE}, the confidence interval is based on the
  standard errors that are robust to  misspecification}
\item{...}{Other arguments when \code{confint} is applied to another classe object}
}

\value{
It returns a matrix with the first column being the lower bound and the second the upper bound.} 


\references{
  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,
  Hansen, L.P. and Heaton, J. and Yaron, A.(1996),
  Finit-Sample Properties of Some Alternative GMM Estimators.
  \emph{Journal of Business and Economic Statistics}, \bold{14}
  262-280.
  Owen, A.B. (2001),
  Empirical Likelihood.
  \emph{Monographs on Statistics and Applied Probability 92, Chapman and
  Hall/CRC}
}


\examples{
#################
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2,0,1), ar = phi, ma = thet, sd = sd)), ncol = 1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H
t0 <- c(0,.5,.5)

resGel <- gel(g, x, t0)

confint(resGel)
confint(resGel, level = 0.90)
confint(resGel, lambda = TRUE)

########################

resGmm <- gmm(g, x)

confint(resGmm)
confint(resGmm, level = 0.90)

## Confidence interval with inversion of the LR, LM or J test.
##############################################################

set.seed(112233)
x <- rt(40, 3)
y <- x+rt(40,3)
# Simple interval on the mean
res <- gel(x~1, ~1, method="Brent", lower=-4, upper=4)
confint(res, type = "invLR")
confint(res)
# Using a Bartlett correction
k <- mean((x-mean(x))^4)/sd(x)^4
s <- mean((x-mean(x))^3)/sd(x)^3
a <- k/2-s^2/3
corr <- 1+a/40
confint(res, type = "invLR", corr=corr)

# Interval on the slope
res <- gel(y~x, ~x)
confint(res, "x", type="invLR")
confint(res, "x")
}

