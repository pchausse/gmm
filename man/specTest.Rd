\name{specTest}
\alias{specTest}
\alias{specTest.gel}
\alias{specTest.gmm}
\alias{print.specTest}

\title{Compute tests of specification}
\description{
 Generic function for testing the specification of estimated models. It computes the J-test from \code{gmm} objects and J-test, LR-test and LM-test from 
\code{gel} objects.
}
\usage{
\method{specTest}{gmm}(x, ...)
\method{specTest}{gel}(x, ...)
\method{print}{specTest}(x, digits = 5, ...)
specTest(x, \dots)
}
\arguments{
 \item{x}{A fitted model object.}
 \item{digits}{The number of digits to be printed.}
 \item{\dots}{Arguments passed to methods.}
}

\value{
Tests and p-values
}

\references{

 Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,

 Smith, R. J. (2004),
 GEL Criteria for Moment Condition Models.
 \emph{CeMMAP working papers, Institute for Fiscal Studies}
 
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
specTest(res)

###################
res <- gmm(g, x)
specTest(res)

}

