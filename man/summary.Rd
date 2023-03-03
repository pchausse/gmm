\name{summary}
\alias{summary.gmm}
\alias{summary.sysGmm}
\alias{summary.gel}
\alias{summary.ategel}
\alias{summary.tsls}
\alias{print.summary.gmm}
\alias{print.summary.sysGmm}
\alias{print.summary.gel}
\alias{print.summary.tsls}
\title{Method for object of class gmm or gel}
\description{
It presents the results from the \code{gmm} or \code{gel} estimation in the same fashion as \code{summary} does for the \code{lm} class objects for example. It also compute the tests for overidentifying restrictions.
}
\usage{
\method{summary}{gmm}(object, ...)
\method{summary}{sysGmm}(object, ...)
\method{summary}{gel}(object, ...)
\method{summary}{ategel}(object, robToMiss = TRUE, ...)
\method{summary}{tsls}(object, vcov = NULL, ...)
\method{print}{summary.gmm}(x, digits = 5, ...)
\method{print}{summary.sysGmm}(x, digits = 5, ...)
\method{print}{summary.gel}(x, digits = 5, ...)
\method{print}{summary.tsls}(x, digits = 5, ...)
}
\arguments{
 \item{object}{An object of class \code{gmm} or \code{gel} returned by the function \code{\link{gmm}} or \code{\link{gel}}}
 \item{x}{An object of class \code{summary.gmm} or \code{summary.gel} returned by the function \code{\link{summary.gmm}} \code{\link{summary.gel}}}
\item{digits}{The number of digits to be printed}
\item{vcov}{An alternative covariance matrix computed with
  \code{vcov.tsls}}
\item{robToMiss}{If \code{TRUE}, it computes the robust to
    misspecification covariance matrix}
\item{...}{Other arguments when summary is applied to another class object}
}

\value{
It returns a list with the parameter estimates and their standard deviations, t-stat and p-values. It also returns the J-test and p-value for the null hypothesis that \eqn{E(g(\theta,X)=0} 
}

\references{
  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,

  Hansen, L.P. and Heaton, J. and Yaron, A.(1996),
  Finit-Sample Properties of Some Alternative GMM Estimators.
  \emph{Journal of Business and Economic Statistics}, \bold{14}
  262-280.
	
Anatolyev, S. (2005), 
GMM, GEL, Serial Correlation, and Asymptotic Bias. 
\emph{Econometrica}, \bold{73}, 
983-1002.

Kitamura, Yuichi (1997), 
Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25},
 2084-2102.

Newey, W.K. and Smith, R.J. (2004), 
Higher Order Properties of GMM and Generalized Empirical Likelihood Estimators. 
\emph{Econometrica}, \bold{72}, 
219-255.
}


\examples{

# GMM #
set.seed(444)
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2,0,1), ar = phi, ma = thet, sd = sd)), ncol = 1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]
ym3 <- x[4:(n-3)]
ym4 <- x[3:(n-4)]
ym5 <- x[2:(n-5)]
ym6 <- x[1:(n-6)]

g <- y ~ ym1 + ym2
x <- ~ym3+ym4+ym5+ym6

res <- gmm(g, x)

summary(res)

# GEL #

t0 <- res$coef
res <- gel(g, x, t0)
summary(res)


# tsls #

res <- tsls(y ~ ym1 + ym2,~ym3+ym4+ym5+ym6)
summary(res)

}

