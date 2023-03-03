\name{bwWilhelm}
\alias{bwWilhelm}
	
\title{Wilhelm (2015) bandwidth selection}

\description{
It computes the optimal bandwidth for the HAC estimation of the
covariance matrix of the moment conditions. The bandwidth was shown by
Wilhelm (2005) to be the one that minimizes the MSE of the GMM estimator. 
}
\usage{
bwWilhelm(x, order.by = NULL, kernel = c("Quadratic Spectral",
  "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
  weights = NULL, prewhite = 1, ar.method = "ols", data = list())
}
\arguments{
\item{x}{An object of class \code{gmm}.}

\item{order.by}{Either a vector 'z' or a formula with a single
  explanatory variable like '~ z'. The observations in the model are ordered by the size of 'z'. If set to 'NULL' (the default)
the observations are assumed to be ordered (e.g., a time
series).}

\item{kernel}{type of kernel used to compute the covariance matrix of the vector of sample moment conditions (see \code{\link{kernHAC}} for more details)}

\item{approx}{A character specifying the approximation method if the bandwidth has to be chosen by \code{bwAndrews}.}

\item{weights}{numeric. A vector of weights used for weighting the
  estimated coefficients of the approximation model (as specified by 'approx'). By default all weights are 1 except that for the
  intercept term (if there is more than one variable)}

\item{prewhite}{logical or integer. Should the estimating functions be prewhitened? If \code{TRUE} or greater than 0 a VAR model of order \code{as.integer(prewhite)} is fitted via \code{ar} with method \code{"ols"} and \code{demean = FALSE}.}

\item{ar.method}{character. The \code{method} argument passed to \code{\link{ar}} for prewhitening.}

\item{data}{an optional data frame containing the variables in the 'order.by' model.}

}

\value{
The function 'bwWilhelm'  returns the optimal bandwidth.
 }


\references{

  Wilhelm, D. (2015), Optimal Bandwidth Selection for Robust Generalized
  Method of Moments Estimation.
  \emph{Econometric Theory}, \bold{31}, 1054--1077

  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \doi{10.18637/jss.v016.i09}.
}

\note{
  The function was written by Daniel Wilhelm and is based on \link{bwAndrews}.
}

\examples{
data(Finance)
f1 <- Finance[1:300, "rm"] 
f2 <- Finance[1:300, "hml"] 
f3 <- Finance[1:300, "smb"] 
y <- Finance[1:300,"WMK"]

## Silly example just to make it over-identified
###############################################
res <- gmm(y ~ f1, ~ f1 + f2 + f3)
summary(res)


## Set the bandwidth using the second step estimate
################################################
bw <- bwWilhelm(res)
res2 <- update(res, bw=bw)
summary(res2)

## Set the bandwidth using the first-step estimate as for bwAndrews
###################################################################
res3 <- gmm(y ~ f1, ~ f1 + f2 + f3, bw=bwWilhelm)
summary(res3)




}

