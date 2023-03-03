\name{FinRes}
\alias{FinRes.baseGmm.res}
\title{Method to finalize the result of the momentEstim method}
\description{
It computes the final results that will be needed to create the object of class \code{gmm}.). 
}
\usage{
\method{FinRes}{baseGmm.res}(z, object, ...)
}
\arguments{
\item{z}{An object of class determined by the method \code{momentEstim}.}
\item{object}{An object produced my \code{getModel}}
\item{...}{Other argument to be passed to other \code{FinRes} methods.}
}

\value{
It returns an object of class \code{gmm}. See \code{\link{gmm}} for more details. 
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
}



