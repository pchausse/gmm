\name{momentEstim}
\alias{momentEstim.baseGmm.twoStep}
\alias{momentEstim.baseGmm.twoStep.formula}
\alias{momentEstim.sysGmm.twoStep.formula}
\alias{momentEstim.tsls.twoStep.formula}
\alias{momentEstim.baseGmm.iterative.formula}
\alias{momentEstim.baseGmm.iterative}
\alias{momentEstim.baseGmm.cue.formula}
\alias{momentEstim.baseGmm.cue}
\alias{momentEstim.baseGmm.eval}
\alias{momentEstim.baseGel.mod}
\alias{momentEstim.baseGel.modFormula}
\alias{momentEstim.baseGel.eval}
\title{Method for estimating models based on moment conditions}
\description{
It estimates a model which is caracterized by the method \code{getModel} (see details). 
}
\usage{
\method{momentEstim}{baseGmm.twoStep}(object, ...)
\method{momentEstim}{baseGmm.twoStep.formula}(object, ...)
\method{momentEstim}{sysGmm.twoStep.formula}(object, ...)
\method{momentEstim}{tsls.twoStep.formula}(object, ...)
\method{momentEstim}{baseGmm.iterative.formula}(object, ...)
\method{momentEstim}{baseGmm.iterative}(object, ...)
\method{momentEstim}{baseGmm.cue.formula}(object, ...)
\method{momentEstim}{baseGmm.cue}(object, ...)
\method{momentEstim}{baseGmm.eval}(object, ...)
\method{momentEstim}{baseGel.mod}(object, ...)
\method{momentEstim}{baseGel.modFormula}(object, ...)
\method{momentEstim}{baseGel.eval}(object, ...)
}
\arguments{
\item{object}{An object created by the method \code{getModel}}
\item{...}{Other arguments when \code{momentEstim} is applied to an other class object}
}

\value{
It returns an object of class determined by the argument "TypeGMM" of \code{\link{gmm}}. By default, it is of class \code{baseGmm.res}. It estimates the model and organize the results that will be finalized by the method \code{FinRes}. More methods can be created in order to use other GMM methods not yet included in the package. 
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



