\name{getModel}
\alias{getModel.baseGmm}
\alias{getModel.sysGmm}
\alias{getModel.baseGel}
\alias{getModel.constGmm}
\alias{getModel.tsls}
\alias{getModel.constGel}
\alias{getModel.ateGel}
\title{Method for setting the properties of a model}
\description{
It collects what is needed by the method \code{momentEstim} (see details). 
}
\usage{
\method{getModel}{baseGmm}(object, ...)
\method{getModel}{sysGmm}(object, ...)
\method{getModel}{baseGel}(object, ...)
\method{getModel}{constGel}(object, ...)
\method{getModel}{constGel}(object, ...)
\method{getModel}{tsls}(object, ...)
\method{getModel}{ateGel}(object, ...)
}
\arguments{
\item{object}{An object of class \code{baseGmm} }
\item{...}{Other arguments when \code{getModel} is applied to another class object}
}

\value{
It returns an object of the right class which determines how the method \code{momentEstim} will treat it. For example, if \code{g} is a formula and \code{type} is set to "cue", it creates an object of class \code{baseGmm.cue.formula}. It this case, \code{momentEstim}, applied to this object, computes the continuously updated GMM of a linear model. It allows more flexibility this way. For example, it could be easy to add a GMM method which is robust in presence of weak identification simply by creating a new class of model and the associated \code{momentEstime} method.
}



