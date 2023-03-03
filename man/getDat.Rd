\name{getDat}
\alias{getDat}
\title{Extracting data from a formula}
\description{
It extract the data from a formula y~z with instrument h and put everything in a matrix. It helps redefine the function \eqn{g(\theta,x)} that is required by \code{\link{gmm}} and \code{\link{gel}}.
}
\usage{
getDat(formula, h, data, error=TRUE) 
}
\arguments{
\item{formula}{A formula that defines the linear model to be estimated (see details).}
\item{h}{A \eqn{n\times nh} matrix of intruments(see details).}
\item{data}{A data.frame or a matrix with colnames (Optionnal).}
\item{error}{If FALSE, the data is generated without giving any error message}
}
\details{The model to be estimated is based on the moment conditions
  \eqn{<h,(y-z\theta)>=0}. It adds a column of ones to z and h by
  default. They are removed if -1 is added to the formula.
  The error argument has been added for \code{\link{sysGmm}} with common
  coefficients because the check is only valid for
  equation by equation identification. 
}
\value{
x: A \eqn{n \times l} matrix, where \eqn{l = ncol(y)+ncol(z)+ncol(h)+2} if "intercept" is TRUE and \eqn{ncol(y)+ncol(z)+xcol(h)} if "intercept" is FALSE. 

nh: dimension of h

k: dimension of z

ny: dimension of y
}

\examples{
n = 500
phi<-c(.2, .7)
thet <- 0.2
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2, 0, 1), ar = phi, ma = thet, sd = sd)), ncol = 1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]
H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])

x <- getDat(y ~ ym1 + ym2, H)

}

