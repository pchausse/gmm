\name{tsls}

\alias{tsls}

\title{Two stage least squares estimation}

\description{
 Function to estimate a linear model by the two stage least squares method.
}
\usage{
tsls(g,x,data)
}
\arguments{
\item{g}{A formula describing the linear regression model (see details below).}

\item{x}{The matrix of instruments (see details below).}

\item{data}{A data.frame or a matrix with column names (Optionnal). }
}

\details{
The function just calls \code{\link{gmm}} with the option vcov="iid". It just simplifies the the implementation of 2SLS. The users don't have to worry about all the options offered in \code{\link{gmm}}. The model is 
\deqn{
Y_i = X_i\beta + u_i
}
In the first step, \code{\link{lm}}  is used to regress \eqn{X_i} on the set of instruments \eqn{Z_i}. The second step also uses \code{\link{lm}} to regress \eqn{Y_i} on the fitted values of the first step. 
}

\value{
'tsls' returns an object of 'class' '"tsls"' which inherits from  class '"gmm"'.

The functions 'summary' is used to obtain and print a summary of the results.  It also compute the J-test of overidentying restriction

The object of class "gmm" is a list containing at least:

\item{coefficients}{\eqn{k\times 1} vector of coefficients}

\item{residuals}{the residuals, that is response minus fitted values if "g" is a formula.}

\item{fitted.values}{the fitted mean values if "g" is a formula.}

\item{vcov}{the covariance matrix of the coefficients}

\item{objective}{the value of the objective function \eqn{\| var(\bar{g})^{-1/2}\bar{g}\|^2}}

\item{terms}{the \code{\link{terms}} object used when g is a formula.}

\item{call}{the matched call.}
 
\item{y}{if requested, the response used (if "g" is a formula).}

\item{x}{if requested, the model matrix used if "g" is a formula or the data if "g" is a function.}

\item{model}{if requested (the default), the model frame used if "g" is a formula.}

\item{algoInfo}{Information produced by either \code{\link{optim}} or \code{\link{nlminb}} related to the convergence if "g" is a function. It is printed by the \code{summary.gmm} method.}
}


\references{

  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,

}
\examples{

n <- 1000
e <- arima.sim(n,model=list(ma=.9))
C <- runif(n,0,5)
Y <- rep(0,n)
Y[1] = 1 + 2*C[1] + e[1]
for (i in 2:n){
Y[i] = 1 + 2*C[i] + 0.9*Y[i-1] + e[i]
}
Yt <- Y[5:n]
X <- cbind(C[5:n],Y[4:(n-1)])
Z <- cbind(C[5:n],Y[3:(n-2)],Y[2:(n-3)],Y[1:(n-4)]) 

res <- tsls(Yt~X,~Z)
res

}

