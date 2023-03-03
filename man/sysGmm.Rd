\name{sysGmm}

\alias{sysGmm}
\alias{five}
\alias{sur}
\alias{randEffect}
\alias{threeSLS}
	
\title{Generalized method of moment estimation for system of equations}

\description{
Functions to estimate a system of equations based on GMM. 
}
\usage{
sysGmm(g, h, wmatrix = c("optimal","ident"),
vcov=c("MDS", "HAC", "CondHom", "TrueFixed"),
       kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
       crit=10e-7,bw = bwAndrews, prewhite = FALSE, ar.method = "ols", approx="AR(1)",
       tol = 1e-7, model=TRUE, X=FALSE, Y=FALSE, centeredVcov = TRUE,
       weightsMatrix = NULL, data, crossEquConst = NULL, commonCoef = FALSE)
five(g, h, commonCoef = FALSE, data = NULL)
threeSLS(g, h, commonCoef = FALSE, data = NULL)
sur(g, commonCoef = FALSE, data = NULL)
randEffect(g, data = NULL)
}
\arguments{
\item{g}{A possibly named list of formulas}

\item{h}{A formula if the same instruments are used in each equation or
  a list of formulas.}

\item{wmatrix}{Which weighting matrix should be used in the objective function. By default, it is the inverse of the covariance matrix of \eqn{g(\theta,x)}. The other choice is the identity matrix.}

\item{vcov}{Assumption on the properties of the moment vector. By
  default, it is a martingale difference sequence. "HAC" is for weakly
  dependent processes and "CondHom" implies conditional
  homoscedasticity. The option "TrueFixed" is used only when the matrix of weights is provided and it is the optimal one.}

\item{kernel}{type of kernel used to compute the covariance matrix of the vector of sample moment conditions (see \code{\link{kernHAC}} for more details)}

\item{crit}{The stopping rule for the iterative GMM. It can be reduce to increase the precision.}

\item{bw}{The method to compute the bandwidth parameter. By default it is \code{\link{bwAndrews}} which is proposed by Andrews (1991). The alternative is \code{\link{bwNeweyWest}} of Newey-West(1994).}

\item{prewhite}{logical or integer. Should the estimating functions be prewhitened? If \code{TRUE} or greater than 0 a VAR model of order \code{as.integer(prewhite)} is fitted via \code{ar} with method \code{"ols"} and \code{demean = FALSE}.}

\item{ar.method}{character. The \code{method} argument passed to \code{\link{ar}} for prewhitening.}

\item{approx}{A character specifying the approximation method if the bandwidth has to be chosen by \code{bwAndrews}.}

\item{tol}{Weights that exceed \code{tol} are used for computing the covariance matrix, all other weights are treated as 0.}

\item{model, X, Y}{logical.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the response) are returned if g is a formula.}

\item{centeredVcov}{Should the moment function be centered when computing its covariance matrix. Doing so may improve inference.}

\item{weightsMatrix}{It allows users to provide \code{gmm} with a fixed weighting matrix. This matrix must be \eqn{q \times q}, symmetric and strictly positive definite. When provided, the \code{type} option becomes irrelevant. }

\item{data}{A data.frame or a matrix with column names (Optional). }

\item{commonCoef}{If true, coefficients accross equations are the same}

\item{crossEquConst}{Only used if the number of regressors are the same
  in each equation. It is a vector which indicates which coefficient are
  constant across equations. The order is 1 for Intercept and 2 to k as
  it is formulated in the formulas \code{g}. Setting it to \code{1:k} is
  equivalent to setting commonCoef to TRUE.}

}

\details{

This set of functions implement the estimation of system of equations as presented in Hayashi (2000)
  
}

\value{
'sysGmm' returns an object of 'class' '"sysGmm"' 

The functions 'summary' is used to obtain and print a summary of the results.  It also compute the J-test of overidentying restriction

The object of class "sysGmm" is a list containing at least:

\item{coefficients}{list of vectors of coefficients for each equation}

\item{residuals}{list of the residuals for each equation.}

\item{fitted.values}{list of the fitted values for each equation.}

\item{vcov}{the covariance matrix of the stacked coefficients}

\item{objective}{the value of the objective function \eqn{\| var(\bar{g})^{-1/2}\bar{g}\|^2}}

\item{terms}{The list of \code{\link{terms}} objects for each equation}

\item{call}{the matched call.}
 
\item{y}{If requested, a list of response variables.}

\item{x}{if requested, a list of the model matrices.}

\item{model}{if requested (the default), a list of the model frames.}

 }


\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \doi{10.18637/jss.v016.i09}.

 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

 Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
 Heteroskedasticity and Autocorrelation Consistent Covariance
 Matrix. \emph{Econometrica}, \bold{55}, 703--708.

 Newey WK & West KD (1994), Automatic Lag Selection in Covariance
 Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.

  Hayashi, F. (2000), Econometrics.
  \emph{Princeton University Press}.
}
\examples{

data(wage)

eq1 <- LW~S+IQ+EXPR
eq2 <- LW80~S80+IQ+EXPR80
g2 <- list(Wage69=eq1, WAGE80=eq2)
h2 <- list(~S+EXPR+MED+KWW, ~S80+EXPR80+MED+KWW)

res <- sysGmm(g2, h2, data=wage, commonCoef=TRUE)
summary(res)

res2 <- sysGmm(g2, h2, data=wage)
summary(res2)

five(g2, h2, data=wage)

threeSLS(g2, h2[[1]], data=wage)

sur(g2, data=wage)

randEffect(g2, data=wage)

## Cross-Equation restrictions
## All but the intercept are assumed to be the same

res <- sysGmm(g2, h2, data=wage, crossEquConst = 2:4)
summary(res)

}

