\name{gmm}

\alias{gmm}
\alias{gmmWithConst}
\alias{evalGmm}
	
\title{Generalized method of moment estimation}

\description{
 Function to estimate a vector of parameters based on moment conditions using the GMM method of Hansen(82). 
}
\usage{
gmm(g,x,t0=NULL,gradv=NULL, type=c("twoStep","cue","iterative"), 
    wmatrix = c("optimal","ident"), vcov=c("HAC","MDS","iid","TrueFixed"), 
    kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
    crit=10e-7,bw = bwAndrews, prewhite = 1, ar.method = "ols", approx="AR(1)",
    tol = 1e-7, itermax=100,optfct=c("optim","optimize","nlminb", "constrOptim"),
    model=TRUE, X=FALSE, Y=FALSE, TypeGmm = "baseGmm", centeredVcov = TRUE, 
    weightsMatrix = NULL, traceIter = FALSE, data, eqConst = NULL, 
    eqConstFullVcov = FALSE, mustar = NULL, onlyCoefficients=FALSE, ...)
evalGmm(g, x, t0, tetw=NULL, gradv=NULL, wmatrix = c("optimal","ident"),
    vcov=c("HAC","iid","TrueFixed"), kernel=c("Quadratic Spectral","Truncated", 
    "Bartlett", "Parzen", "Tukey-Hanning"),crit=10e-7,bw = bwAndrews,
    prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7,
    model=TRUE, X=FALSE, Y=FALSE,  centeredVcov = TRUE, weightsMatrix = NULL,
    data, mustar = NULL)
gmmWithConst(obj, which, value)
}
\arguments{
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build the q sample moment conditions. It can also be a formula if the model is linear (see details below).}

\item{x}{The matrix or vector of data from which the function \eqn{g(\theta,x)} is computed. If "g" is a formula, it is an \eqn{n \times Nh} matrix of instruments or a formula (see details below).}

\item{t0}{A \eqn{k \times 1} vector of starting values. It is required
  only when "g" is a function because only then a numerical algorithm is
  used to minimize the objective function. If the dimension of
  \eqn{\theta} is one, see the argument "optfct".}

\item{tetw}{A \eqn{k \times 1} vector to compute the weighting matrix.}

\item{gradv}{A function of the form \eqn{G(\theta,x)} which returns a \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)} with respect to \eqn{\theta}. By default, the numerical algorithm \code{numericDeriv} is used. It is of course strongly suggested to provide this function when it is possible. This gradient is used to compute the asymptotic covariance matrix of \eqn{\hat{\theta}} and to obtain the analytical gradient of the objective function if the method is set to "CG" or "BFGS" in \code{\link{optim}} and if "type" is not set to "cue". If "g" is a formula, the gradiant is not required (see the details below).}

\item{type}{The GMM method: "twostep" is the two step GMM proposed by Hansen(1982) and the "cue" and "iterative" are respectively the continuous updated and the iterative GMM proposed by Hansen, Eaton et Yaron (1996)}

\item{wmatrix}{Which weighting matrix should be used in the objective function. By default, it is the inverse of the covariance matrix of \eqn{g(\theta,x)}. The other choice is the identity matrix which is usually used to obtain a first step estimate of \eqn{\theta} }

\item{vcov}{Assumption on the properties of the random vector x. By default, x is a weakly dependant process. The "iid" option will avoid using the HAC matrix which will accelerate the estimation if one is ready to make that assumption. The option "TrueFixed" is used only when the matrix of weights is provided and it is the optimal one.}

\item{kernel}{type of kernel used to compute the covariance matrix of the vector of sample moment conditions (see \code{\link{kernHAC}} for more details)}

\item{crit}{The stopping rule for the iterative GMM. It can be reduce to increase the precision.}

\item{bw}{The method to compute the bandwidth parameter in the HAC
  weighting matrix. The default is \code{link{bwAndrews}} (as proposed in Andrews
  (1991)), which minimizes the MSE of the weighting matrix. Alternatives
  are \code{link{bwWilhelm}} (as proposed in Wilhelm
  (2015)), which minimizes the mean-square error (MSE) of the resulting
  GMM estimator, and \code{link{bwNeweyWest}} (as proposed in Newey-West(1994)).}

\item{prewhite}{logical or integer. Should the estimating functions be prewhitened? If \code{TRUE} or greater than 0 a VAR model of order \code{as.integer(prewhite)} is fitted via \code{ar} with method \code{"ols"} and \code{demean = FALSE}.}

\item{ar.method}{character. The \code{method} argument passed to \code{\link{ar}} for prewhitening.}

\item{approx}{A character specifying the approximation method if the bandwidth has to be chosen by \code{bwAndrews}.}

\item{tol}{Weights that exceed \code{tol} are used for computing the covariance matrix, all other weights are treated as 0.}

\item{itermax}{The maximum number of iterations for the iterative GMM. It is unlikely that the algorithm does not converge but we keep it as a safety.}

\item{optfct}{Only when the dimension of \eqn{\theta} is 1, you can choose between the algorithm \code{\link{optim}} or \code{\link{optimize}}. In that case, the former is unreliable. If \code{\link{optimize}} is chosen, "t0" must be \eqn{1\times 2} which represents the interval in which the algorithm seeks the solution. It is also possible to choose the \code{\link{nlminb}} algorithm. In that case, boundaries for the coefficients can be set by the options \code{upper=} and \code{lower=}. The \code{\link{constrOptim}} is only available for nonlinear models for now. The standard errors may have to be corrected if the estimtes reach the boundary set by ui and ci.}

\item{model, X, Y}{logical.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the response) are returned if g is a formula.}

\item{TypeGmm}{The name of the class object created by the method \code{getModel}. It allows developers to extend the package and create other GMM methods.}

\item{centeredVcov}{Should the moment function be centered when computing its covariance matrix. Doing so may improve inference.}

\item{weightsMatrix}{It allows users to provide \code{gmm} with a fixed weighting matrix. This matrix must be \eqn{q \times q}, symmetric and strictly positive definite. When provided, the \code{type} option becomes irrelevant. }

\item{traceIter}{Tracing information for GMM of type "iter"}

\item{data}{A data.frame or a matrix with column names (Optional). }

\item{eqConst}{Either a named vector (if "g" is a function), a simple vector for the nonlinear case indicating which of the \eqn{\theta_0} is restricted, or a qx2 vector defining equality constraints of the form \eqn{\theta_i=c_i}. See below for an example.}

\item{which, value}{The equality constraint is of the form which=value. "which" can be a vector of type characters with the names of the coefficients being constrained, or a vector of type numeric with the position of the coefficient in the whole vector.}

\item{obj}{Object of class "gmm"}

\item{eqConstFullVcov}{If FALSE, the constrained coefficients are assumed to be fixed and only the covariance of the unconstrained coefficients is computed. If TRUE, the covariance matrix of the full set of coefficients is computed.}

\item{mustar}{If not null, it must be a vector with the number of
  elements being equal to the number of moment conditions. In that case,
  the vector is subtracted from the sample moment vector before
  minimizing the objective function. It is useful to do a bootstrap
  procedure.
}

\item{onlyCoefficients}{If set to \code{TRUE}, the function only returns
  the coefficient estimates. It may be of interest when the standard
  errors are not needed}

\item{...}{More options to give to \code{\link{optim}}.}
}

\details{
If we want to estimate a model like \eqn{Y_t = \theta_1 + X_{2t} \theta_2 + \cdots + X_{k}\theta_k + \epsilon_t} using the moment conditions \eqn{Cov(\epsilon_tH_t)=0}, where \eqn{H_t} is a vector of \eqn{Nh} instruments, than we can define "g" like we do for \code{\link{lm}}. We would have \eqn{g = y ~\tilde{}~ x2+x3+ \cdots +xk} and the argument "x" above would become the matrix H of instruments. As for \code{\link{lm}}, \eqn{Y_t} can be a \eqn{Ny \times 1} vector which would imply that \eqn{k=Nh \times Ny}. The intercept is included by default so you do not have to add a column of ones to the matrix \eqn{H}. You do not need to provide the gradiant in that case since in that case it is embedded in \code{\link{gmm}}. The intercept can be removed by adding -1 to the formula. In that case, the column of ones need to be added manually to H. It is also possible to express "x" as a formula. For example, if the instruments are \eqn{\{1,z_1,z_2,z_3\}}, we can set "x" to \eqn{\tilde{} z1+z2+z3}. By default, a column of ones is added. To remove it, set "x" to \eqn{\tilde{}z1+z2+z3-1}. 

The following explains the last example bellow. Thanks to Dieter Rozenich, a student from the Vienna University of Economics and Business Administration. He suggested that it would help to understand the implementation of the Jacobian.  

For the two parameters of a normal distribution \eqn{(\mu,\sigma)} we have the following three moment conditions:
\deqn{
m_{1} = \mu - x_{i}
}
\deqn{ 
m_{2} = \sigma^2 - (x_{i}-\mu)^2
}
\deqn{ 
m_{3} = x_{i}^{3} - \mu (\mu^2+3\sigma^{2}) 
}
\eqn{m_{1},m_{2}} can be directly obtained by the definition of \eqn{(\mu,\sigma)}. The third moment condition comes from the third derivative of the moment generating function (MGF)

\deqn{
    M_{X}(t) = exp\Big(\mu t + \frac{\sigma^{2}t^{2}}{2}\Big)
}

evaluated at \eqn{(t=0)}.

Note that we have more equations (3) than unknown parameters (2).

The Jacobian of these two conditions is (it should be an array but I can't make it work):

\deqn{ 1~~~~~~~~~~ 0 }
\deqn{ -2\mu+2x ~~~~~ 2\sigma }
\deqn{-3\mu^{2}-3\sigma^{2} ~~~~ -6\mu\sigma}

\code{gmmWithConst()} re-estimates an unrestricted model by adding an
equality constraint.
\code{evalGmm()} creates an object of class '"gmm"' for a given
parameter vector. If no vector "tetw" is provided and the weighting
matrix needs to be computed, "t0" is used.,  
}

\value{
'gmm' returns an object of 'class' '"gmm"' 

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
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \doi{10.18637/jss.v016.i09}.

 Pierre Chausse (2010), Computing Generalized Method of Moments and Generalized Empirical Likelihood with R.
 \emph{Journal of Statistical Software}, \bold{34}(11), 1--35.
 URL \doi{10.18637/jss.v034.i11}.

 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

 Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
 Heteroskedasticity and Autocorrelation Consistent Covariance
 Matrix. \emph{Econometrica}, \bold{55}, 703--708.

 Newey WK & West KD (1994), Automatic Lag Selection in Covariance
 Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.

  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,

  Hansen, L.P. and Heaton, J. and Yaron, A.(1996),
  Finite-Sample Properties of Some Alternative GMM Estimators.
  \emph{Journal of Business and Economic Statistics}, \bold{14}
  262-280.
}
\examples{

## CAPM test with GMM
data(Finance)
r <- Finance[1:300, 1:10]
rm <- Finance[1:300, "rm"]
rf <- Finance[1:300, "rf"]

z <- as.matrix(r-rf)
t <- nrow(z)
zm <- rm-rf
h <- matrix(zm, t, 1)
res <- gmm(z ~ zm, x = h)
summary(res)

## linear tests can be performed using linearHypothesis from the car package
## The CAPM can be tested as follows:

library(car)
linearHypothesis(res,cbind(diag(10),matrix(0,10,10)),rep(0,10))

# The CAPM of Black
g <- function(theta, x) {
	e <- x[,2:11] - theta[1] - (x[,1] - theta[1]) \%*\% matrix(theta[2:11], 1, 10)
	gmat <- cbind(e, e*c(x[,1]))
	return(gmat) }

x <- as.matrix(cbind(rm, r))
res_black <- gmm(g, x = x, t0 = rep(0, 11))

summary(res_black)$coefficients


## APT test with Fama-French factors and GMM

f1 <- zm
f2 <- Finance[1:300, "hml"] 
f3 <- Finance[1:300, "smb"] 
h <- cbind(f1, f2, f3)
res2 <- gmm(z ~ f1 + f2 + f3, x = h)
coef(res2)
summary(res2)$coefficients

## Same result with x defined as a formula:

res2 <- gmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3)
coef(res2)

## The following example has been provided by Dieter Rozenich (see details).
# It generates normal random numbers and uses the GMM to estimate 
# mean and sd.
#-------------------------------------------------------------------------------
# Random numbers of a normal distribution
# First we generate normally distributed random numbers and compute the two parameters:
n <- 1000
x <- rnorm(n, mean = 4, sd = 2)
# Implementing the 3 moment conditions
g <- function(tet, x)
        {
        m1 <- (tet[1] - x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
        f <- cbind(m1, m2, m3)
        return(f)
        }
# Implementing the jacobian
Dg <- function(tet, x)
        {
        jacobian <- matrix(c( 1, 2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
			   -6*tet[1]*tet[2]), nrow=3,ncol=2)
        return(jacobian)
        }
# Now we want to estimate the two parameters using the GMM.
gmm(g, x, c(0, 0), grad = Dg)

# Two-stage-least-squares (2SLS), or IV with iid errors.
# The model is:
# Y(t) = b[0] + b[1]C(t) + b[2]Y(t-1) + e(t)
# e(t) is an MA(1)
# The instruments are Z(t)={1 C(t) y(t-2) y(t-3) y(t-4)}

getdat <- function(n) {
e <- arima.sim(n,model=list(ma=.9))
C <- runif(n,0,5)
Y <- rep(0,n)
Y[1] = 1 + 2*C[1] + e[1]
for (i in 2:n){
Y[i] = 1 + 2*C[i] + 0.9*Y[i-1] + e[i]
}
Yt <- Y[5:n]
X <- cbind(1,C[5:n],Y[4:(n-1)])
Z <- cbind(1,C[5:n],Y[3:(n-2)],Y[2:(n-3)],Y[1:(n-4)]) 
return(list(Y=Yt,X=X,Z=Z))
}

d <- getdat(5000)
res4 <- gmm(d$Y~d$X-1,~d$Z-1,vcov="iid")
res4

### Examples with equality constraint 
######################################

# Random numbers of a normal distribution

## Not run:
# The following works but produces warning message because the dimension of coef is 1
# Brent should be used

# without named vector
# Method Brent is used because the problem is now one-dimensional
gmm(g, x, c(4, 0), grad = Dg, eqConst=1, method="Brent", lower=-10,upper=10)
# with named vector
gmm(g, x, c(mu=4, sig=2), grad = Dg, eqConst="sig", method="Brent", lower=-10,upper=10)

## End(Not run)

gmm(g, x, c(4, 0), grad = Dg, eqConst=1,method="Brent",lower=0,upper=6)
gmm(g, x, c(mu=4, sig=2), grad = Dg, eqConst="sig",method="Brent",lower=0,upper=6)

# Example with formula
# first coef = 0 and second coef = 1
# Only available for one dimensional yt

z <- z[,1]
res2 <- gmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3, eqConst = matrix(c(1,2,0,1),2,2))
res2

# CUE with starting t0 requires eqConst to be a vector

res3 <- gmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3, t0=c(0,1,.5,.5), type="cue", eqConst = c(1,2))
res3

### Examples with equality constraints, where the constrained coefficients is used to compute
### the covariance matrix.
### Useful when some coefficients have been estimated before, they are just identified in GMM 
### and don't need to be re-estimated.
### To use with caution because the covariance won't be valid if the coefficients do not solve
### the GMM FOC.
######################################

res4 <- gmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3, t0=c(0,1,.5,.5), eqConst = c(1,2), 
  	    eqConstFullVcov=TRUE)
summary(res4)


### Examples with equality constraint using gmmWithConst
###########################################################

res2 <- gmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3)
gmmWithConst(res2,c("f2","f3"),c(.5,.5))
gmmWithConst(res2,c(2,3),c(.5,.5))

## Creating an object without estimation for a fixed parameter vector
###################################################################

res2_2 <- evalGmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3,
                  t0=res2$coefficients, tetw=res2$coefficients)
summary(res2_2)

}

