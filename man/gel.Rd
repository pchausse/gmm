\name{gel}

\alias{gel}
\alias{evalGel}

\title{Generalized Empirical Likelihood estimation}

\description{
Function to estimate a vector of parameters based on moment conditions using the GEL method as presented by Newey-Smith(2004) and Anatolyev(2005).
}
\usage{
gel(g, x, tet0 = NULL, gradv = NULL, smooth = FALSE,
    type = c("EL","ET","CUE","ETEL","HD","ETHD","RCUE"), 
    kernel = c("Truncated", "Bartlett"), bw = bwAndrews, 
    approx = c("AR(1)", "ARMA(1,1)"), prewhite = 1, ar.method = "ols", 
    tol_weights = 1e-7, tol_lam = 1e-9, tol_obj = 1e-9, tol_mom = 1e-9, 
    maxiterlam = 100, constraint = FALSE, optfct = c("optim", "optimize", 
    "nlminb"), optlam = c("nlminb", "optim", "iter", "Wu"), data,
    Lambdacontrol = list(), model = TRUE, X = FALSE, Y = FALSE,
    TypeGel = "baseGel", alpha = NULL, eqConst = NULL,
    eqConstFullVcov = FALSE, onlyCoefficients=FALSE, ...)
evalGel(g, x, tet0, gradv = NULL, smooth = FALSE,
        type = c("EL", "ET", "CUE", "ETEL", "HD", "ETHD","RCUE"),
        kernel = c("Truncated", "Bartlett"), bw = bwAndrews,
        approx = c("AR(1)", "ARMA(1,1)"), prewhite = 1,
        ar.method = "ols", tol_weights = 1e-7, tol_lam = 1e-9, tol_obj = 1e-9, 
        tol_mom = 1e-9, maxiterlam = 100, optlam = c("nlminb", "optim",
        "iter", "Wu"), data, Lambdacontrol = list(),
        model = TRUE, X = FALSE, Y = FALSE, alpha = NULL, ...)
}
\arguments{
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build the q sample moment conditions. It can also be a formula if the model is linear (see details below).  }

\item{tet0}{A \eqn{k \times 1} vector of starting values. If the
  dimension of \eqn{\theta} is one, see the argument "optfct". In the
  linear case, if tet0=NULL, the 2-step gmm estimator is used as starting
  value. However, it has to be provided when eqConst is not NULL}

\item{x}{The matrix or vector of data from which the function \eqn{g(\theta,x)} is computed. If "g" is a formula, it is an \eqn{n \times Nh} matrix of instruments (see details below).}

\item{gradv}{A function of the form \eqn{G(\theta,x)} which returns a \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)} with respect to \eqn{\theta}. By default, the numerical algorithm \code{numericDeriv} is used. It is of course strongly suggested to provide this function when it is possible. This gradiant is used compute the asymptotic covariance matrix of \eqn{\hat{\theta}}. If "g" is a formula, the gradiant is not required (see the details below).}

\item{smooth}{If set to TRUE, the moment function is smoothed as proposed by Kitamura(1997)}

\item{type}{"EL" for empirical likelihood, "ET" for exponential tilting,
  "CUE" for continuous updated estimator, "ETEL" for exponentially
  tilted empirical likelihood of Schennach(2007), "HD" for Hellinger
  Distance of Kitamura-Otsu-Evdokimov (2013), and "ETHD" for the
  exponentially tilted Hellinger distance of Antoine-Dovonon
  (2015). "RCUE" is a restricted version of "CUE" in which the
  probabilities are bounded below by zero. In that case, an analytical
  Kuhn-Tucker method is used to find the solution.}

\item{kernel}{type of kernel used to compute the covariance matrix of the vector of sample moment conditions (see \code{\link{kernHAC}} for more details) and to smooth the moment conditions if "smooth" is set to TRUE. Only two types of kernel are available. The truncated implies a Bartlett kernel for the HAC matrix and the Bartlett implies a Parzen kernel (see Smith 2004).}

\item{bw}{The method to compute the bandwidth parameter. By default it is \code{\link{bwAndrews}} which is proposed by Andrews (1991). The alternative is \code{\link{bwNeweyWest}} of Newey-West(1994).}

\item{prewhite}{logical or integer. Should the estimating functions be prewhitened? If \code{TRUE} or greater than 0 a VAR model of order \code{as.integer(prewhite)} is fitted via \code{ar} with method \code{"ols"} and \code{demean = FALSE}.}

\item{ar.method}{character. The \code{method} argument passed to \code{\link{ar}} for prewhitening.}

\item{approx}{a character specifying the approximation method if the bandwidth has to be chosen by \code{bwAndrews}.}

\item{tol_weights}{numeric. Weights that exceed \code{tol} are used for computing the covariance matrix, all other weights are treated as 0.}

\item{tol_lam}{Tolerance for \eqn{\lambda} between two iterations. The algorithm stops when \eqn{\|\lambda_i -\lambda_{i-1}\|} reaches \code{tol_lamb} (see \code{\link{getLamb}}) }

\item{maxiterlam}{The algorithm to compute \eqn{\lambda} stops if there is no convergence after "maxiterlam" iterations (see \code{\link{getLamb}}).}

\item{tol_obj}{Tolerance for the gradiant of the objective function to compute \eqn{\lambda} (see \code{\link{getLamb}}).}

\item{optfct}{Only when the dimension of \eqn{\theta} is 1, you can choose between the algorithm \code{\link{optim}} or \code{\link{optimize}}. In that case, the former is unreliable. If \code{\link{optimize}} is chosen, "t0" must be \eqn{1\times 2} which represents the interval in which the algorithm seeks the solution.It is also possible to choose the \code{\link{nlminb}} algorithm. In that case, borns for the coefficients can be set by the options \code{upper=} and \code{lower=}.}

\item{constraint}{If set to TRUE, the constraint optimization algorithm is used. See \code{\link{constrOptim}} to learn how it works. In particular, if you choose to use it, you need to provide "ui" and "ci" in order to impose the constraint \eqn{ui \theta - ci \geq 0}.}

\item{tol_mom}{It is the tolerance for the moment condition \eqn{\sum_{t=1}^n p_t g(\theta(x_t)=0}, where \eqn{p_t=\frac{1}{n}D\rho(<g_t,\lambda>)} is the implied probability. It adds a penalty if the solution diverges from its goal.}

\item{optlam}{Algorithm used to solve for the lagrange multiplier in
  \code{\link{getLamb}}. The algorithm Wu is only for
  \code{type="EL"}. The value of \code{optlam} is ignored for "CUE"
  because in that case, the analytical solution exists.}

\item{data}{A data.frame or a matrix with column names (Optional). }

\item{Lambdacontrol}{Controls for the optimization of the vector of Lagrange multipliers used by either \code{\link{optim}}, \code{\link{nlminb}} or \code{\link{constrOptim}}}

\item{model, X, Y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the response) are returned if g is a formula.}

\item{TypeGel}{The name of the class object created by the method \code{getModel}. It allows developers to extand the package and create other GEL methods.}

\item{alpha}{Regularization coefficient for discrete CGEL estimation (experimental). By setting alpha to any value, the model is estimated by CGEL of type specified    by the option \code{type}. See Chausse (2011)}

\item{eqConst}{Either a named vector (if "g" is a function), a simple vector for the nonlinear case indicating which of the \eqn{\theta_0} is restricted, or a qx2 vector defining equality constraints of the form \eqn{\theta_i=c_i}. See \code{\link{gmm}} for an example.}

\item{eqConstFullVcov}{If FALSE, the constrained coefficients are assumed to be fixed and only the covariance of the unconstrained coefficients is computed. If TRUE, the covariance matrix of the full set of coefficients is computed.}

\item{onlyCoefficients}{If \code{TRUE}, only the vector of coefficients
  and Lagrange multipliers are returned}

\item{...}{More options to give to \code{\link{optim}}, \code{\link{optimize}} or \code{\link{constrOptim}}.}

}

\details{
If we want to estimate a model like \eqn{Y_t = \theta_1 + X_{2t}\theta_2 + ... + X_{k}\theta_k + \epsilon_t} using the moment conditions \eqn{Cov(\epsilon_tH_t)=0}, where \eqn{H_t} is a vector of \eqn{Nh} instruments, than we can define "g" like we do for \code{\link{lm}}. We would have \code{g = y~x2+x3+...+xk} and the argument "x" above would become the matrix H of instruments. As for \code{\link{lm}}, \eqn{Y_t} can be a \eqn{Ny \times 1} vector which would imply that \eqn{k=Nh \times Ny}. The intercept is included by default so you do not have to add a column of ones to the matrix \eqn{H}. You do not need to provide the gradiant in that case since in that case it is embedded in \code{\link{gel}}. The intercept can be removed by adding -1 to the formula. In that case, the column of ones need to be added manually to H.

If "smooth" is set to TRUE, the sample moment conditions \eqn{\sum_{t=1}^n g(\theta,x_t)} is replaced by:
\eqn{\sum_{t=1}^n g^k(\theta,x_t)},
where
\eqn{g^k(\theta,x_t)=\sum_{i=-r}^r k(i) g(\theta,x_{t+i})},
where \eqn{r} is a truncated parameter that depends on the bandwidth and \eqn{k(i)} are normalized weights so that they sum to 1.

The method solves 
\eqn{\hat{\theta} = \arg\min \left[\arg\max_\lambda \frac{1}{n}\sum_{t=1}^n \rho(<g(\theta,x_t),\lambda>) - \rho(0) \right]}

\code{\link{evalGel}} generates the object of class "gel" for a fixed vector of
parameters. There is no estimation for \eqn{\theta}, but the optimal
vector of Lagrange multipliers \eqn{\lambda} is computed. The objective
function is then the profiled likelihood for a given \eqn{\theta}. It
can be used to construct a confidence interval by inverting  the
likelihood ratio test. 
}
 
\value{
'gel' returns an object of 'class' '"gel"' 

The functions 'summary' is used to obtain and print a summary of the results. 

The object of class "gel" is a list containing at least the following:

\item{coefficients}{\eqn{k\times 1} vector of parameters}

\item{residuals}{the residuals, that is response minus fitted values if "g" is a formula.}

\item{fitted.values}{the fitted mean values if "g" is a formula.}

\item{lambda}{\eqn{q \times 1} vector of Lagrange multipliers.}

\item{vcov_par}{the covariance matrix of "coefficients"}

\item{vcov_lambda}{the covariance matrix of "lambda"}

\item{pt}{The implied probabilities}

\item{objective}{the value of the objective function}

\item{conv_lambda}{Convergence code for "lambda" (see \code{\link{getLamb}})}

\item{conv_mes}{Convergence message for "lambda" (see \code{\link{getLamb}})}

\item{conv_par}{Convergence code for "coefficients" (see \code{\link{optim}}, \code{\link{optimize}} or \code{\link{constrOptim}})}

\item{terms}{the \code{\link{terms}} object used when g is a formula.}

\item{call}{the matched call.}
 
\item{y}{if requested, the response used (if "g" is a formula).}

\item{x}{if requested, the model matrix used if "g" is a formula or the data if "g" is a function.}

\item{model}{if requested (the default), the model frame used if "g" is a formula.}

}

\references{
Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Andrews DWK (1991), Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
\emph{Econometrica}, \bold{59}, 817--858.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Kitamura, Y. and Otsu, T. and Evdokimov, K. (2013), Robustness,
Infinitesimal Neighborhoods and Moment Restrictions.
\emph{Econometrica}, \bold{81}, 1185-1201.

Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

Smith, R.J. (2004), GEL Criteria for Moment Condition Models. \emph{Working paper, CEMMAP}.

Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
Heteroskedasticity and Autocorrelation Consistent Covariance
Matrix. \emph{Econometrica}, \bold{55}, 703--708.

Newey WK & West KD (1994), Automatic Lag Selection in Covariance
Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.

Schennach, Susanne, M. (2007), Point Estimation with Exponentially Tilted Empirical Likelihood.
\emph{Econometrica}, \bold{35}, 634-672.

Wu, C. (2005), Algorithms and R codes for the pseudo empirical
likelihood method in survey sampling.
\emph{Survey Methodology}, \bold{31}(2), page 239.

Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
\emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
URL \doi{10.18637/jss.v016.i09}.

Chausse (2010), Computing Generalized Method of Moments and Generalized Empirical Likelihood with R.
 \emph{Journal of Statistical Software}, \bold{34}(11), 1--35.
 URL \doi{10.18637/jss.v034.i11}.

Chausse (2011), Generalized Empirical likelihood for a continumm of moment conditions.
 \emph{Working Paper}, 
 \emph{Department of Economics},
 \emph{University of Waterloo}.

}
\examples{
# First, an exemple with the fonction g()

g <- function(tet, x)
	{
	n <- nrow(x)
	u <- (x[7:n] - tet[1] - tet[2]*x[6:(n-1)] - tet[3]*x[5:(n-2)])
	f <- cbind(u, u*x[4:(n-3)], u*x[3:(n-4)], u*x[2:(n-5)], u*x[1:(n-6)])
	return(f)
	}

Dg <- function(tet,x)
	{
	n <- nrow(x)
	xx <- cbind(rep(1, (n-6)), x[6:(n-1)], x[5:(n-2)])
        H  <- cbind(rep(1, (n-6)), x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
	f <- -crossprod(H, xx)/(n-6)
	return(f)
	}
n = 200
phi<-c(.2, .7)
thet <- 0.2
sd <- .2
set.seed(123)
x <- matrix(arima.sim(n = n, list(order = c(2, 0, 1), ar = phi, ma = thet, sd = sd)), ncol = 1)

res <- gel(g, x, c(0, .3, .6), grad = Dg)
summary(res)

# The same model but with g as a formula....  much simpler in that case

y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H

res <- gel(g, x, c(0, .3, .6))
summary(res)

# Using evalGel to create the object without estimation

res <- evalGel(g, x, res$coefficients)

}

