\name{getLamb}

\alias{getLamb}

\title{Solving for the Lagrange multipliers of Generalized Empirical Likelihood (GEL) }

\description{
 It computes the vector of Lagrange multipliers, which maximizes the GEL objective function, using an iterative Newton method.
}
\usage{
getLamb(gt, l0, type = c("EL","ET","CUE", "ETEL", "HD","ETHD","RCUE"),
        tol_lam = 1e-7, maxiterlam = 100, 
	tol_obj = 1e-7, k = 1, method = c("nlminb", "optim", "iter", "Wu"),
        control = list())
}
\arguments{
\item{gt}{A \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)}}

\item{l0}{Vector of starting values for lambda}

\item{type}{"EL" for empirical likelihood, "ET" for exponential tilting,
  "CUE" for continuous updated estimator, and "HD" for Hellinger
  Distance. See details for "ETEL" and "ETHD". "RCUE" is a restricted
  version of "CUE" in which the probabilities are bounded below by
  zero. In that case, an analytical Kuhn-Tucker method is used to find
  the solution.}

\item{tol_lam}{Tolerance for \eqn{\lambda} between two iterations. The
algorithm stops when \eqn{\|\lambda_i -\lambda_{i-1}\|} reaches
\code{tol_lam} }

\item{maxiterlam}{The algorithm stops if there is no convergence after
"maxiterlam" iterations.}

\item{tol_obj}{Tolerance for the gradiant of the objective function. The
algorithm returns a non-convergence message if \eqn{\max(|gradiant|)}
does not reach \code{tol_obj}. It helps the \code{gel} algorithm to
select the right space to look for \eqn{\theta}}

\item{k}{It represents the ratio k1/k2, where
\eqn{k1=\int_{-\infty}^{\infty} k(s)ds} and
\eqn{k2=\int_{-\infty}^{\infty} k(s)^2 ds}. See Smith(2004).}

\item{method}{The iterative procedure uses a Newton method for solving
  the FOC. It i however recommended to use \code{optim} or
  \code{nlminb}. If type is set to "EL" and method to "optim",
  \code{\link{constrOptim}} is called to prevent \eqn{log(1-gt'\lambda)}
  from producing NA. The gradient and hessian is provided to
  \code{nlminb} which speed up the convergence. The latter is therefore
  the default value. "Wu" is for "EL" only. It uses the algorithm of Wu
  (2005). The value of \code{method} is ignored for "CUE" because in
  that case, the analytical solution exists.}

\item{control}{Controls to send to \code{\link{optim}},
\code{\link{nlminb}} or \code{\link{constrOptim}}} } \details{ It solves
the problem \eqn{\max_{\lambda} \frac{1}{n}\sum_{t=1}^n
\rho(gt'\lambda)}. For the type "ETEL", it is only used by
\code{\link{gel}}. In that case \eqn{\lambda} is obtained by maximizing
\eqn{\frac{1}{n}\sum_{t=1}^n \rho(gt'\lambda)}, using
\eqn{\rho(v)=-\exp{v}} (so ET) and \eqn{\theta} by minimizing the same
equation but with \eqn{\rho(v)-\log{(1-v)}}. To avoid NA's,
\code{\link{constrOptim}} is used with the restriction \eqn{\lambda'g_t
< 1}. The type "ETHD" is experimental and proposed by Antoine-Dovonon
(2015). The paper is not yet available.  }

\value{
lambda: A \eqn{q\times 1} vector of Lagrange multipliers which solve the system of equations given above.
\code{conv}: Details on the type of convergence.
 }

\references{
Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

Smith, R.J. (2004), GEL Criteria for Moment Condition Models. \emph{Working paper, CEMMAP}.

Wu, C. (2005), Algorithms and R codes for the pseudo empirical
likelihood method in survey sampling.
\emph{Survey Methodology}, \bold{31}(2), page 239.

}

\examples{
g <- function(tet,x)
	{
	n <- nrow(x)
	u <- (x[7:n] - tet[1] - tet[2]*x[6:(n-1)] - tet[3]*x[5:(n-2)])
	f <- cbind(u, u*x[4:(n-3)], u*x[3:(n-4)], u*x[2:(n-5)], u*x[1:(n-6)])
	return(f)
	}
n = 500
phi<-c(.2, .7)
thet <- 0.2
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2, 0, 1), ar = phi, ma = thet, sd = sd)), ncol = 1)
gt <- g(c(0,phi),x)
getLamb(gt, type = "EL",method="optim")
}


