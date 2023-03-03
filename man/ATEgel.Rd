\name{ATEgel}

\alias{ATEgel}
\alias{checkConv}

\title{ATE with Generalized Empirical Likelihood estimation}

\description{
Function to estimate the average treatment effect with the sample being
balanced by GEL.
}
\usage{
ATEgel(g, balm, w=NULL, y=NULL, treat=NULL, tet0=NULL,momType=c("bal","balSample","ATT"),
                   popMom = NULL, family=c("linear","logit", "probit"),
                   type = c("EL", "ET", "CUE", "ETEL", "HD", "ETHD", "RCUE"),
                   tol_lam = 1e-9, tol_obj = 1e-9, tol_mom = 1e-9, maxiterlam = 100,
                   optfct = c("optim", "nlminb"), 
                   optlam = c("nlminb", "optim", "iter", "Wu"), data=NULL,
                   Lambdacontrol = list(),
                   model = TRUE, X = FALSE, Y = FALSE, ...)
checkConv(obj, tolConv=1e-4, verbose=TRUE, ...)
}
\arguments{
\item{g}{A formula as \code{y~z}, where code{y} is the response and
  \code{z} the treatment indicator. If there is more than one
  treatment, more indicators can be added or \code{z} can be set as a
  factor. It can also be of the form
  \code{g(theta, y, z)} for non-linear models. It is however, not
  implemented yet.}

\item{obj}{Object of class \code{"ategel"} produced y \code{ATEgel}}

\item{balm}{A formula for the moments to be balanced between the treated
  and control groups (see details)}

\item{y}{The response variable when \code{g} is a function. Not
  implemented yet}

\item{treat}{The treatment indicator when \code{g} is a function. Not
  implemented yet}

\item{w}{A formula to add covariates to the main regression. When
  \code{NULL}, the default value, the main regression only include
  treatment indicators.}

\item{tet0}{A \eqn{3 \times 1} vector of starting values. If not
  provided, they are obtained using an OLS regression}

\item{momType}{How the moments of the covariates should be balanced. By
  default, it is simply balanced without restriction. Alternatively,
  moments can be set equal to the sample moments of the whole sample, or
  to the sample  moments of the treated group. The later will produce
  the average treatment effect of the treated (ATT)}

\item{popMom}{A vector of population moments to use for balancing. It
  can be used of those moments are available from a census, for
  example. When available, it greatly improves efficiency.}

\item{family}{By default, the outcome is linearly related to the
  treatment indicators. If the outcome is binary, it is possible to use
  the estimating equations of either the logit or probit model.}

\item{type}{"EL" for empirical likelihood, "ET" for exponential tilting,
  "CUE" for continuous updated estimator, "ETEL" for exponentially
  tilted empirical likelihood of Schennach(2007), "HD" for Hellinger
  Distance of Kitamura-Otsu-Evdokimov (2013), and "ETHD" for the
  exponentially tilted Hellinger distance of Antoine-Dovonon
  (2015). "RCUE" is a restricted version of "CUE" in which the
  probabilities are bounded below by zero. In that case, an analytical
  Kuhn-Tucker method is used to find the solution.}

\item{tol_lam}{Tolerance for \eqn{\lambda} between two iterations. The
  algorithm stops when \eqn{\|\lambda_i -\lambda_{i-1}\|} reaches
  \code{tol_lamb} (see \code{\link{getLamb}}) }

\item{maxiterlam}{The algorithm to compute \eqn{\lambda} stops if there
  is no convergence after "maxiterlam" iterations (see
  \code{\link{getLamb}}).}

\item{tol_obj}{Tolerance for the gradiant of the objective function to
  compute \eqn{\lambda} (see \code{\link{getLamb}}).}

\item{optfct}{Algorithm used for the parameter estimates}

\item{tol_mom}{It is the tolerance for the moment condition
  \eqn{\sum_{t=1}^n p_t g(\theta(x_t)=0}, where
  \eqn{p_t=\frac{1}{n}D\rho(<g_t,\lambda>)} is the implied probability. It
  adds a penalty if the solution diverges from its goal.}

\item{optlam}{Algorithm used to solve for the lagrange multiplier in
  \code{\link{getLamb}}. The algorithm Wu is only for
  \code{type="EL"}. The value of \code{optlam} is ignored for "CUE"
  because in that case, the analytical solution exists.}

\item{data}{A data.frame or a matrix with column names (Optional). }

\item{Lambdacontrol}{Controls for the optimization of the vector of
  Lagrange multipliers used by either \code{\link{optim}},
  \code{\link{nlminb}} or \code{\link{constrOptim}}}

\item{model, X, Y}{logicals.  If \code{TRUE} the corresponding
  components of the fit (the model frame, the model matrix, the response)
  are returned if g is a formula.}

\item{verbose}{If TRUE, a summary of the convergence is printed}

\item{tolConv}{The tolerance for comparing moments between groups}

\item{...}{More options to give to \code{\link{optim}} or
  \code{\link{nlminb}}. In \code{checkConv}, they are options passed to
  \code{\link{getImpProb}}.}

}

\details{
We want to estimate the model \eqn{Y_t = \theta_1 + \theta_2 treat +
  \epsilon_t}, where \eqn{\theta_2} is the treatment effect.  GEL is
used to balance the sample based on the argument \code{x} above.

For example, if we want the sample mean of \code{x1} and \code{x2} to be
balanced between the treated and control, we set \code{x} to
\code{~x1+x2}. If we want the sample mean of \code{x1}, \code{x2},
\code{x1*x2}, \code{x1^2} and \code{x2^2}, we set \code{x} to
\code{~x1*x2 + I(x1^2) + I(x2^2)}.
}
 
\value{
'gel' returns an object of 'class' '"ategel"'  

The functions 'summary' is used to obtain and print a summary of the results. 

The object of class "ategel" is a list containing the same elements
contained in objects of class \code{\link{gel}}.

}

\references{
Lee, Seojeong (2016), Asymptotic refinements of misspecified-robust
bootstrap for GEL estimators, \emph{Journal of Econometrics},
\bold{192}, 86--104.
  
Schennach, Susanne, M. (2007), Point Estimation with Exponentially Tilted Empirical Likelihood.
\emph{Econometrica}, \bold{35}, 634-672.

Wu, C. (2005), Algorithms and R codes for the pseudo empirical
likelihood method in survey sampling.
\emph{Survey Methodology}, \bold{31}(2), page 239.

Chausse, P. (2010), Computing Generalized Method of Moments and Generalized Empirical Likelihood with R.
 \emph{Journal of Statistical Software}, \bold{34}(11), 1--35.
 URL \doi{10.18637/jss.v034.i11}.

Chausse, P. and Giurcanu, M. and Luta, G. (2021) Estimating the Average Causal Effect
 using Generalized Empirical Likelihood Methods, Work in progress. 
}

\examples{
data(nsw)
# Scale income
nsw$re78 <- nsw$re78/1000
nsw$re75 <- nsw$re75/1000
res <- ATEgel(re78~treat, ~age+ed+black+hisp+married+nodeg+re75,
data=nsw,type="ET")
summary(res)
chk <- checkConv(res)

res2 <- ATEgel(re78~treat, ~age+ed+black+hisp+married+nodeg+re75,
data=nsw,type="ET", momType="balSample")
summary(res2)
chk2 <- checkConv(res2)

}
