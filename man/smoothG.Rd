\name{smoothG}
\alias{smoothG}
\title{Kernel smoothing of a matrix of time series}
\description{
 It applies the required kernel smoothing to the moment function in order for the GEL estimator to be valid. It is used by the \code{gel} function.}
\usage{
smoothG(x, bw = bwAndrews, prewhite = 1, ar.method = "ols", weights = weightsAndrews,
	kernel = c("Bartlett", "Parzen", "Truncated", "Tukey-Hanning"), 
	approx = c("AR(1)", "ARMA(1,1)"), tol = 1e-7) 
}
\arguments{
 \item{x}{a \eqn{n\times q} matrix of time series, where n is the sample size.}
 \item{bw}{The method to compute the bandwidth parameter. By default, it uses the bandwidth proposed by Andrews(1991). As an alternative, we can choose bw=bwNeweyWest (without "") which is proposed by Newey-West(1996).}
 \item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}.}
\item{ar.method}{character. The \code{method} argument passed to
   \code{\link{ar}} for prewhitening.}
\item{weights}{The smoothing weights can be computed by \code{\link{weightsAndrews}} of it can be provided manually. If provided, it has to be a \eqn{r\times 1}vector (see details). }
\item{approx}{a character specifying the approximation method if the
    bandwidth has to be chosen by \code{bwAndrews}.}
\item{tol}{numeric. Weights that exceed \code{tol} are used for computing
   the covariance matrix, all other weights are treated as 0.}
\item{kernel}{The choice of kernel}
}


\details{
The sample moment conditions \eqn{\sum_{t=1}^n g(\theta,x_t)} is replaced by:
\eqn{\sum_{t=1}^n g^k(\theta,x_t)}, where \eqn{g^k(\theta,x_t)=\sum_{i=-r}^r k(i) g(\theta,x_{t+i})},
where \eqn{r} is a truncated parameter that depends on the bandwidth and \eqn{k(i)} are normalized weights so that they sum to 1.

If the vector of weights is provided, it gives only one side weights. For exemple, if you provide the vector (1,.5,.25), \eqn{k(i)} will become \eqn{(.25,.5,1,.5,.25)/(.25+.5+1+.5+.25) =  (.1,.2,.4,.2,.1)}
}

\value{
smoothx: A \eqn{q \times q} matrix containing an estimator of the asymptotic variance of \eqn{\sqrt{n} \bar{x}}, where \eqn{\bar{x}} is \eqn{q\times 1}vector with typical element \eqn{\bar{x}_i = \frac{1}{n}\sum_{j=1}^nx_{ji}}. This function is called by \code{\link{gel}} but can also be used by itself.

\code{kern_weights}: Vector of weights used for the smoothing.
}

\references{
Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Andrews DWK (1991),
Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
\emph{Econometrica}, \bold{59},
817--858.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
\emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
URL \doi{10.18637/jss.v016.i09}.
}

\examples{
g <- function(tet, x)
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
gt <- g(c(0, phi), x) 
sgt <- smoothG(gt)$smoothx
plot(gt[,1])
lines(sgt[,1])
}


