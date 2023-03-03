\name{bread}
\alias{bread.gmm}
\alias{bread.gel}
\alias{bread.tsls}
\title{Bread for sandwiches}
\description{
Computes the bread of the sandwich covariance matrix
}
\usage{
\method{bread}{gmm}(x, ...)
\method{bread}{gel}(x, ...)
\method{bread}{tsls}(x, ...)
}
\arguments{
\item{x}{A fitted model of class \code{gmm} or \code{gel}.}
\item{...}{Other arguments when \code{bread} is applied to another class object}
}

\details{
When the weighting matrix is not the optimal one, the covariance matrix of the estimated coefficients is:
\eqn{(G'WG)^{-1} G'W V W G(G'WG)^{-1}}, 
where \eqn{G=d\bar{g}/d\theta}, \eqn{W} is the matrix of weights, and \eqn{V} is the covariance matrix of the moment function. Therefore, the bread is \eqn{(G'WG)^{-1}}, which is the second derivative of the objective function. 

The method if not yet available for \code{gel} objects.
}

\value{
A \eqn{k \times k} matrix (see details). 
}

\references{
   Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \doi{10.18637/jss.v016.i09}.

}

\examples{
# See \code{\link{gmm}} for more details on this example.
# With the identity matrix 
# bread is the inverse of (G'G)

n <- 1000
x <- rnorm(n, mean = 4, sd = 2)
g <- function(tet, x)
        {
        m1 <- (tet[1] - x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
        f <- cbind(m1, m2, m3)
        return(f)
        }
Dg <- function(tet, x)
        {
        jacobian <- matrix(c( 1, 2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
				-6*tet[1]*tet[2]), nrow=3,ncol=2)
        return(jacobian)
        }

res <- gmm(g, x, c(0, 0), grad = Dg,weightsMatrix=diag(3))
G <- Dg(res$coef, x)
bread(res)
solve(crossprod(G))
}












