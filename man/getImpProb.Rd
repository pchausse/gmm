\name{getImpProb}
\alias{getImpProb}
\alias{getImpProb.gel}

\title{Implied Probabilities} \description{ It computes the implied
probabilities from objects of class \code{gel} with additional options.
} 

\usage{ 
\method{getImpProb}{gel}(object, posProb=TRUE, normalize=TRUE,
                         checkConv=FALSE,...)
}



\arguments{ \item{object}{Object of class
\code{gel}.}  \item{posProb}{Should the implied probabilities be
transformed into positive probabilities?}  \item{normalize}{Should we
normalize the probabilities so that they sum to one?}
\item{checkConv}{Should we add the attribute convergence to check the
  sum of the probabilities and the weighted sum of the moment conditions?}
\item{...}{Additional arguments to pass to other methods}
}

\value{
  A vector af implied probabilities.
}

\references{
  Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and
  Generalized Empirical Likelihood Estimators. \emph{Econometrica},
  \bold{72}, 219-255. 
}



\examples{

#################
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H
t0 <- c(0,.5,.5)

res <- gel(g, x, t0)
pt <- getImpProb(res)
}

