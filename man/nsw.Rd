\name{nsw}
\docType{data}
\alias{nsw}
\title{Lalonde subsample of the National Supported Work Demonstration Data
(NSW)}
\description{
This data was collected to evaluate the National Supported Work
     (NSW) Demonstration project in Lalonde (1986).
}
\usage{data(nsw)}
\format{
A data frame containing 9 variables. 
\describe{
\item{treat}{Treatment assignment} 
\item{age}{Age} 
\item{ed}{Years of Education} 
\item{black}{1 if Black, 0 otherwise} 
\item{hisp}{1 if Hispanic 0 otherwise} 
\item{married}{1 if married 0 otherwise} 
\item{nodeg}{1 if no college degree 0 otherwise} 
\item{re75}{1975 earnings} 
\item{re78}{1978 earnings}
}
}

\details{
  The dataset was obtained from the ATE package (see reference).
  }
  

\source{
"NSW Data Files" from Rajeev Dehejia's website. URL:
\url{http://users.nber.org/~rdehejia/data/.nswdata2.html}

"National Supported Work Evaluation Study, 1975-1979: Public Use Files."
from the Interuniversity Consortium for Political and Social
Research. URL: \url{http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/7865}
}

\references{

  Lalonde, R. (1986). "Evaluating the Econometric Evaluations of Training Programs," American Economic Review, 76(4), 604-620.

  Dehejia R. and Wahba S. (1999).  "Causal Effects in Non-Experimental Studies: Re-Evaluating the Evaluation of Training Programs," JASA 94 (448), 1053-1062.

  Asad Haris and Gary Chan (2015). ATE: Inference for Average Treatment
  Effects using Covariate Balancing. R package version 0.2.0.
  \url{https://CRAN.R-project.org/package=ATE}

}


\keyword{datasets}
