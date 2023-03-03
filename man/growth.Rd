\name{Growth}
\docType{data}
\alias{Growth}
\title{Growth Data}
\description{
Panel of Macroeconomic data for 125 countries from 1960 to 1985
constructed by Summers and Heston (1991))
}
\usage{data(Growth)}
\format{
A data frame containing 9 vectors. 
\describe{
  \item{Country_ID}{Country identification number}
  \item{COM}{1 if the country is in a communist regime, 0 otherwise} 
  \item{OPEC}{1 if the country is part of the OPEC, 0 otherwise} 
  \item{Year}{Year} 
  \item{GDP}{Per capita GDP (in thousands) in 1985 U.S. dollars.} 
  \item{LagGDP}{GDP of the previous period}
  \item{SavRate}{Saving rate measured as the ratio of real investment to
  real GDP}
  \item{LagSavRate}{SavRate of the previous period}
  \item{Country}{Country names}
  \item{Pop}{Population in thousands}
  \item{LagPop}{Population of the previous period}
}
}

\source{\url{http://fhayashi.fc2web.com/datasets.htm}}
\keyword{datasets}
