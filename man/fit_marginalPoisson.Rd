% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitMarginalDistributions.R
\name{fit_marginalPoisson}
\alias{fit_marginalPoisson}
\title{Fits a Poisson distribution as marginal law}
\usage{
fit_marginalPoisson(x, LM, plotdiag = FALSE)
}
\arguments{
\item{x}{vector of equidistant time series data}

\item{LM}{Lebesgue measure of the estimated trawl}

\item{plotdiag}{binary variable specifying whether or not diagnostic plots
should be provided}
}
\value{
v: the rate parameter in the Poisson marginal distribution
}
\description{
Fits a Poisson distribution as marginal law
}
\details{
The moment estimator for the Poisson rate parameter is given by
 \deqn{\hat v = \mbox{E}(X)/\widehat{ \mbox{LM}}.}
}
