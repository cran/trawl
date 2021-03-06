% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitCrossTerms.R
\name{fit_trawl_intersection_LM}
\alias{fit_trawl_intersection_LM}
\title{Finds the intersection of two long memory (LM) trawl sets}
\usage{
fit_trawl_intersection_LM(alpha1, H1, alpha2, H2, LM1, LM2, plotdiag = FALSE)
}
\arguments{
\item{alpha1, H1, alpha2, H2}{parameters of the two long memory trawls}

\item{LM1}{Lebesgue measure of the first trawl}

\item{LM2}{Lebesgue measure of the second trawl}

\item{plotdiag}{binary variable specifying whether or not diagnostic plots
should be provided}
}
\value{
the Lebesgue measure of the intersection of the two trawl sets
}
\description{
Finds the intersection of two long memory (LM) trawl sets
}
\details{
Computes \eqn{R_{12}(0)=\mbox{Leb}(A_1 \cap A_2)} based on two trawl
 functions \eqn{g_1} and \eqn{g_2}.
}
