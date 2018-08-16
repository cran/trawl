#'Fits a Poisson distribution as marginal law
#'@name fit_marginalPoisson
#'@param x vector of equidistant time series data
#'@param LM Lebesgue measure of the estimated trawl
#'@param plotdiag binary variable specifying whether or not diagnostic plots
#'  should be provided
#'@return v: the rate parameter in the Poisson marginal distribution
#'@details The moment estimator for the Poisson rate parameter is given by
#'  \deqn{\hat v = \mbox{E}(X)/\widehat{ \mbox{LM}}.}
#'@examples
#'\donttest{
#'#Simulate a univariate trawl process and fit the exponential trawl function
#'#and the marginal Poisson law
#'set.seed(1)
#'t <- 1000
#'Delta <- 1
#'v <- 250
#'lambda <- 0.25
#'#Simulate a univariate trawl process with exponential trawl function and
#'#Poisson marginal law
#'trawl <- sim_UnivariateTrawl(t,Delta,burnin=50,marginal =c("Poi"),trawl
#'="Exp",v=v, lambda1=lambda)
#'#Fit the exponential trawl function to the simulated data
#'fittrawlfct <- fit_Exptrawl(trawl,Delta, plotacf=TRUE,lags=500)
#'#Fit the Poisson marginal law
#'fitmarginallaw <- fit_marginalPoisson(trawl, fittrawlfct$LM, plotdiag=TRUE)
#'#Print the results
#'print(paste("lambda: estimated:", fittrawlfct$lambda, ", theoretical:",
#'lambda))
#'print(paste("v: estimated:", fitmarginallaw$v, ", theoretical:", v))
#'}

#'@export
fit_marginalPoisson <- function(x,LM, plotdiag=FALSE){
   v <- base::mean(x)/LM
   n <- base::NROW(x)

  if(plotdiag){

    graphics::par(mfrow=base::c(1,2))
    MASS::truehist(x,xlab="")
    graphics::lines(stats::dpois(1:n, v*LM))

    stats::qqplot(x,stats::rpois(1:n, v*LM), xlab="", ylab="Poisson")
    graphics::abline(0,1)
    graphics::par(mfrow=base::c(1,1))
  }

  return (base::list("v"=v))
}


#'Fist a negative binomial distribution as marginal law
#'@name fit_marginalNB
#'@param x vector of equidistant time series data
#'@param LM Lebesgue measure of the estimated trawl
#'@param plotdiag binary variable specifying whether or not diagnostic plots
#'  should be provided
#'@return m: parameter in the negative binomial marginal distribution
#'@return theta: parameter in the negative binomial marginal distribution
#'@return a: Here \eqn{a=\theta/(1-\theta)}. This is given for an alternative
#'  parametrisation of the  negative binomial marginal distribution
#'@details The moment estimator for the parameters of the negative binomial
#'  distribution are given by \deqn{\hat \theta = 1-\mbox{E}(X)/\mbox{Var}(X),}
#'  and \deqn{\hat m = \mbox{E}(X)(1-\hat \theta)/(\widehat{ \mbox{LM}} \hat
#'  \theta).}
#' @examples
#'\donttest{
#'#Simulate a univariate trawl process and fit the exponential trawl function
#'#and the marginal negative binomial law
#'set.seed(1)
#'t <- 1000
#'Delta <- 1
#'m<-200
#'theta<-0.5
#'lambda <- 0.25
#'#Simulate a univariate trawl process with exponential trawl function and
#'#negative binomial marginal law
#'trawl <- sim_UnivariateTrawl(t,Delta,burnin=50,marginal =c("NegBin"),trawl
#'="Exp",m=m, theta=theta, lambda1=lambda)
#'#Fit the exponential trawl function to the simulated data
#'fittrawlfct <- fit_Exptrawl(trawl,Delta, plotacf=TRUE,lags=500)
#'#Fit the Poisson marginal law
#'fitmarginallaw <- fit_marginalNB(trawl, fittrawlfct$LM, plotdiag=TRUE)
#'#Print the results
#'print(paste("lambda: estimated:", fittrawlfct$lambda, ", theoretical:",
#'lambda))
#'print(paste("m: estimated:", fitmarginallaw$m, ", theoretical:", m))
#'print(paste("theta: estimated:", fitmarginallaw$theta, ", theoretical:",
#'theta))
#'}

#'@export
fit_marginalNB <- function(x,LM, plotdiag=FALSE){
  n <- base::NROW(x)

  theta <- 1 -base::mean(x)/stats::var(x)
  m     <- base::mean(x)/LM*(1-theta)/theta
  a     <- theta/(1-theta)

  if(plotdiag){

    graphics::par(mfrow=base::c(1,2))
    MASS::truehist(x,xlab="")
    graphics::lines(stats::dnbinom(1:n, size=m*LM,  mu=m*LM*theta/(1-theta), log = FALSE))

    stats::qqplot(x,stats::rnbinom(1:n, size=m*LM,  mu=m*LM*theta/(1-theta)), xlab="", ylab="Negative Binomial")
    graphics::abline(0,1)
    graphics::par(mfrow=base::c(1,1))
  }

  return (base::list("m"=m, "theta"=theta, "a"=a))
}


