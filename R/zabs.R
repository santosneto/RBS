#'@name ZARBS
#'
#'@aliases ZARBS
#'@aliases dZARBS
#'@aliases pZARBS
#'@aliases qZARBS
#'@aliases rZARBS
#'@aliases plotZARBS
#'@aliases  meanZARBS
#'
#'@title Zero-Adjusted Reparameterized Birnbaum-Saunders (ZARBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{ZARBS()} defines the ZARBS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}. The zero adjusted Birnbaum-Saunders distribution is similar
#'to the Birnbaum-Saunders distribution but allows zeros as y values. The extra parameter
#'models the probabilities at zero. The functions dZARBS, pZARBS, qZARBS and rZARBS define
#'the density, distribution function, quantile function and random generation for
#'the ZARBS. plotZARBS can be used to plot the distribution. meanZARBS calculates the expected
#'value of the response for a fitted model.
#'
#'@usage ZARBS(mu.link = "log", sigma.link = "log", nu.link = "logit")
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles.
#' @param mu vector of scale parameter values.
#' @param sigma vector of shape parameter values.
#' @param nu vector of mixture parameter values.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from.
#' @param to up to where to plot the distribution.
#' @param obj a fitted ZARBS object.
#' @param ... other graphical parameters for plotting.
#' @param title title of the plot.
#' @param nu.link link function of the parameter nu.
#'
#'
#' @details  The parametrization of the zero adjusted reparameterized Birnbaum-Saunders distribution given in the function ZARBS() is
#'
#'\deqn{ f_{Y}(y;\mu,\delta,p) =\frac{[1-p]\sqrt{\delta+1}}{4\,y^{3/2}\,\sqrt{\pi\mu}}\left[y+\frac{\delta\mu}{\delta+1} \right]\exp\left(-\frac{\delta}{4}\left[\frac{y[\delta+1]}{\delta\mu}+\frac{\delta\mu}{y[\delta+1]}-2\right]\right) I_{(0, \infty)}(y)+  pI_{\{0\}}(y).}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function ZARBS(), mu is the mean and sigma is the precision parameter and nu is the proportion of zeros of the ZARBS distribution.
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2016) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
#'\emph{Applied Stochastic Models in Business and Industry.}, 32(1), 74--89. doi:\email{10.1002/asmb.2124}.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@import gamlss
#'
#'@examples
#'data(oil)
#'fit1 = gamlss::gamlss(Demand~1,sigma.formula=~1, nu.formula=~1,
#'family=ZARBS(mu.link="log",sigma.link = "identity",nu.link = "identity"),
#'method=CG(),data=oil)
#'summary(fit1)
#'
#'@importFrom gamlss.dist checklink
#'@importFrom pracma harmmean
#'@export
#'
ZARBS <- function(mu.link = "log", sigma.link = "log", nu.link = "logit")
{
  mstats <- checklink("mu.link", "ZARBS", substitute(mu.link), c("sqrt", "log", "identity","own"))
  dstats <- checklink("sigma.link", "ZARBS", substitute(sigma.link),c("sqrt", "log", "identity", "own"))
  vstats <- checklink("nu.link", "ZARBS", substitute(nu.link),c("logit", "probit", "cloglog", "own"))
  structure(
                 list(family = c("ZARBS", "Zero Adjusted RBS"),
                 parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
                 nopar = 3,
                 type = "Mixed",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 nu.link = as.character(substitute(nu.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 nu.linkfun = vstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 nu.linkinv = vstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 nu.dr = vstats$mu.eta,
                 dldm = function(y,mu,sigma)
                 {
                   mustart = 1/(2*mu)
                   ystart =  ((sigma+1)*y)/(4*mu*mu) - (sigma^2)/(4*(sigma+1)*y) + (sigma/(sigma+1))*(1/( y + ((mu*sigma)/(sigma+1)) ))
                   dldm <- ifelse((y == 0) , 0, ystart - mustart)
                   dldm
                 },
                 d2ldm2 = function(y,mu,sigma) {
                   d2ldm2 = ifelse((y == 0) , 0, - sigma/(2*mu*mu) - ((sigma/(sigma+1))^2)*Ims(mu,sigma) )
                   d2ldm2
                 },
                 dldd = function(y,mu,sigma) {   #first derivate log-density respect to sigma
                   sigmastart  = -(sigma+2)/(2*(sigma+1))
                   ystart.sigma   = (mu/((sigma+1)^2))*(1/( y + ((mu*sigma)/(sigma+1)) )) - y/(4*mu) - (mu*sigma*(sigma+2))/(4*y*((sigma+1)^2))
                   dldd  = ifelse((y==0),0,(ystart.sigma - sigmastart))
                   dldd
                 },
                 d2ldd2 = function(y,mu,sigma) {
                   lss =  ((sigma^2) + 3*sigma + 1)/(2*((sigma+1)^2)*(sigma^2))
                   d2ldd2 <- ifelse((y==0),0, -lss - ((mu^2)/((sigma+1)^4))*Ims(mu,sigma) )
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu,sigma) {
                   lms = 1/(2*mu*(sigma+1))
                   d2ldmdd <- ifelse((y == 0), 0, - lms - ((mu*sigma)/((sigma+1)^3))*Ims(mu,sigma) )
                   d2ldmdd
                 },
                 dldv = function(y,nu) ifelse(y == 0, 1/nu, -1/(1 - nu)),
                 d2ldv2 = function(nu) -1/(nu*(1 - nu)),
                 d2ldmdv = function(y) rep(0,length(y)),
                 d2ldddv = function(y) rep(0,length(y)),
                 G.dev.incr  = function(y,mu,sigma,nu,...)
                   -2*dZARBS(y,mu,sigma,nu,log=TRUE),
                 rqres = expression(rqres(pfun="pZARBS", type="Mixed",  mass.p=0,
                                          prob.mp=nu, y=y, mu=mu, sigma=sigma, nu=nu)),
                 mu.initial =  expression(mu <- (y+mean(y))/2),
                 #mu.initial =  expression({mu <- rep(median(y),length(y)) }),
                 #sigma.initial =  expression({sigma <- rep( sd(y),length(y)) }),
                 sigma.initial =  expression({sigma <- rep( 1/(sqrt(mean(y[y>0])/harmmean(y[y>0]))-1),length(y)) }),
                 #sigma.initial =  expression(sigma <- rep(1,length(y))),
                 nu.initial =  expression(   nu <- rep(mean(1*(y==0)), length(y))),
                 mu.valid = function(mu) TRUE ,
                 sigma.valid = function(sigma)  all(sigma > 0),
                 nu.valid = function(nu) all(nu > 0) && all(nu < 1),
                 y.valid = function(y)  all(y>=0),
                 mean =  function(mu, sigma, nu) (1 - nu) * mu,
                 variance = function(mu, sigma, nu) (1 - nu)*(mu*mu)*((2*sigma+5)/((sigma+1)^2))
                 ),
                 class = c("gamlss.family","family"))

}

#'@rdname ZARBS
#'
#'@export

dZARBS<-function(x, mu=1, sigma=1, nu=.1, log=FALSE)
{
  if (any(mu < 0) || any(is.na(mu)))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0) || any(is.na(sigma)))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(x < 0))  stop(paste("x must be positive", "\n", ""))
  log.lik <- ifelse(x==0, log(nu), log(1-nu) +  0.5*(sigma - log(mu) + log(sigma+1) - log(16*pi)) - (3/2)*log(x) - ((sigma+1)/(4*mu))*x - ((mu*sigma*sigma)/(4*(sigma+1)))*(1/x)  + log(x + ((mu*sigma)/(sigma+1))))
  if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
  fy
}

#'@rdname ZARBS
#'
#'@export

pZARBS <- function(q, mu=1, sigma=1, nu=0.1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
  a <- sqrt(2/sigma)
  b <- (mu*sigma)/(sigma+1)
  cdf1 <- pnorm((1/a)*(sqrt(q/b) - sqrt(b/q)))
  cdf <- cdf1

  ## the problem with this approximation is that it is not working with
  ## small sigmas and produce NA's. So here it is a solution
  if (any(is.na(cdf)))
  {
    index <- seq(along=q)[is.na(cdf)]
    for (i in index)
    {
      cdf[i] <- integrate(function(x)
        dcbs(x, alpha = a[i], beta = b[i], log=FALSE), 0.001, q[i] )$value
    }
  }
  cdf <- ifelse((q==0), nu, nu+(1-nu)*cdf)
  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf
}


#'@rdname ZARBS
#'
#'@export

qZARBS <- function (p, mu = 0.5, sigma = 1, nu = 0.1, lower.tail = TRUE,
                   log.p = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive ", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)|any(nu >= 1))
    stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  a <- sqrt(2/sigma)
  b <- (mu*sigma)/(sigma+1)

  suppressWarnings(q <- ifelse((nu >= p),0, qcbs((p - nu)/(1-nu),alpha = a, beta = b, lower.tail = TRUE, log.p = FALSE)))
  q
}


#'@rdname ZARBS
#'@importFrom stats runif
#'@export

rZARBS <- function (n, mu = 0.5, sigma = 1, nu = 0.1)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  #In this parametrization  sigma = delta
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = p
    stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qZARBS(p, mu = mu, sigma = sigma, nu = nu)
  r
}

#'@rdname ZARBS
#'@importFrom graphics plot
#'@export

plotZARBS = function (mu = .5, sigma = 1, nu = 0.1, from = 0, to = 0.999, n = 101, title="title",
                     ...)
{
  y <- seq(from = 0.001, to = to, length.out = n)
  pdf <- dZARBS(y, mu = mu, sigma = sigma, nu = nu)
  pr0 <- c(dZARBS(0, mu = mu, sigma = sigma, nu = nu))
  po <- c(0)
  plot(pdf ~ y, main=title, ylim = c(0, max(pdf,pr0)), type = "l",lwd=3)
  points(po, pr0, type = "h",lwd=3)
  points(po, pr0, type = "p", col = "red",lwd=3)
}


#'@rdname ZARBS
#'
#'@importFrom stats fitted
#'@export

meanZARBS <- function(obj)
{
  if (obj$family[1] != "ZARBS")
    stop("the object do not have a ZABS distribution")
  meanofY <- (1 - fitted(obj, "nu")) * fitted(obj, "mu")
  meanofY
}


