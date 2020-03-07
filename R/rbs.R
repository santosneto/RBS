#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'
dcbs <- function(x, alpha = 1, beta = 1, log = FALSE){
  if(!is.numeric(x)||!is.numeric(alpha)||!is.numeric(beta)){
    stop("non-numeric argument to mathematical function")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  x   <- x
  c   <-(1 / sqrt(2 * pi))
  u   <- (alpha^(-2)) * ((x / beta) + (beta / x) - 2)
  e   <- exp((-1 / 2 ) * u)
  du  <- (x^(-3 / 2) * (x + beta)) /
    (2 * alpha * sqrt (beta))
  pdf <- c * e * du
  if (log==TRUE){pdf <-log(pdf)}
  return(pdf)
}

#--------------------------------------------------------------------------------------------
#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'

pcbs <- function(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if(!is.numeric(q)||!is.numeric(alpha)||!is.numeric(beta)){
    stop("non-numeric argument to mathematical function")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  x   <- q
  s   <- (x / beta)
  a   <- ((1 / alpha) * (s^(1 / 2) - s^(-1 / 2)))
  cdf <- pnorm(a, 0, 1)
  if (lower.tail == FALSE){cdf <-(1 - cdf)}
  if (log.p == TRUE){cdf <-log(cdf)}
  return(cdf)
}

#--------------------------------------------------------------------------------------------
#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'
qcbs <- function(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  if (log.p == TRUE){p  <- log(p)}
  if (lower.tail == FALSE){p <- (1 - p)}
  q   <- beta * (((alpha * qnorm(p, 0, 1) / 2) +
                    sqrt(((alpha * qnorm(p, 0, 1) / 2)^2) +
                           1)))^2
  return(q)
}

#--------------------------------------------------------------------------------------------
#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'
rcbs <- function(n, alpha = 1, beta = 1)
{
  if (!is.numeric(n)||!is.numeric(alpha)||!is.numeric(beta))
  {stop("non-numeric argument to mathematical function")}
  if (n == 0){stop("value of n must be greater or equal then 0")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  z   <- rnorm(n, 0, 1)
  t   <- beta * (1 + (((alpha^2) * (z^2)) / 2) +
                   (alpha * z *sqrt((((alpha^2) * (z^2))/ 4)+ 1)))
  return(t)
}


#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
RBS <- function (mu.link = "log" , sigma.link="log")
{
  mstats = checklink("mu.link", "BirnbaumSaunders", substitute(mu.link),
                     c("sqrt","log","identity"))
  dstats = checklink("sigma.link", "BirnbaumSaunders", substitute(sigma.link)
                     ,c("sqrt", "log", "identity"))
  structure(list(family = c("RBS","BirnbaumSaunders"),
                 parameters = list(mu=TRUE,sigma=TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 #the first derivative of the likelihood with respect to the location parameter mu
                 dldm = function(y,mu,sigma) #first derivate of log-density respect to mu
                 {
                   mustart = 1/(2*mu)
                   ystart =  ((sigma+1)*y)/(4*mu*mu) - (sigma^2)/(4*(sigma+1)*y) + sigma/((sigma*y) + y + (sigma*mu))

                   dldm = ystart-mustart
                   dldm
                 },
                 #the expected second derivative of the likelihood with respect to the location parameter mu
                 d2ldm2 = function(mu,sigma) {        #expected of second derivate of log-density respect to mu
                   Ims = mapply(esp,mu,sigma)
                   d2ldm2 =  - sigma/(2*mu*mu) - ((sigma/(sigma+1))^2)*Ims
                   d2ldm2
                 },

                 #the first derivative of the likelihood with respect to the scale parameter sigma
                 dldd = function(y,mu,sigma) {      #first derivate log-density respect to sigma
                   sigmastart  = -(sigma)/(2*(sigma+1))
                   y2start   = (y+mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (mu*sigma*(sigma+2))/(4*y*((sigma+1)^2))
                   dldd  = y2start-sigmastart
                   dldd
                 },
                 #the expected second derivative of the likelihood with respect to the scale parameter sigma ok
                 d2ldd2 = function(mu,sigma) {      #expected of second derivate log-density respect to sigma
                   Ims = mapply(esp,mu,sigma)
                   lss =  ((sigma^2) + (3*sigma) + 1)/(2*sigma*sigma*((sigma+1)^2))
                   d2ldd2 = -lss - ((mu^2)/((sigma+1)^4))*Ims
                   d2ldd2
                 },
                 #the expected cross derivative of the likelihood with respect to both the location mu and scale parameter sigma
                 d2ldmdd = function(mu,sigma) {   #expected of partial derivate of log-density respect to mu and sigma
                   Ims = mapply(esp,mu,sigma)
                   lms = 1/(2*mu*(sigma+1))
                   d2ldmdd = - lms - ((mu*sigma)/((sigma+1)^3))*Ims
                   d2ldmdd
                 },


                 G.dev.incr = function(y,mu,sigma,...) -2*dRBS(y,mu,sigma,log=TRUE),

                 rqres = expression(resrbs(y=y,mu=mu,sigma=sigma)),

                 mu.initial = expression({mu = mean(y)}),
                 sigma.initial = expression({sigma = rep(sigmatil(y),length(y))}),
                 mu.valid = function(mu) all(mu>0) ,
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0),
                 mean = function(mu,sigma) mu,
                 variance = function(mu, sigma) (mu*mu)*((2*sigma+5)/((sigma+1)^2)) ),
            class = c("gamlss.family","family"))
}


#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
sigmatil=function(y)
{
  s = mean(y)
  r = 1/mean(1/y)
  alphatil = (2*( (s/r)^(1/2)  - 1))^(1/2)
  dest = 2/(alphatil^2 )
  return(dest)
}

#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
esp = function(mu=1,sigma=1)
{

  integral=function(aest)
  {
    fu=function(u)
    {
      w1 = (1 / ((1 +u^2)*(u^2)))
      w2 = (exp((-1 /(2*aest^2) )*((u - 1/u)^2)))
      (w1*w2)
    }
    return(integrate(fu,0,Inf)$value)
  }

  const = function(alpha,beta)
  {
    const = 1/(alpha*beta*beta*sqrt(2*pi))
    return(const)
  }

  alpha = sqrt(2/sigma)
  beta = (mu*sigma)/(sigma+1)
  e = const(alpha,beta)*integral(alpha)
  return(e)
}

#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
resrbs=function(y,mu,sigma)
{
  Ims = mapply(esp,mu,sigma)
  z = -1/(2*mu) - (sigma^2)/(4*(sigma+1)*y) + ((sigma+1)*y)/(4*mu*mu) + sigma/((sigma*y) + y + (sigma*mu))
  v = sigma/(2*mu*mu) + ((sigma*sigma)/((sigma+1)*(sigma+1)))*Ims
  res = z/sqrt(v)
  return(res)
}


#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
dRBS<-function(x, mu=1, sigma=1, log=FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0))  stop(paste("x must be positive", "\n", ""))
  log.lik =  0.5*(sigma - log(mu) + log(sigma+1) - log(16*pi)) - (3/2)*log(x) - ((sigma+1)/(4*mu))*x - ((mu*sigma*sigma)/(4*(sigma+1)))*(1/x)  + log(x + ((mu*sigma)/(sigma+1)))
  if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
  fy
}

#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
pRBS <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
{       if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
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

  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf
}

#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
qRBS = function (p, mu = 0.5, sigma = 1, lower.tail = TRUE,
                 log.p = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive ", "\n", ""))
  if (any(sigma < 0))  #In this parametrization  sigma = phi
    stop(paste("sigma must be positive", "\n", ""))


  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p <= 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)

  suppressWarnings(q <- qcbs(p ,alpha = a, beta = b, lower.tail = TRUE, log.p = FALSE))
  q
}


#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'

rRBS = function(n, mu=1, sigma=1)
{
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  r = rcbs(n, alpha = a, beta = b)
  r
}


#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'
plotRBS = function(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, title="title", ...)
{
  y = seq(from = 0.001, to = to, length.out = n)
  pdf = dRBS(y, mu = mu, sigma = sigma)
  plot(pdf ~ y, main=title, ylim = c(0, max(pdf)), type = "l",lwd=3)

}

#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'

meanRBS = function (obj)
{
  if (obj$family[1] != "RBS")
    stop("the object do not have a RBS distribution")
  meanofY = fitted(obj, "mu")
  meanofY
}


#--------------------------------------------------------------------------------------------
#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\sigma)=\frac{\exp\left(\sigma/2\right)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\sigma \mu}{\sigma+1}\right] \exp\left(-\frac{\sigma}{4}
#'                                                     \left[\frac{y\{\sigma+1\}}{\sigma\mu}+\frac{\sigma\mu}{y\{\sigma+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2012) On new parameterizations of the Birnbaum-Saunders distribution. \emph{PAK J STAT}, v. 28, p. 1-26, 2012.
#'
#'Santos-Neto, M., Cysneiros, F.J.A, Leiva, V., Barros, M. (2014) On a reparameterized Birnbaum-Saunders distribution and its moments, estimation and application. \emph{Revstat Statistical Journal}, v. 12, p. 247-272, 2014.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'

est.rbs <- function(x,xi=0.95)
{

n <- length(x)

ic.bs <- function(modelo=NULL,estimates = NULL,n=NULL,method = "lr",level=0.95)
{

  zlevel <- abs(qnorm((1-level)/2))
  est <- estimates
  if(is.null(est) == TRUE && is.null(modelo) == FALSE)
  {


    if(method == "lr")
    {
      est.mu <- modelo$mu.coefficients
      est.delta <- modelo$sigma.coefficients
      var. <- diag(vcov(modelo))
      se.mu <- sqrt(var.[1])
      se.delta <- sqrt(var.[2])
      li.mu <- est.mu - se.mu*zlevel
      li.delta <- est.delta- se.delta*zlevel
      ls.mu <- est.mu + se.mu*zlevel
      ls.delta <- est.delta + se.delta*zlevel

    }
    else
    {
      est.mu <- coef(modelo)[1]
      est.delta <- coef(modelo)[2]
      var. <- diag(vcov(modelo))
      se.mu <- sqrt(var.[1])
      se.delta <- sqrt(var.[2])
      li.mu <- est.mu - se.mu*zlevel
      li.delta <- est.delta- se.delta*zlevel
      ls.mu <- est.mu + se.mu*zlevel
      ls.delta <- est.delta + se.delta*zlevel
    }

  }
  else{
    if(method == "mm" )
    {

      est.mu <- est[1]
      est.delta <- est[2]

      var.mu  <- (1/n)*(((est.mu^2)*(2*est.delta + 5))/((est.delta+1)^2))
      var.delta <- (1/n)*( (2*(est.delta^4) + 28*(est.delta^3) + 122*(est.delta^2) +  126*est.delta + 57)/((est.delta+4)^2))

      se.mu <- sqrt(var.mu)
      se.delta <- sqrt(var.delta)
      li.mu <- est.mu - se.mu*zlevel
      li.delta <- est.delta - se.delta*zlevel
      ls.mu <- est.mu + se.mu*zlevel
      ls.delta <- est.delta + se.delta*zlevel

    }

    else{

      est.mu <- est[1]
      est.delta <- est[2]

      var.mu  <- (1/n)*(((est.mu^2)*(2*est.delta + 5))/((est.delta+1)^2))
      var.delta <- (1/n)*(2*(est.delta^2))

      se.mu <- sqrt(var.mu)
      se.delta <- sqrt(var.delta)
      li.mu <- est.mu - se.mu*zlevel
      li.delta <- est.delta - se.delta*zlevel
      ls.mu <- est.mu + se.mu*zlevel
      ls.delta <- est.delta + se.delta*zlevel

    }

  }

  ic.mu <- round(c(li.mu,ls.mu),2)
  ic.delta <- round(c(li.delta,ls.delta),2)

  ics <- list(ic.mu = ic.mu,ic.delta = ic.delta)

  return(ics)
}


g1 <- function(vP,x)
{
  mu <- vP[1];delta=vP[2]
  m1 <- (mu - x)
  m2 <- (((mu^2)*((2*delta)+5))/((delta+1)^2) - (x - mu)^2)
  m3 <- (((delta+1)^2)/(mu*(delta^2)) - 1/x)
  f <- cbind(m1,m2,m3)
  return(f)
}

Dg1 <- function(vP,x)
{
  mu <- vP[1];delta=vP[2]
  g11 <- 1 ;g12 = 0
  g21 <- (2*mu*(2*delta+5))/((delta+1)^2) - 2*mu + 2*mean(x)
  g22 <-  -(2*(mu^2)*(delta+4))/((delta+1)^3)
  g31 <- -((delta+1)^2)/((mu*delta)^2)
  g32 <- -(2*(delta+1))/(mu*(delta^3))

  G <- matrix(c(g11,g21,g31,g12,g22,g32),nrow=3,ncol=2)
  return(G)
}

mum <- mean(x)
xbarh <- 1/mean(1/x)
deltamm <- 1/(sqrt(mum/xbarh) - 1)
s2 <- ((n-1)/n)*var(x)
deltam <- (mum^2 - s2 + sqrt((mum^4) + 3*(mum^2)*s2))/s2
con <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
estmv <- gamlss(x~1,family=RBS(mu.link="identity",sigma.link="identity"),control=con,method=CG())
res <- gmm(g1,x,c(mu = mum, delta = deltamm),gradv = Dg1)
teste <- specTest(res,theta0=c(mu,delta))
pvalue <- teste$test[2]

vpm <- round(c(mum,deltam),2)
vpmm <- round(c(mum,deltamm),2)
vpmg <- round(res$coefficients,2)
vpmv <- round(c(estmv$mu.coefficients,estmv$sigma.coefficients),2)

ic.mv  <- ic.bs(modelo=estmv,method="lr",level=xi)
ic.mm  <- ic.bs(estimates=c(mum,deltam),n = n, method="mm",level=xi)
ic.mmm <- ic.bs(estimates=c(mum,deltamm),n=n, method="mmm",level=xi)
ic.mmg <- ic.bs(modelo=res,method="gmm",level=xi)

r1<-rbind(c(ic.mm$ic.mu[1],ic.mm$ic.mu[2]),c(ic.mm$ic.delta[1],ic.mm$ic.delta[2]))
r2<-rbind(c(ic.mmm$ic.mu[1],ic.mmm$ic.mu[2]),c(ic.mmm$ic.delta[1],ic.mmm$ic.delta[2]))
r3<-rbind(c(ic.mmg$ic.mu[1],ic.mmg$ic.mu[2]),c(ic.mmg$ic.delta[1],ic.mmg$ic.delta[2]))
r4<-rbind(c(ic.mv$ic.mu[1],ic.mv$ic.mu[2]),c(ic.mv$ic.delta[1],ic.mv$ic.delta[2]))

result <- cbind(vpm,r1,vpmm,r2,vpmg,r3,vpmv,r4)
colnames(result) <- c("MO","Lower","Upper","MM","Lower","Upper","GMM","Lower","Upper","MLE","Lower","Upper")
rownames(result) <- c("mu","sigma")
return(result)

}





