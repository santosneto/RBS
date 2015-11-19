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



#'Brown and Miller's biaxial fatigue data
#'
#'Several data sets related to the BS distribution are available in the RBS package, which have been taken from the literature on this topic.
#'
#'@docType data
#'@keywords datasets
#'@name biaxial
#'@usage data(biaxial)
#'@format A vector containing 46 observations.
#'@details Rieck and Nedelman (1991) reported a data set from Brown and Miller (1978). These data represent results of fatigue tests on 1 \eqn{\%} Cr-Mo-V steel. Cylindrical specimens were subjected to combined torsional and axial loads over constant-amplitude cycles until failure.
#'
#'@references
#'
#'Rieck, J. R. and Nedelman, J. (1991). A Log-Linear Model for the Birnbaum-Saunders Distribution. Technometrics. 33, 51-60.
#'
#'Brown, M. W. and Miller, K. J. (1978). Biaxial Fatigue Data. Report CEMR1/78. University of Sheffield, Dept. of Mechanical Engineering.
#'
#'@examples
#'## Load data sets
#'data(biaxial)
#'## Histogram for biaxial
#'hist(biaxial)
#'
NULL


#'Fatigue lifetime data
#'
#'Several data sets related to the BS distribution are available in the RBS package, which have been taken from the literature on this topic.
#'
#'@docType data
#'@keywords datasets
#'@name bearings
#'@usage data(bearings)
#'@format A vector containing 10 observations.
#'@details McCool (1974) gives data related to the fatigue life in hours of 10 bearings of a certain type.
#'
#'@references
#'
#'McCool, J. I. (1974). Inferential techniques for Weibull populations. Aerospace Research Laboratories Report ARL TR74-0180, Wright-Patterson Air Force Base, Dayton, OH.
#'
#'@examples
#'## Load data sets
#'data(bearings)
#'## Histogram for bearings
#'hist(bearings)
#'
NULL


#'Lifetimes of aluminum specimens exposed to a maximum stress of 21,000 psi.
#'
#'Several data sets related to the BS distribution are available in the RBS package, which have been taken from the literature on this topic.
#'
#'@docType data
#'@keywords datasets
#'@name psi21
#'@usage data(psi21)
#'@format A vector containing 101 observations.
#'@details psi21, psi26 and psi31 were taken from Birnbaum and Saunders (1969),
#'who reported fatigue life data correspond to the cycles (\eqn{\times 10^{-3}}) of aluminum
#'specimens of type 6061-T6. These specimens were cut in a parallel angle to the direction of
#'rotation and oscillating at 18 cycles per seconds. They were exposed to a pressure with maximum
#'stress of 21,000, 26,000, and 31,000 psi (pounds per square inch), for \eqn{n = 101, 102,}
#'and \eqn{101} specimens, respectively.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). Estimation for a family of life distributions with applications to fatigue. J. Appl. Probab. 6(2): 328-347.
#'
#'@examples
#'## Load data sets
#'data(psi21)
#'## Histogram for psi21
#'hist(psi21)
#'
NULL

#'Lifetimes of aluminum specimens exposed to a maximum stress of 26,000 psi.
#'
#'Several data sets related to the BS distribution are available in the RBS package, which have been taken from the literature on this topic.
#'
#'@docType data
#'@keywords datasets
#'@name psi26
#'@usage data(psi26)
#'@format A vector containing 102 observations.
#'@details psi21, psi26 and psi31 were taken from Birnbaum and Saunders (1969),
#'who reported fatigue life data correspond to the cycles (\eqn{\times 10^{-3}}) of aluminum
#'specimens of type 6061-T6. These specimens were cut in a parallel angle to the direction of
#'rotation and oscillating at 18 cycles per seconds. They were exposed to a pressure with maximum
#'stress of 21,000, 26,000, and 31,000 psi (pounds per square inch), for \eqn{n = 101, 102,}
#'and \eqn{101} specimens, respectively.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). Estimation for a family of life distributions with applications to fatigue. J. Appl. Probab. 6(2): 328-347.
#'
#'@examples
#'## Load data sets
#'data(psi26)
#'## Histogram for psi26
#'hist(psi26)
#'
NULL


#'Lifetimes of aluminum specimens exposed to a maximum stress of 31,000 psi.
#'
#'Several data sets related to the BS distribution are available in the RBS package, which have been taken from the literature on this topic.
#'
#'@docType data
#'@keywords datasets
#'@name psi31
#'@usage data(psi31)
#'@format A vector containing 101 observations.
#'@details psi21, psi26 and psi31 were taken from Birnbaum and Saunders (1969),
#'who reported fatigue life data correspond to the cycles (\eqn{\times 10^{-3}}) of aluminum
#'specimens of type 6061-T6. These specimens were cut in a parallel angle to the direction of
#'rotation and oscillating at 18 cycles per seconds. They were exposed to a pressure with maximum
#'stress of 21,000, 26,000, and 31,000 psi (pounds per square inch), for \eqn{n = 101, 102,}
#'and \eqn{101} specimens, respectively.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). Estimation for a family of life distributions with applications to fatigue. J. Appl. Probab. 6(2): 328-347.
#'
#'@examples
#'## Load data sets
#'data(psi31)
#'## Histogram for psi31
#'hist(psi31)
#'
NULL

