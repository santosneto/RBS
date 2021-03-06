#'Residuals
#'
#'@description \code{residuals} is a function which extracts model residuals from objects returned by modeling functions.
#'@usage residuals(model,residual = "deviance")
#'
#' @param model an object for which the extraction of model residuals is meaningful.
#' @param residual type of residual to be used ("pearson","score", "deviance" or "quantile"). Default is the deviance.
#'
#'@return Residuals extracted from the object object.
#'
#'@author
#'Manoel Santos-Neto \url{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \url{cysneiros@de.ufpe.br}, Victor Leiva \url{victorleivasanchez@gmail.com} and Michelli Barros \url{michelli.karinne@gmail.com}
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'@export
residuals <- function(model,residual = "deviance")
{

  rP <- function(fit)
  {

    mu.hat <- fit$mu.fv
    sigma.hat <- fit$sigma.fv
    y <- model$y

    dif <- (y-mu.hat)
    phi.hat <- sqrt((2*((sigma.hat+1)^2))/((2*sigma.hat) +5))
    raiz <- sqrt(2*mu.hat*mu.hat)
    res <- (dif*phi.hat)/raiz
    return(res)
  }

  rS <- function(fit)
  {

    mu.hat <- fit$mu.fv
    sigma.hat <- fit$sigma.fv
    y <- model$y
    intemod <-  Ims(mu.hat,sigma.hat)

    z <- -1/(2*mu.hat) - (sigma.hat^2)/(4*(sigma.hat+1)*y) + ((sigma.hat+1)*y)/(4*mu.hat*mu.hat) + sigma.hat/((sigma.hat*y) + y + (sigma.hat*mu.hat))
    v <- sigma.hat/(2*mu.hat*mu.hat) + ((sigma.hat*sigma.hat)/((sigma.hat+1)*(sigma.hat+1)))*intemod
    res <- z/sqrt(v)
    return(res)
  }


  rD <- function(fit)
  {
    mu.hat <- fit$mu.fv
    sigma.hat <- fit$sigma.fv
    y <- model$y
    dif <- y - mu.hat
    sinal <- sign(dif)
    const <- sqrt(2)
    term0 <- log(2)
    term1 <- -(sigma.hat/2)
    term21 <- (y*(sigma.hat+1))/(4*mu.hat)
    term22 <- ((sigma.hat^2)*mu.hat)/(4*y*(sigma.hat+1))
    term23<-(1/2)*log((y*sigma.hat*mu.hat*(sigma.hat+1))/((sigma.hat*y+y+sigma.hat*mu.hat)^2))
    res <- sinal*const*sqrt(term0 + term1 + term21 + term22+ term23)
    return(res)
  }


  rQ<- function(fit)
  {
    mu.hat <- fit$mu.fv
    sigma.hat <- fit$sigma.fv
    y <- model$y
    F <- pRBS(y,mu=mu.hat,sigma=sigma.hat)
    res <- qnorm(F)
    return(res)
  }

  output <- switch(residual,
                   "pearson" = rP(model),
                   "score" = rS(model),
                   "deviance" = rD(model),
                   "quantile" = rQ(model))
  return(output)
}
