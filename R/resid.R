#'Residuals
#'
#'@description \code{residuals} is a function which extracts model residuals from objects returned by modeling functions.
#'@usage residuals(model,residual = "deviance")
#'
#' @param model an object for which the extraction of model residuals is meaningful.
#' @param residual type of residual to be used.
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
#'@examples
#'##
#'library(alr3)
#'data(landrent)
#'attach(landrent)
#'resp <- I(Y/X1)
#'y1 <-  split(resp, X4)$"1"
#'x21 <-  split(X2, X4)$"1"
#'
#'##Fixed Precision
#'fit0 <- gamlss(y1 ~ x21, family=RBS(mu.link="identity"),method=CG()  )
#'plot(fitted(fit0),residuals(fit0),xlab="fitted values",ylab="Deviance")
#'##Varying Precision
#'fit1 <- gamlss(y1 ~ x21,sigma.formula = y1 ~x21, family=RBS(mu.link="identity",sigma.link="sqrt"),method=CG()  )
#'plot(fitted(fit1),residuals(fit1),xlab="fitted values",ylab="Deviance")
#'
#'
#'data(Snacks,package="ssym")
#'head(Snacks)
#'attach(Snacks)
#'type <- factor(type,labels = c("A","B","C","D","E"))
#'w1 <- week
#'w2 <- I(week^2)
#'fit. <- gamlss(texture~type+w1+w2,~type, family=RBS(mu.link="log",sigma.link = "log"),method=CG())
#'plotreg(fit., custom.model.names = "",custom.note = "CI")
#'summary(fit.)
#'plot(fit.$mu.fv,residuals(fit.,residual="deviance"),ylab="rd",xlab="m",pch=19,lwd=2,ylim=c(-4,4))
#'abline(h=2,lwd=2,lty=2)
#'abline(h=-2,lwd=2,lty=2)
#'text(fit.$mu.fv[c(91)],residuals(fit.,residual="deviance")[c(91)]-.2,c(91))
#'
#'
#'
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
    intemod <-  mapply(esp,mu.hat,sigma.hat)

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
