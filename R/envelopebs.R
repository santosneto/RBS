#'@name envelope
#'
#'@aliases envelope
#'@aliases envelope.ZARBS
#'@aliases envelope.ZAGA
#'
#'@title Envelopes
#'
#'@description A normal plot with simulated envelope of the residual is produced.
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param res type of residuals to be extracted. Default is deviance.
#' @param precision If \code{precision = "fixed"} a model with precision fixed is used;
#' else a model with precision variable is used.
#' @param dist The function RBS() defines the RBS distribution.
#' @param xlabel a label for the x axis.
#' @param color the color of the envelope.
#' @param ylabel a label for the y axis.
#' @param font the font used in x and y axis.
#'
#'
#'@return A simulated envelope of the class RBS, ZARBS and ZAGA.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'@importFrom graphics par points polygon
#'@importFrom stats qqnorm
#'@importFrom gamlss gamlss gamlss.control
#'@import ggplot2
#'@export
envelope <- function(model, k = 100, res="deviance", precision = "fixed", dist = RBS(mu.link ="identity", sigma.link = "identity"), color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{
  if(precision != "fixed")
  {
    n <- model$N
    td  <- residuals(model,residual=res)
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    re <- matrix(0,n,k)
    X <- model$mu.x
    Z <- model$sigma.x
    p<- ncol(X)
    q<- ncol(Z)

    for(i in 1:k)
    {

      y1 <- mapply(rRBS,n=1,mu=mu,sigma=sigma)
      nresp <- y1

      if(p==1){form <- nresp ~ 0 + X} else {form <- nresp ~ X[,-1]}
      if(q==1){form1 <- nresp ~ 0 + Z} else {form1 <- nresp ~ Z[,-1]}

      conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form, sigma.formula = form1 ,family=dist, method=CG(),control = conh0)
      rdf <- residuals(model1,residual=res)

      re[,i]=sort(rdf)
    }

    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    xab <- NULL
    emin.e10 <- NULL
    emax.e20 <- NULL
    emin.e11 <- NULL
    emax.e21 <- NULL
    emin.e12 <- NULL
    emax.e22 <- NULL

    for(l in 1:n)
    {
      eo = sort(re[l,])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }

    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb <- apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x
    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=10,family=font))
    } else{
    n<-model$N
    td  <- residuals(model, residual= res)
    sigma <- model$sigma.fv[1]
    mu <- model$mu.fv
    alpha<-sqrt(2/sigma)
    bet.a<- (sigma*mu)/(sigma+1)
    re <- matrix(0,n,k)
    X <- model$mu.x
    p <- ncol(X)
    mulink <- model$mu.link


    for(i in 1:k)
    {

      y1 <- mapply(rRBS,n=1,mu=mu,sigma=sigma)
      nresp <- y1

      if(p==1){form <- nresp ~ 0 + X} else{form <- nresp ~ X[,-1]}

      conh0 <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form ,family=dist,method=CG(),control = conh0)
      rdf <- residuals(model1,residual=res)
      re[,i] <- sort(rdf)
    }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    xab <- NULL
    emin.e10 <- NULL
    emax.e20 <- NULL
    emin.e11 <- NULL
    emax.e21 <- NULL
    emin.e12 <- NULL
    emax.e22 <- NULL
    for(l in 1:n)
    {
      eo = sort(re[l,])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }

    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb <- apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x
    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=10,family=font))
  }
}

#'@rdname envelope
#'
#'@importFrom graphics par points polygon
#'@importFrom stats qqnorm
#'@importFrom gamlss gamlss gamlss.control
#'@import ggplot2
#' @export
envelope.ZARBS <- function (model, k = 100, precision = "fixed", color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{
  if (precision != "fixed") {
    n <- model$N
    td <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    nu <- model$nu.fv
    re <- matrix(0, n, k)
    X <- model$mu.x
    Z <- model$sigma.x
    W <- model$nu.x
    p <- ncol(X)
    q <- ncol(Z)
    s <- ncol(W)
    for (i in 1:k) {
      y1 <- rnorm(n)
      re[, i] = sort(y1)
    }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    xab <- NULL
    emin.e10 <- NULL
    emax.e20 <- NULL
    emin.e11 <- NULL
    emax.e21 <- NULL
    emin.e12 <- NULL
    emax.e22 <- NULL

    for (l in 1:n) {
      eo = sort(re[l, ])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }
    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb = apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x
    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td)) + geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point(size=0.5) + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=10,family=font))
  }
  else {

    n <- model$N
    td <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    nu <- model$nu.fv
    re <- matrix(0, n, k)
    X <- model$mu.x
    W <- model$nu.x
    p <- ncol(X)
    s <- ncol(W)
    mulink <- model$mu.link
    for (i in 1:k) {
      y1 <- rnorm(n)
      re[, i] <- sort(y1)
    }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    xab <- NULL
    emin.e10 <- NULL
    emax.e20 <- NULL
    emin.e11 <- NULL
    emax.e21 <- NULL
    emin.e12 <- NULL
    emax.e22 <- NULL

    for (l in 1:n) {
      eo = sort(re[l, ])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }
    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb <- apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x
    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=10,family=font))
  }
}


#'@rdname envelope
#'
#'@importFrom graphics par points polygon
#'@importFrom stats qqnorm
#'@importFrom gamlss gamlss gamlss.control
#'@importFrom gamlss.dist rZAGA ZAGA
#'@import ggplot2
#' @export

envelope.ZAGA <- function (model, k = 100, precision = "fixed", color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{
  if (precision != "fixed") {
    n <- model$N
    td <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    nu <- model$nu.fv
    re <- matrix(0, n, k)
    X <- model$mu.x
    Z <- model$sigma.x
    W <- model$nu.x
    p <- ncol(X)
    q <- ncol(Z)
    s <- ncol(W)
    for (i in 1:k) {
      y1 <- rnorm(n)
      re[, i] <- sort(y1)
      }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    xab <- NULL
    emin.e10 <- NULL
    emax.e20 <- NULL
    emin.e11 <- NULL
    emax.e21 <- NULL
    emin.e12 <- NULL
    emax.e22 <- NULL

    for (l in 1:n) {
      eo = sort(re[l, ])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }
    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb = apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x
    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td)) + geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=10,family=font))
  }
  else {

    n <- model$N
    td <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    nu <- model$nu.fv
    re <- matrix(0, n, k)
    X <- model$mu.x
    W <- model$nu.x
    p <- ncol(X)
    s <- ncol(W)
    mulink <- model$mu.link
    for (i in 1:k) {
      y1 <- rnorm(n)
      re[, i] <- sort(y1)
    }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    xab <- NULL
    emin.e10 <- NULL
    emax.e20 <- NULL
    emin.e11 <- NULL
    emax.e21 <- NULL
    emin.e12 <- NULL
    emax.e22 <- NULL

    for (l in 1:n) {
      eo = sort(re[l, ])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }

    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb <- apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x
    df <- data.frame(r=r, xab=a ,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td)) + geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20), fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=10,family=font))
  }
}

