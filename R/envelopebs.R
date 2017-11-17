#' Simulated Evelope
#' @description A normal plot with simulated envelope of the residual is produced.
#'
#'@usage envelope(model,k=19,alpha=0.05,res="deviance", precision = c("fixed","varying"),
#' dist = RBS(mu.link = "identity",sigma.link = "identity"),...)
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param alpha value giving the size of the envelope. Default is 0.05 which is equivalent to a 95\% band.
#' @param res type of residuals to be extracted. Default is deviance.
#' @param precision If \code{precision = "fixed"} a model with precision fixed is used;
#' else a model with precision variable is used.
#' @param dist The function RBS() defines the RBS distribution.
#'
#'@return A simulated envelope of the class RBS.
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
#'
#' @examples
#'
#'## Fixed Precision
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'set.seed(2015)
#'envelope(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'set.seed(2015)
#'envelope(model)
#'
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
#'summary(fit0)
#'set.seed(2015)
#'envelope(fit0,alpha=0.01, precision="fixed",res="quantile",dist=RBS(mu.link="identity"))
#'##Varying Precision
#'fit1 <- gamlss(y1 ~ x21,sigma.formula = y1 ~x21, family=RBS(mu.link="identity",sigma.link="sqrt"),method=CG()  )
#'summary(fit11)
#'set.seed(2015)
#'envelope(fit1,alpha=0.01, precision="varying",res="quantile",dist=RBS(mu.link="identity",sigma.link="sqrt"))
#'
#'
#'
#'library(ssym)
#'data(Snacks,package="ssym")
#'head(Snacks)
#'attach(Snacks)
#'type <- factor(type,labels = c("A","B","C","D","E"))
#'w1 <- week
#'w2 <- I(week^2)
#'fit. <- gamlss(texture~type+w1+w2,~type, family=RBS(mu.link="log",sigma.link = "log"),method=CG())
#'summary(fit.)
#'envelope(fit.,dist=RBS(mu.link="log",sigma.link = "log"),precision = "varying",color="white",xlabel = "perc",ylabel = "rd",border.col = "black")
#' @export


envelope <- function(model,k=19,alpha=0.05,res="deviance", precision = c("fixed","varying"), dist = RBS(mu.link = "identity",sigma.link = "identity") ,color="gray80",xlabel="Theorical Quantile",
                     ylabel="Empirical Quantile",main="",pch=20, border.col = "gray")
{

  if(precision!= "fixed")
  {
    alfa1 <- ceiling(k*alpha)
    alfa2 <- ceiling(k*(1-alpha))
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

      if(p == 1) form <- nresp ~ 0 + X else form <- nresp ~ X[,-1]

      if(q == 1) form1 <- nresp ~ 0 + Z else  form1 <- nresp ~ Z[,-1]

      conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form,sigma.formula = form1 ,family=dist,method=CG(),control = conh0)
      rdf <- residuals(model1,residual=res)

      re[,i]=sort(rdf)
    }
    e1 = numeric(n)
    e2 = numeric(n)
    for(l in 1:n)
    {
      eo = sort(re[l,])
      e1[l] = eo[alfa1]
      e2[l] = eo[alfa2]
    }
    a<-  qqnorm(e1,plot.it=FALSE)$x
    a1<-  qqnorm(e1,plot.it=FALSE)$y
    b<-  qqnorm(e2,plot.it=FALSE)$x
    b1<-  qqnorm(e2,plot.it=FALSE)$y
    r<-  qqnorm(td,plot.it=FALSE)$x
    r1<-  qqnorm(td,plot.it=FALSE)$y
    xx <- c(a,rev(b))
    yy <- c(a1,rev(b1))

    xb = apply(re,1,mean)
    faixa = range(td,e1,e2)
    par(mar=c(4., 4., 0.1, 0.1))
    plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
    par(new=TRUE)
    polygon(xx,yy,col=color,border=NA)
    par(new=TRUE)
    qqnorm(e1,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(e2,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(xb,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=2,lwd=1,col="black")
    par(new=TRUE)
    qqnorm(td,xlab=xlabel,main=main,ylab=ylabel,ylim=faixa,pch=pch,cex=1,lwd=1)
  } else{
    alfa1 <- ceiling(k*alpha)
    alfa2 <- ceiling(k*(1-alpha))
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
      if(p == 1) form <- nresp ~ 0 + X else  form <- nresp ~ X[,-1];

      conh0 <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form ,family=dist,method=CG(),control = conh0)
      rdf <- residuals(model1,residual=res)
      re[,i] <- sort(rdf)
    }
    e1 = numeric(n)
    e2 = numeric(n)
    for(l in 1:n)
    {
      eo = sort(re[l,])
      e1[l] = eo[alfa1]
      e2[l] = eo[alfa2]
    }

    a<-  qqnorm(e1,plot.it=FALSE)$x
    a1<-  qqnorm(e1,plot.it=FALSE)$y
    b<-  qqnorm(e2,plot.it=FALSE)$x
    b1<-  qqnorm(e2,plot.it=FALSE)$y
    r<-  qqnorm(td,plot.it=FALSE)$x
    r1<-  qqnorm(td,plot.it=FALSE)$y
    xx <- c(a,rev(b))
    yy <- c(a1,rev(b1))

    xb = apply(re,1,mean)
    faixa = range(td,e1,e2)
    par(mar=c(4., 4., 0.1, 0.1))
    plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
    par(new=TRUE)
    polygon(xx,yy,col=color,border=NA)
    par(new=T)
    qqnorm(e1,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(e2,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(xb,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=2,lwd=1,col="black")
    par(new=TRUE)
    qqnorm(td,xlab=xlabel,main=main,ylab=ylabel,ylim=faixa,pch=pch,cex=1,lwd=1)
  }
}


#' Simulated Evelope
#' @description A normal plot with simulated envelope of the residual is produced.
#'
#'@usage envelope(model,k=19,alpha=0.05,res="deviance", precision = c("fixed","varying"),
#' dist = RBS(mu.link = "identity",sigma.link = "identity"),...)
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param alpha value giving the size of the envelope. Default is 0.05 which is equivalent to a 95\% band.
#' @param res type of residuals to be extracted. Default is deviance.
#' @param precision If \code{precision = "fixed"} a model with precision fixed is used;
#' else a model with precision variable is used.
#' @param dist The function RBS() defines the RBS distribution.
#'
#'@return A simulated envelope of the class RBS.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
#'regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
#'
#'Santos-Neto, M., Tomazella, V., Perreira, G.H., Nobre, J. (2017+) A general class of Birnbaum-Saunders regression models for data containing zeros.
#'
#'
#' @examples
#'
#'## Fixed Precision
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'set.seed(2015)
#'envelope(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'set.seed(2015)
#'envelope(model)
#'
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
#'summary(fit0)
#'set.seed(2015)
#'envelope(fit0,alpha=0.01, precision="fixed",res="quantile",dist=RBS(mu.link="identity"))
#'##Varying Precision
#'fit1 <- gamlss(y1 ~ x21,sigma.formula = y1 ~x21, family=RBS(mu.link="identity",sigma.link="sqrt"),method=CG()  )
#'summary(fit11)
#'set.seed(2015)
#'envelope(fit1,alpha=0.01, precision="varying",res="quantile",dist=RBS(mu.link="identity",sigma.link="sqrt"))
#'
#'
#'
#'library(ssym)
#'data(Snacks,package="ssym")
#'head(Snacks)
#'attach(Snacks)
#'type <- factor(type,labels = c("A","B","C","D","E"))
#'w1 <- week
#'w2 <- I(week^2)
#'fit. <- gamlss(texture~type+w1+w2,~type, family=RBS(mu.link="log",sigma.link = "log"),method=CG())
#'summary(fit.)
#'envelope(fit.,dist=RBS(mu.link="log",sigma.link = "log"),precision = "varying",color="white",xlabel = "perc",ylabel = "rd",border.col = "black")
#' @export

envelope.ZARBS <- function(model,k=19,alpha=0.05, precision = c("fixed","varying"), dist = ZARBS(mu.link = "identity",sigma.link = "identity",nu.link="logit"),color="gray80",xlabel="Theorical Quantile",ylabel="Empirical Quantile",main="",pch=20, border.col = "gray")
{
  if(precision!="fixed"){
    alfa1 <- ceiling(k*alpha)
    alfa2 <- ceiling(k*(1-alpha))
    n <- model$N
    td  <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    nu <- model$nu.fv
    re <- matrix(0,n,k)
    X <- model$mu.x
    Z <- model$sigma.x
    W <- model$nu.x
    p<- ncol(X)
    q<- ncol(Z)
    s<- ncol(W)

    for(i in 1:k)
    {

      y1 <- mapply(rZARBS,n=1,mu=mu,sigma=sigma,nu=nu)
      nresp <- y1

      if(p == 1){ form <- nresp ~ 0 + X} else{ form <- nresp ~ X[,-1]}

      if(q == 1){ form1 <- nresp ~ 0 + Z} else  form1 <- nresp ~ Z[,-1]

      if(s == 1){ form2 <- nresp ~ 0 + W} else  form2 <- nresp ~ W[,-1]

      conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form,sigma.formula = form1,nu.formula=form2,family=dist,method=RS(),control = conh0)
      rdf <- model1$residuals

      re[,i]=sort(rdf)
    }
    e1 = numeric(n)
    e2 = numeric(n)
    for(l in 1:n)
    {
      eo = sort(re[l,])
      e1[l] = eo[alfa1]
      e2[l] = eo[alfa2]
    }
    a<-  qqnorm(e1,plot.it=FALSE)$x
    a1<-  qqnorm(e1,plot.it=FALSE)$y
    b<-  qqnorm(e2,plot.it=FALSE)$x
    b1<-  qqnorm(e2,plot.it=FALSE)$y
    r<-  qqnorm(td,plot.it=FALSE)$x
    r1<-  qqnorm(td,plot.it=FALSE)$y
    xx <- c(a,rev(b))
    yy <- c(a1,rev(b1))

    xb = apply(re,1,mean)
    faixa = range(td,e1,e2)
    par(mar=c(4., 4., 0.1, 0.1))
    plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
    par(new=TRUE)
    polygon(xx,yy,col=color,border=NA)
    par(new=TRUE)
    qqnorm(e1,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(e2,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(xb,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=2,lwd=1,col="black")
    par(new=TRUE)
    qqnorm(td,xlab=xlabel,main=main,ylab=ylabel,ylim=faixa,pch=pch,cex=1,lwd=1)
  }
  else{
    alfa1 <- ceiling(k*alpha)
    alfa2 <- ceiling(k*(1-alpha))
    n<-model$N
    td  <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    nu <- model$nu.fv
    re <- matrix(0,n,k)
    X <- model$mu.x
    W <- model$nu.x
    p <- ncol(X)
    s <- ncol(W)
    mulink <- model$mu.link


    for(i in 1:k)
    {

      y1 <- mapply(rZARBS,n=1,mu=mu,sigma=sigma,nu=nu)
      nresp <- y1
      if(p == 1){form <- nresp~0+X} else{form <- nresp ~ X[,-1]}
      if(s == 1){form1 <- nresp ~ 0 + W} else{form1 <- nresp ~ W[,-1]}


      conh0 <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form,nu.formula=form1,family=dist,method=RS(),control = conh0)
      rdf <- model1$residuals
      re[,i] <- sort(rdf)
    }
    e1 = numeric(n)
    e2 = numeric(n)
    for(l in 1:n)
    {
      eo = sort(re[l,])
      e1[l] = eo[alfa1]
      e2[l] = eo[alfa2]
    }

    a<-  qqnorm(e1,plot.it=FALSE)$x
    a1<-  qqnorm(e1,plot.it=FALSE)$y
    b<-  qqnorm(e2,plot.it=FALSE)$x
    b1<-  qqnorm(e2,plot.it=FALSE)$y
    r<-  qqnorm(td,plot.it=FALSE)$x
    r1<-  qqnorm(td,plot.it=FALSE)$y
    xx <- c(a,rev(b))
    yy <- c(a1,rev(b1))

    xb = apply(re,1,mean)
    faixa = range(td,e1,e2)
    par(mar=c(4., 4., 0.1, 0.1))
    plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
    par(new=TRUE)
    polygon(xx,yy,col=color,border=NA)
    par(new=T)
    qqnorm(e1,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(e2,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col=border.col)
    par(new=TRUE)
    qqnorm(xb,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=2,lwd=1,col="black")
    par(new=TRUE)
    qqnorm(td,xlab=xlabel,main=main,ylab=ylabel,ylim=faixa,pch=pch,cex=1,lwd=1)
  }
}

