#'@name diag.RBS
#'
#'@aliases diag.RBS
#'
#'@title Diagnostic Analysis - Local Influnce
#'
#'@description Diagnostics for the RBS model
#'
#'@param model Object of class \code{gamlss} holding the fitted model.
#'@param scheme Default is "case.weight". But, can be "response", "location" or "precision".
#'@param mu.link  Defines the mu.link, with "identity" link as the default for the mu parameter.
#'@param sigma.link Defines the sigma.link, with "identity" link as the default for the sigma parameter
#'@param lx Used in the scheme 'location'.
#'@param lz Used in the scheme 'precision'.
#'
#'@return Local influence measures.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@examples
#'data(cpd,package='faraway')
#'attach(cpd)
#'fit = gamlss::gamlss(actual ~ 0+projected,
#'family=RBS(mu.link="identity",sigma.link="identity"),method=CG())
#'Cib <- diag.RBS(fit)$Ci.beta
#'plot(Cib,ylim=c(0,1),pch=19,ylab=expression(C[i](beta)),xlab="Index")
#'abline(h=2*mean(Cib),lty=2)
#'Cia <- diag.RBS(fit)$Ci.alpha
#'plot(Cia,ylim=c(0,1),pch=19,ylab=expression(C[i](beta)),xlab="Index")
#'abline(h=2*mean(Cia),lty=2)
#'
#'
#'@importFrom pracma hessian ones
#'@importFrom stats make.link sd
#'@export

diag.RBS=function(model,mu.link = "identity",sigma.link = "identity",scheme="case.weight",lx=NULL,lz=NULL)
{

  x <- model$mu.x
  z  <- model$sigma.x
  y <- model$y
  p<-ncol(x)
  q<-ncol(z)

  linkstr <- mu.link
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta  <- linkobj$mu.eta

  sigma_linkstr <- sigma.link
  sigma_linkobj <- make.link(sigma_linkstr)
  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta  <- sigma_linkobj$mu.eta


  B=function(Delta,I,M)
  {
    B=(t(Delta)%*%(I-M)%*%Delta)
    return(B)
  }

  loglik <- function(vP)
  {
    betab = vP[1:p]
    alpha = vP[-(1:p)]
    eta   = as.vector(x%*%betab)
    tau   = as.vector(z%*%alpha)
    mu    = linkinv(eta)
    sigma = sigma_linkinv(tau)


    f <- 0.5*sigma -0.5*log((sigma+1)) - 0.5*log(mu) - 1.5*log(y) + log((sigma*y) + y + (sigma*mu)) - y*(sigma+1)/(4*mu) - (sigma*sigma*mu)/(4*y*(sigma+1)) - 0.5*log(16*pi)
    return(sum(f))
  }

  muest <- model$mu.coefficients
  sigmaest <- model$sigma.coefficients
  x0<- c(muest,sigmaest)
  h0 <- hessian(loglik,x0)

  Ldelta= h0[(p+1):(p+q),(p+1):(p+q)]
  Lbeta=h0[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, q))
  b12=cbind(matrix(0, q, p), solve(Ldelta))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(solve(Lbeta), matrix(0, p, q))
  b212= cbind(matrix(0, q, p), matrix(0, q, q))
  B2=rbind(b211,b212)  # parameter delta

  b311 =cbind(matrix(0, p, p), matrix(0, p, q))
  b312= cbind(matrix(0, q, p), matrix(0, q, q))
  B3=rbind(b311,b312)  # parameter theta

  if(scheme=="case.weight")
  {
  ############################Case Weight####################################

  mu <- model$mu.fv
  sigma <- model$sigma.fv
  eta <- linkfun(mu)
  ai <- mu.eta(eta)
  dmu <- (-1/(2*mu)) + sigma/((y*sigma) + y + (sigma*mu)) +  ((sigma+1)*y)/(4*(mu^2)) - (sigma^2)/(4*y*(sigma+1)) #ok!
  Deltamu <- crossprod(x,diag(ai*dmu))


  tau <- sigma_linkfun(sigma)
  bi <- sigma_mu.eta(tau)
  dsigma <- (y+ mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (sigma*(sigma+2)*mu)/(4*(sigma+1)*(sigma+1)*y) + sigma/(2*(sigma+1))
  Deltasigma <- crossprod(z,diag(bi*dsigma))

  Delta <- rbind(Deltamu,Deltasigma)

  ##################theta#########################
  BT<-B(Delta,solve(h0),B3)
  autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
  vetorpcthetaPC<- eigen(BT,symmetric=TRUE)$vec[,1]
  dmaxG.theta<-abs(vetorpcthetaPC)
  vCithetaPC<-2*abs(diag(BT))
  Cb0<-vCithetaPC
  Cb.theta<-Cb0/sum(Cb0)
  ######################betas########################
  BM<-B(Delta,solve(h0),B1)
  autovmaxbetaPC<-eigen(BM,symmetric=TRUE)$val[1]
  vetorpcbetaPC<-eigen(BM,symmetric=TRUE)$vec[,1]
  dmaxG.beta<-abs(vetorpcbetaPC)
  vCibetaPC<-2*abs(diag(BM))
  Cb1<-vCibetaPC
  Cb.beta<-Cb1/sum(Cb1)
  ####################alphas#########################
  BD<-B(Delta,solve(h0),B2)
  autovmaxdeltaPC<-eigen(BD,symmetric=TRUE)$val[1]
  vetordeltaPC<-eigen(BD,symmetric=TRUE)$vec[,1]
  dmaxG.alpha=abs(vetordeltaPC)
  vCideltaPC=2*abs(diag(BD))
  Cb2=vCideltaPC
  Cb.alpha=Cb2/sum(Cb2)

  result <- list(dmax.beta = dmaxG.beta,
                 dmax.alpha = dmaxG.alpha,
                 dmax.theta = dmaxG.theta,
                 Ci.beta = Cb.beta,
                 Ci.alpha = Cb.alpha,
                 Ci.theta = Cb.theta)
  return(result)
  }

  if(scheme=="response")
  {
  ############################Response####################################
  mu <- model$mu.fv
  sigma <- model$sigma.fv
  eta <- linkfun(mu)
  ai <- mu.eta(eta)
  tau <- sigma_linkfun(sigma)
  bi <- sigma_mu.eta(tau)
  phi<- ((2*sigma)+5)/((sigma+1)^2)
  sy<- sqrt((mu^2)*phi)

  dymu <- -(sigma*(sigma+1))/( ((sigma*y)  + y  + (sigma*mu))^2) + (sigma+1)/(4*(mu^2)) +  ((sigma^2)/(4*(sigma+1)*(y^2)))
  Deltamu <- crossprod(x,diag(ai*dymu*sy))
  p<-ncol(x);q<-ncol(z)
  dysigma <- (-mu/( ((sigma*y) + y + (sigma*mu))^2))  -  1/(4*mu) + (sigma*(sigma+2)*mu)/(4*(y^2)*((sigma+1)^2))
  Deltasigma <- crossprod(z,diag(bi*dysigma*sy))
  Delta <- rbind(Deltamu,Deltasigma)

  ###############thetas###########################
  BT<-B(Delta,solve(h0),B3)
  autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
  vetorthetaRP<- eigen(BT,symmetric=TRUE)$vec[,1]
  dmaxG.theta<-abs(vetorthetaRP)
  vCithetaRP<-2*abs(diag(BT))
  Cb0<-vCithetaRP
  Cb.theta<-Cb0/sum(Cb0)

  #################betas##########################
  BM=B(Delta,solve(h0),B1)
  autovmaxbetaRP <- eigen(BM,symmetric=TRUE)$val[1]
  vetorbetaRP <- eigen(BM,symmetric=TRUE)$vec[,1]
  dmaxG.beta <- abs(vetorbetaRP)
  vCibetaRP <- 2*abs(diag(BM))
  Cb1 <- vCibetaRP
  Cb.beta <- Cb1/sum(Cb1)
  ####################alpha#######################
  BD=B(Delta,solve(h0),B2)
  autovmaxdeltaRP <- eigen(BD,symmetric=TRUE)$val[1]
  vetordeltaRP <- eigen(BD,symmetric=TRUE)$vec[,1]
  dmaxG.alpha <- abs(vetordeltaRP)
  vCideltaRP <- 2*abs(diag(BD))
  Cb2 <- vCideltaRP
  Cb.alpha <- Cb2/sum(Cb2)


  result <- list(dmax.beta = dmaxG.beta,
                 dmax.alpha = dmaxG.alpha,
                 dmax.theta = dmaxG.theta,
                 Ci.beta = Cb.beta,
                 Ci.alpha = Cb.alpha,
                 Ci.theta = Cb.theta)
  return(result)
  }

  if(scheme=="location")
  {
  ############################Location Predictor####################################
  l <- lx
  mu <- model$mu.fv
  sigma <- model$sigma.fv
  eta <- linkfun(mu)
  ai <- mu.eta(eta)
  tau <- sigma_linkfun(sigma)
  bi <- sigma_mu.eta(tau)
  bl <- coef(model,what="mu")[l]
  sxl <- sd(x[,l])
  dmu <- (-1/(2*mu)) + sigma/((y*sigma) + y + (sigma*mu)) +  ((sigma+1)*y)/(4*(mu^2)) - (sigma^2)/(4*y*(sigma+1))
  ai1 <- rep(0,length(mu))
  dmu2 <- 1/(2*mu*mu) - (sigma^2)/(((y*sigma) + y + (sigma*mu))^2) - (y*(sigma+1))/(2*mu*mu*mu)
  ci <- dmu2*(ai^2) + dmu*ai1*ai
  xaux <- matrix(0,nrow=nrow(x),ncol=ncol(x))
  xaux[,l] <- ones(nrow(x),1)
  x0 <- xaux

  Deltamu <- bl*sxl*crossprod(x,diag(ci)) + sxl*crossprod(x0,diag(dmu*ai))

  dmusigma <- y/(((y*sigma) + y + (sigma*mu))^2) + y/(4*mu*mu) - (sigma*(sigma+2))/(4*(sigma+1)*(sigma+1)*y)
  mi <- dmusigma*bi*ai
  Deltasigma <- bl*sxl*crossprod(z,diag(mi))

  Delta <- rbind(Deltamu,Deltasigma)


  ###############thetas###########################
  BT<-B(Delta,solve(h0),B3)
  autovmaxthetaPX<- eigen(BT,symmetric=TRUE)$val[1]
  vetorthetaPX<- eigen(BT,symmetric=TRUE)$vec[,1]
  dmaxG.theta<-abs(vetorthetaPX)
  vCithetaPX<-2*abs(diag(BT))
  Cb0<-vCithetaPX
  Cb.theta<-Cb0/sum(Cb0)

  #################betas##########################
  BM=B(Delta,solve(h0),B1)
  autovmaxbetaPX <- eigen(BM,symmetric=TRUE)$val[1]
  vetorbetaPX <- eigen(BM,symmetric=TRUE)$vec[,1]
  dmaxG.beta <- abs(vetorbetaPX)
  vCibetaPX <- 2*abs(diag(BM))
  Cb1 <- vCibetaPX
  Cb.beta <- Cb1/sum(Cb1)
  ####################alpha#######################
  BD=B(Delta,solve(h0),B2)
  autovmaxdeltaPX <- eigen(BD,symmetric=TRUE)$val[1]
  vetordeltaPX <- eigen(BD,symmetric=TRUE)$vec[,1]
  dmaxG.alpha <- abs(vetordeltaPX)
  vCideltaPX <- 2*abs(diag(BD))
  Cb2 <- vCideltaPX
  Cb.alpha <- Cb2/sum(Cb2)

  result <- list(dmax.beta = dmaxG.beta,
                 dmax.alpha = dmaxG.alpha,
                 dmax.theta = dmaxG.theta,
                 Ci.beta = Cb.beta,
                 Ci.alpha = Cb.alpha,
                 Ci.theta = Cb.theta)
  return(result)
  }

  if(scheme=="precision")
  {
  ############################Precision Predictor####################################
  mu <- model$mu.fv
  sigma <- model$sigma.fv
  eta <- linkfun(mu)
  ai <- mu.eta(eta)
  tau <- sigma_linkfun(sigma)
  bi <- sigma_mu.eta(tau)
  k <- lz
  ak <- coef(model,what="sigma")[k]
  szk <- sd(z[,k])
  dsigma <- (y+ mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (sigma*(sigma+2)*mu)/(4*(sigma+1)*(sigma+1)*y) + sigma/(2*(sigma+1))
  Deltasigma <- crossprod(z,diag(bi*dsigma))

  bi1 <- rep(2,length(sigma))
  dsigma2 <- 1/(2*(sigma+1)*(sigma+1)) - ((y+mu)^2)/(((y*sigma) + y + (sigma*mu))^2) - mu/(2*(sigma+1)*(sigma+1)*(sigma+1)*y)
  wi <- dsigma2*(bi^2) + dsigma*bi1*bi
  zaux <- matrix(0,nrow=nrow(z),ncol=ncol(z))
  zaux[,k] <- ones(nrow(z),1)
  z0 <- zaux

  Deltasigma <- ak*szk*crossprod(z,diag(wi)) + szk*crossprod(z0,diag(dsigma*bi))

  dmusigma <- y/(((y*sigma) + y + (sigma*mu))^2) + y/(4*mu*mu) - (sigma*(sigma+2))/(4*(sigma+1)*(sigma+1)*y)
  mi <- dmusigma*bi*ai
  Deltamu <- ak*szk*crossprod(x,diag(mi))

  Delta <- rbind(Deltamu,Deltasigma)


  ###############thetas###########################
  BT<-B(Delta,solve(h0),B3)
  autovmaxthetaPZ<- eigen(BT,symmetric=TRUE)$val[1]
  vetorthetaPZ<- eigen(BT,symmetric=TRUE)$vec[,1]
  dmaxG.theta<-abs(vetorthetaPZ)
  vCithetaPZ<-2*abs(diag(BT))
  Cb0<-vCithetaPZ
  Cb.theta<-Cb0/sum(Cb0)

  #################betas##########################
  BM=B(Delta,solve(h0),B1)
  autovmaxbetaPZ <- eigen(BM,symmetric=TRUE)$val[1]
  vetorbetaPZ <- eigen(BM,symmetric=TRUE)$vec[,1]
  dmaxG.beta <- abs(vetorbetaPZ)
  vCibetaPZ <- 2*abs(diag(BM))
  Cb1 <- vCibetaPZ
  Cb.beta <- Cb1/sum(Cb1)
  ####################alpha#######################
  BD=B(Delta,solve(h0),B2)
  autovmaxdeltaPZ <- eigen(BD,symmetric=TRUE)$val[1]
  vetordeltaPZ <- eigen(BD,symmetric=TRUE)$vec[,1]
  dmaxG.alpha <- abs(vetordeltaPZ)
  vCideltaPZ <- 2*abs(diag(BD))
  Cb2 <- vCideltaPZ
  Cb.alpha <- Cb2/sum(Cb2)

  result <- list(dmax.beta = dmaxG.beta,
                 dmax.alpha = dmaxG.alpha,
                 dmax.theta = dmaxG.theta,
                 Ci.beta = Cb.beta,
                 Ci.alpha = Cb.alpha,
                 Ci.theta = Cb.theta)
  return(result)
  }
  ############################Joint Predictor####################################

  }

#'@name diag.ZARBS
#'
#'
#'@aliases dias.ZARBS
#'
#'
#'@title Diagnostic Analysis - Local Influence
#'
#'@description Diagnostics for the ZARBS model
#'
#'@param model object of class \code{gamlss} holding the fitted model.
#'@param links link function used.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@export
#'@importFrom pracma hessian

diag.ZARBS <- function(model,links=c("log","identity","probit"))
{
  x <- model$mu.x
  z <- model$sigma.x
  w <- model$nu.x
  y <- model$y
  p <-ncol(x)
  q <-ncol(z)
  s <-ncol(w)

  linkstr <- links[1]
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta  <- linkobj$mu.eta

  sigma_linkstr <- links[2]
  sigma_linkobj <- make.link(sigma_linkstr)
  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta  <- sigma_linkobj$mu.eta

  nu_linkstr <- links[3]
  nu_linkobj <- make.link(nu_linkstr)
  nu_linkfun <- nu_linkobj$linkfun
  nu_linkinv <- nu_linkobj$linkinv
  nu_mu.eta  <- nu_linkobj$mu.eta


  B <- function(Delta,I,M)
  {
    B <- (t(Delta)%*%(I-M)%*%Delta)

    B
  }

  loglik <- function(vP,y){
    betab <- vP[1:p]
    alpha <- vP[(p+1):(p+q)]
    gaMMa <- vP[(p+q+1):(p+q+s)]
    eta   <- as.vector(x%*%betab)
    tau   <- as.vector(z%*%alpha)
    xi    <- as.vector(w%*%gaMMa)
    mu    <- linkinv(eta)
    sigma <- sigma_linkinv(tau)
    nu    <- nu_linkinv(xi)
    f <- ifelse(y==0,log(nu),(log(1-nu)+ 0.5*sigma -0.5*log((sigma+1)) - 0.5*log(mu) - 1.5*log(y) + log((sigma*y) + y + (sigma*mu)) - y*(sigma+1)/(4*mu) - (sigma*sigma*mu)/(4*y*(sigma+1)) - 0.5*log(16*pi)))
    sum(f)
  }

  betaest <- model$mu.coefficients
  alphaest <- model$sigma.coefficients
  gammaest <- model$nu.coefficients
  x0<- c(betaest,alphaest,gammaest)
  h0 <- hessian(loglik,x0,y=y)

  Lbeta <- h0[1:p,1:p]
  Lalpha <-  h0[(p+1):(p+q),(p+1):(p+q)]
  Lbetaalpha <- h0[1:(p+1),1:(p+1)]
  LgaMMa <- h0[(p+q+1):(p+q+s),(p+q+1):(p+q+s)]

  B1 <- matrix(0,p+q+s,p+q+s)
  B1[(p+q):(p+q+s-1),(p+q):(p+q+s-1)] <- -solve(Lalpha)
  B1[(p+q+1):(p+q+s),(p+q+1):(p+q+s)] <- -solve(LgaMMa)

  B2 <- matrix(0,p+q+s,p+q+s)
  B2[(1:p),1:p] <- -solve(Lbeta)
  B2[(p+q+1):(p+q+s),(p+q+1):(p+q+s)] <- -solve(LgaMMa)

  B3 <- matrix(0,p+q+s,p+q+s)
  B3[1:(p+1),1:(p+1)] <- -solve(Lbetaalpha)

  B4 <- matrix(0,p+q+s,p+q+s)

  mu <- model$mu.fv
  sigma <- model$sigma.fv
  nu <- model$nu.fv
  ki <- (1-(y==0))
  eta <- linkfun(mu)
  tau <- sigma_linkfun(sigma)
  xi <- nu_linkfun(nu)
  ai <- mu.eta(eta)
  bi <- sigma_mu.eta(tau)
  ci <- nu_mu.eta(tau)


  dmu <- ifelse(y>0,(-1/(2*mu)) + sigma/((y*sigma) + y + (sigma*mu)) +  ((sigma+1)*y)/(4*(mu^2)) - (sigma^2)/(4*y*(sigma+1)),0) #ok!
  kdmu <- ki*dmu
  Delta.beta <- crossprod(x,diag(ai*kdmu))



  dsigma <- ifelse(y>0,(y+ mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (sigma*(sigma+2)*mu)/(4*(sigma+1)*(sigma+1)*y) + sigma/(2*(sigma+1)),0)
  kdsigma <- ki*dsigma
  Delta.alpha <- crossprod(z,diag(bi*kdsigma))


  dnu <- ifelse(y==0,(1/(nu*(1-nu))),0) - 1/(1-nu)
  Delta.nu <- crossprod(w,diag(ci*dnu))

  Delta <- rbind(Delta.beta,Delta.alpha,Delta.nu)

  ##################theta#########################
  BT <- B(Delta,solve(h0),B4)
  autovmaxthetaPC <- eigen(BT,symmetric=TRUE)$val[1]
  vetorpcthetaPC <- eigen(BT,symmetric=TRUE)$vec[,1]
  dmaxG.theta <- abs(vetorpcthetaPC)
  vCithetaPC <- 2*abs(diag(BT))
  Cb0 <- vCithetaPC
  Cb.theta <-Cb0/sum(Cb0)
  ######################beta's########################
  BM<-B(Delta,solve(h0),B1)
  autovmaxbetaPC<-eigen(BM,symmetric=TRUE)$val[1]
  vetorpcbetaPC<-eigen(BM,symmetric=TRUE)$vec[,1]
  dmaxG.beta<-abs(vetorpcbetaPC)
  vCibetaPC<-2*abs(diag(BM))
  Cb1<-vCibetaPC
  Cb.beta<-Cb1/sum(Cb1)
  ####################alpha's#########################
  BD<-B(Delta,solve(h0),B2)
  autovmaxdeltaPC<-eigen(BD,symmetric=TRUE)$val[1]
  vetordeltaPC<-eigen(BD,symmetric=TRUE)$vec[,1]
  dmaxG.alpha=abs(vetordeltaPC)
  vCideltaPC=2*abs(diag(BD))
  Cb2=vCideltaPC
  Cb.alpha=Cb2/sum(Cb2)
  ####################gamma's#########################
  BP<-B(Delta,solve(h0),B3)
  autovmaxgammaPC<-eigen(BP,symmetric=TRUE)$val[1]
  vetorgammaPC<-eigen(BP,symmetric=TRUE)$vec[,1]
  dmaxG.gamma=abs(vetorgammaPC)
  vCigammaPC=2*abs(diag(BP))
  Cb3=vCigammaPC
  Cb.gamma=Cb3/sum(Cb3)


  result <- list(dmax.beta = dmaxG.beta,
                 dmax.alpha = dmaxG.alpha,
                 dmax.gamma = dmaxG.gamma,
                 dmax.theta = dmaxG.theta,
                 Ci.beta = Cb.beta,
                 Ci.alpha = Cb.alpha,
                 Ci.gamma = Cb.gamma,
                 Ci.theta = Cb.theta)
  return(result)
}




