% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resid.R
\name{residuals}
\alias{residuals}
\title{Residuals}
\usage{
residuals(model,residual = "deviance")
}
\arguments{
\item{model}{an object for which the extraction of model residuals is meaningful.}

\item{residual}{type of residual to be used ("pearson","score", "deviance" or "quantile"). Default is the deviance.}
}
\value{
Residuals extracted from the object object.
}
\description{
\code{residuals} is a function which extracts model residuals from objects returned by modeling functions.
}
\references{
Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.

Santos-Neto, M., Cysneiros, F.J.A., Leiva, V., Barros, M. (2016) Reparameterized Birnbaum-Saunders
regression models with varying precision. \emph{Electronic Journal of Statistics}, 10, 2825--2855. doi: \email{10.1214/16-EJS1187}.
}
\author{
Manoel Santos-Neto \url{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \url{cysneiros@de.ufpe.br}, Victor Leiva \url{victorleivasanchez@gmail.com} and Michelli Barros \url{michelli.karinne@gmail.com}
}
