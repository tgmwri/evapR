% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainFun.R
\name{gimme_eq}
\alias{gimme_eq}
\title{Return model equation}
\usage{
gimme_eq(y, x, type = "lm", eqOnly = FALSE, message = TRUE)
}
\arguments{
\item{y}{dependent variable}

\item{x}{regressor}

\item{type}{type of trend to use: linear - "lm", exponential - "exp" and power trend "power" or "all". Default is "lm".}

\item{eqOnly}{default is FALSE, if TRUE, returns only equation without coefficients.}

\item{message}{whether to print equation to the console (default is TRUE).}
}
\description{
Function that returns model equation and it's coefficients.
}
\examples{
gimme_eq(est$evaporation, eat$radiation, type = "lm", eqOnly = T)
}
