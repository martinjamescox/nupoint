\name{mle}
\alias{mle}
\title{Maximum likelihood estimation (wrapper for 'optim')}
\usage{
  mle(x, dist, start = NULL, method = "Nelder-Mead")
}
\arguments{
  \item{x}{data}

  \item{dist}{distribution being fitted}

  \item{start}{starting values for parameters}

  \item{method}{passed through to optim for optimisation
  method to use}
}
\value{
  list containing AIC values, VC matrix, parameter
  estimates, log likelihood...
}
\description{
  Maximum likelihood estimation (wrapper for 'optim')
}

