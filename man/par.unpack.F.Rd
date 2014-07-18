\name{par.unpack.F}
\alias{par.unpack.F}
\title{Unpack parameters so that optim can run using a single parameter vector.}
\usage{
  par.unpack.F(grad.type, pars, n)
}
\arguments{
  \item{grad.type}{gradient type (Normal etc).
  c('NORM','BETA','LOGNORM','UNIFORM','TWEEDIE','MNORM')}

  \item{pars}{parameter vector containing all parameters to
  be estimated by optim.  pars vector NB structure for
  MNORM
  \eqn{\mu_1,\sigma_1,\alpha_1},...,\eqn{\mu_n,\sd_n},then
  detection function parameters.  there is no alpha_n.}

  \item{n}{number of distributions in a multinomial
  distribution (default NULL).}
}
\value{
  list([1]=density gradient parameters; [2]=detection
  function parameters).
}
\description{
  Unpack parameters so that optim can run using a single
  parameter vector.
}

