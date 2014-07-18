\name{tweedief}
\alias{tweedief}
\title{Density of a scaled Tweedie distribution.}
\usage{
  tweedief(x, pars, zlim)
}
\arguments{
  \item{x}{observation (seabed depth in the case of the
  whales; depth of swarm in case of multibeam)}

  \item{pars}{parameter vector for the distribution (does
  not include the detection function) pars[1] mu; pars[2]
  phi, and pars[3] power.}

  \item{zlim}{z-dimension (whale) or y-dimension
  (multi-beam) truncation distance - used to scale the
  density.}
}
\value{
  scaled density of the tweedie distribution for parameters
  par
}
\description{
  Density of a scaled Tweedie distribution.
}
\note{
  this function is required because passing an NA into
  dtweedie causes an error.  This function returns an NA
  when it is passed an NA in the x argument.
}

