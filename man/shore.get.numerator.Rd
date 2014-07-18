\name{shore.get.numerator}
\alias{shore.get.numerator}
\title{Numerator of the environmental variable (depth in case of whale) log-likelihood. NB this is an internal function.}
\usage{
  shore.get.numerator(pars, z, dzdy, rd, minz, wx, wz,
    grad.type, det.type, n = NULL, diagnostic = FALSE)
}
\arguments{
  \item{pars}{model parameters, paramter vector optimised
  using optim. e.g.: pars[1]=depth gradient par 1;
  pars[2]=depth gradient par 2; pars[3]=half normal rng
  detection function. NB pars unpacked using par.unpack.F}

  \item{z}{seabed depth under cue sighting.}

  \item{dzdy}{rate of change of seabed depth at cue
  sighting.}

  \item{rd}{radial distance from observer to cue.}

  \item{minz}{minimum z (depth) NB assumes z is positive.}

  \item{wx}{x-dimension truncation distance (for
  x-direction uniform sampling).}

  \item{wz}{seabed depth truncation distance. NB assumes z
  is positive}

  \item{grad.type}{cue depth disribution type e.g. "NORM" -
  this is currently coded for}

  \item{a}{radial detection function only}

  \item{det.type}{radial detection function form.}

  \item{n}{number of distributions in a multinomial
  distribution (default NULL).}

  \item{diagnostic}{print diagnostic data.}
}
\value{
  single value of the log-likelihood numerator.
}
\description{
  Numerator of the environmental variable (depth in case of
  whale) log-likelihood. NB this is an internal function.
}

