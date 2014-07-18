\name{shore.get.denominator}
\alias{shore.get.denominator}
\title{Denominator of the environmental variable (depth in case of whale) log-likelihood. NB this is an internal function.}
\usage{
  shore.get.denominator(pars, num.detects, wx, wy, minz,
    wz, grad.type, z.mat, dzdy.mat, rd.mat, det.type,
    n = NULL, diagnostic = FALSE)
}
\arguments{
  \item{pars}{model parameters, paramter vector optimised
  using optim. e.g.: pars[1]=depth gradient par 1;
  pars[2]=depth gradient par 2; pars[3]=half normal rng
  detection function. NB pars unpacked using par.unpack.F}

  \item{num.detects}{number of detections.}

  \item{z.mat}{gridded seabed depth throughout survey area
  (regular grid; positive depths).}

  \item{dzdy.mat}{rate of change of seabed depth throughout
  survey area.}

  \item{rd.mat}{matrix of radial distances from observer to
  each grid cell in survey area.}

  \item{minz}{minimum z (depth) NB assumes z is positive.}

  \item{wx}{x-dimension truncation distance (for
  x-direction uniform sampling).}

  \item{wy}{y-dimension truncation distance.}

  \item{wz}{seabed depth truncation distance. NB assumes z
  is positive}

  \item{grad.type}{cue depth disribution type e.g. "NORM" -
  this is currently coded for a radial detection function
  only}

  \item{n}{number of distributions in a multinomial
  distribution (default NULL).}

  \item{diagnostic}{print diagnostic data.}
}
\value{
  single value of the log-likelihood denominator.
}
\description{
  Denominator of the environmental variable (depth in case
  of whale) log-likelihood. NB this is an internal
  function.
}

