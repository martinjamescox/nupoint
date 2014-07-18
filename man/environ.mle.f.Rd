\name{environ.mle.f}
\alias{environ.mle.f}
\title{Log-likelihood forenvironmental variable (depth in case of whale) NB this function is optimised by optim.}
\usage{
  environ.mle.f(pars, z, rd, dzdy, z.mat, dzdy.mat, rd.mat,
    minz, wx, wy, wz, g.type, det.type, verbose = FALSE,
    n = NULL, diagnostic = TRUE)
}
\arguments{
  \item{pars}{model parameters, paramter vector optimised
  using optim. e.g.: pars[1]=depth gradient par 1;
  pars[2]=depth gradient par 2; pars[3]=half normal rng
  detection function. NB (i) pars unpacked using
  par.unpack.F, (ii) structure for MNORM
  \eqn{\mu_1,\sigma_1,\alpha_1,...,\mu_n,\sd_n}, there is
  no \eqn{alpha_n}}

  \item{z}{seabed depth under cue sighting.}

  \item{dzdy}{rate of change of seabed depth at cue
  sighting.}

  \item{rd}{radial distance from observer to cue.}

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

  \item{diagnostic}{print diagnostic data from the
  hore.get.numerator and shore.get.denominator.}

  \item{verbose}{prints the log-likelihood and parameter
  estimates for s,}
}
\value{
  single value of the log-likelihood denominator.
}
\description{
  Log-likelihood forenvironmental variable (depth in case
  of whale) NB this function is optimised by optim.
}

