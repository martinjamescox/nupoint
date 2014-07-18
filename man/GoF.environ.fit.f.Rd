\name{GoF.environ.fit.f}
\alias{GoF.environ.fit.f}
\title{Chi-sq GoF for whale data (this is the 1D version)}
\usage{
  GoF.environ.fit.f(pars, r.mat, z.mat, minz, wz, z.obs,
    grad.type, det.type, nbr.intervals, plot = FALSE,
    n = NULL, dzdy.mat, breaks = NULL)
}
\arguments{
  \item{pars}{= model parameters, paramter vector optimised
  using optim. e.g.:}

  \item{r.mat}{radial distance from observer to centre of
  each integration grid cell.}

  \item{z.mat}{seabed depth at each integration grid.}

  \item{minz}{minimum z (depth) NB assumes z is positive.}

  \item{wz}{seabed depth truncation distance. NB assumes z
  is positive}

  \item{z.obs}{vector of depths at each observed cue,}

  \item{grad.type}{cue depth disribution type e.g. "NORM" -
  this is currently coded for a radial detection function
  only}

  \item{det.type}{detection function form.}

  \item{nbr.intervals}{number of intervals for use in GoF.}

  \item{n}{number of distributions in a multinomial
  distribution (default NULL).}

  \item{plot}{TRUE/FALSE 4 panel plot of GoF}

  \item{dzdy.mat}{matrix of derivatives of z with respct to
  y}
}
\value{
  single value of the log-likelihood numerator.
}
\description{
  Chi-sq GoF for whale data (this is the 1D version)
}
\details{
  Calls par.unpack.F; pi.z.f; package fields for plotting.
}
\note{
  Requires fields package, functions: histline; detectF;
  pi.z.f
}

