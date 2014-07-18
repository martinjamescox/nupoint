\name{shore.get.denominator.mbe}
\alias{shore.get.denominator.mbe}
\title{Numerical integration for likelihood denominator}
\usage{
  shore.get.denominator.mbe(pars, num.detects, w,
    grad.type, det.type, n = NULL, mask, xV, yV,
    attenuation, angularDetect, xmax)
}
\arguments{
  \item{pars}{model parameters, pars[1]=depth gradient par
  1; pars[2]=depth gradient par 2; pars[3]=half normal rng
  detection function}

  \item{num.detects}{number of animals detected}

  \item{w}{truncation distance (y-dimension)}

  \item{grad.type}{parametric form of the environmental
  gradient.  See environ.fit.f.}

  \item{det.type}{parametric form of the detection
  function.  See detectF.}

  \item{n}{number of normal mixture distributions in the
  MNORM function.}

  \item{mask}{area within the point transect where
  observation is possible e.g. MBE swath. created using
  swathInOutF in mbe.fit.f}

  \item{xV}{vector of x dimension coordinates for
  integration grid.}

  \item{yV}{vector of y dimension coordinates for
  integration grid.}

  \item{attenuation}{Logical TRUE - use attenuation
  function (NB attenuation function must be called
  atten.f); FALSE - no attenuation.}

  \item{angularDetect}{logical: FALSE. Should angular
  detection be considered (not yet implemented).}
}
\value{
  likelihood denominator.
}
\description{
  Numerical integration for likelihood denominator
}

