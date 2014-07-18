\name{shore.calc.integral.grid}
\alias{shore.calc.integral.grid}
\title{Calculate value at a given numerical integration cell in the likelihood denominator.}
\usage{
  shore.calc.integral.grid(x, y, pars, trunc.dist,
    grad.type, det.type, n = NULL, attenuation,
    angularDetect)
}
\arguments{
  \item{x}{and y coordinates}

  \item{pars}{model parameters, pars[1]=depth gradient par
  1; pars[2]=depth gradient par 2; pars[3]=half normal rng
  detection function}

  \item{trunc.dist}{truncation distance (y-dimension)}

  \item{grad.type}{parametric form of the environmental
  gradient.  See environ.fit.f.}

  \item{det.type}{parametric form of the detection
  function.  See detectF.}

  \item{n}{number of normal mixture distributions in the
  MNORM function.}

  \item{attenuation}{Logical TRUE - use attenuation
  function (NB attenuation function must be called
  atten.f); FALSE - no attenuation.}

  \item{angularDetect}{logical: FALSE. Should angular
  detection be considered (not yet implemented).}
}
\value{
  value at denominator numerical integration grid cell.
}
\description{
  Calculate value at a given numerical integration cell in
  the likelihood denominator.
}

