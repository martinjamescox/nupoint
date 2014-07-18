\name{shore.only}
\alias{shore.only}
\title{Log-likelihood for the density gradient with respect to a linear feasture.}
\usage{
  shore.only(pars, sight.x, swarm.z, grid.density, td,
    sb.z, g.type, det.type, n = NULL, verbose = TRUE, mask,
    xV, yV, angularDetect)
}
\arguments{
  \item{pars}{a single vector of starting parameters to be
  estimated.  Parameters must be ordered: environment
  density gradient first, followed by detection function.
  See environ.fit.f'}

  \item{sight.x}{observation vector: target x coordinate.
  This is cross-track distance for multi-beam
  echosounder(MBE) observations.'}

  \item{swarm.z}{observation vector: target y coordinate.
  This is water depth to the centre of a swarm for MBE
  observations.'}

  \item{grid.density}{numerical integration grid density
  (number of elements in each (x and y) dimensions.'}

  \item{td}{truncation distance in the y-dimension (depth
  for MBE)'}

  \item{sb.z}{observation vector: default NULL if no
  attenuation function is required (all transects extend to
  w). FOr MBE observations, This is water depth to the
  centre of a swarm.'}

  \item{g.type}{parametric form of the environmental
  gradient.  See environ.fit.f.'}

  \item{det.type}{parametric form of the detection
  function.  See detectF.'}

  \item{n}{number of normal mixture distributions in the
  MNORM function.}

  \item{angularDetect}{logical: FALSE. Should angular
  detection be considered (not yet implemented).'}

  \item{verbose}{Logical. Default FALSE. If TRUE parameter
  estimates and log-likelihood are printed to the console
  for each optim search iteration.'}

  \item{mask}{area within the point transect where
  observation is possible e.g. MBE swath. created using
  swathInOutF in mbe.fit.f}
}
\value{
  -log lik.
}
\description{
  Log-likelihood for the density gradient with respect to a
  linear feasture.
}

