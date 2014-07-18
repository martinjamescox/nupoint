\name{shore.get.numerator.mbe}
\alias{shore.get.numerator.mbe}
\title{Calculate log-likelihood numerator}
\usage{
  shore.get.numerator.mbe(pars, sight.x, swarm.depth, w,
    seabed.z, grad.type, det.type, n = NULL, angularDetect,
    x.max)
}
\arguments{
  \item{pars}{model parameters, pars[1]=depth gradient par
  1; pars[2]=depth gradient par 2; pars[3]=half normal rng
  detection function}

  \item{sight.x}{observation vector: target x coordinate.
  This is cross-track distance for multi-beam
  echosounder(MBE) observations.'}

  \item{swarm.depth}{depth at geometric swarm centre (m)}

  \item{w}{truncation distance (m)}

  \item{seabed.z}{seabed depth under geometric swarm centre
  (m)}

  \item{grad.type}{parametric form of the environmental
  gradient.  See environ.fit.f.}

  \item{det.type}{parametric form of the detection
  function.  See detectF.}

  \item{n}{number of normal mixture distributions in the
  MNORM function.}

  \item{angularDetect}{logical: FALSE. Should angular
  detection be considered (not yet implemented).}

  \item{x.max}{maximum distance in the x-dimension (w *
  sin(theta.max))}
}
\value{
  numerator for the log-likelihood
}
\description{
  Calculate log-likelihood numerator
}

