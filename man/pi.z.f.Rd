\name{pi.z.f}
\alias{pi.z.f}
\title{Scaled density for a variety of distributions}
\usage{
  pi.z.f(g.type, pars, z, z.lim, seabedAtObs = NULL,
    n = NULL, mnormPlotFlag = FALSE,
    mnormVerboseFlag = FALSE)
}
\arguments{
  \item{g.type}{gradient type (Normal etc).
  c('NORM','BETA','LOGNORM','UNIFORM','TWEEDIE','MNORM')}

  \item{pars}{parameter vector containing only the gradient
  distribution parameters. Within the likelihood this is
  provided by par.unpack.F}

  \item{z}{depth}

  \item{z.lim}{vector (min,max) seabed or dimension of
  interest.}

  \item{seabedAtObs}{for multibeam (krill) only. this
  allows density to be rescaled on an}

  \item{n}{number of multivariate normal distributions}

  \item{mnormPlotFlag}{argument to be passed mnorm.f
  (diagnostic) - should mvnorm be plotted}

  \item{mnormVerboseFlag}{argument to be passed mnorm.f
  (diagnostic)}
}
\value{
  scaled density for a given depth (dimension of interest).
  CALLS: tweedief; mnorm.f
}
\description{
  Scaled density for a variety of distributions
}
\details{
  the scaling of density could be removed for the by
  observation calculation of the integration grid.  that
  way the integration grid could be calculated for each
  parameter estimate, then trimmed to seabed depth and then
  scaled. variable maxz can be overwritten by seabed at
  each observation in multibeam case
}

