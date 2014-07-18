\name{est.abundance.whales}
\alias{est.abundance.whales}
\title{Estimate animal abundance from actual data set}
\usage{
  est.abundance.whales(environ.sim.dat, trunc.prop = 0.9)
}
\arguments{
  \item{fit.obj}{fitted object}

  \item{truncation}{distance proportion (default 0.9) such
  that sightings beyond 0.9*max.r are deleted}
}
\value{
  list containing abundance estimate within covered region
  and abundance estimate for entire study area (assuming
  grid cells are unit square in area)
}
\description{
  Code as of 2 August 2013 is incomplete and not in a
  functional state.
}
\details{
  This function estimates animal abundance within the study
  area (grid) by calculating density \eqn{\pi (z(x,y))} as
  a function of covariate for each grid cell.

  Calls to appropriate distribution (normal, lognormal,
  beta, uniform, mixture of normals) in association with
  the parameters estimated by the likelihood routine
  (\code{nupoint.env.fit}) are summed to produce estimate.

  Should your grid cell sizes not be unit square, then
  multiply the value returned by this function by the grid
  cell size to produce abundance estimate in the units
  appropriate for your study.
}
\author{
  Eric Rexstad
}
\references{
  M.J. Cox, D.L. Borchers, D.A. Demer, G.R. Cutter, and
  A.S. Brierley. 2011. Estimating the density of Antarctic
  krill (Euphausia superba) from multi-beam echo-sounder
  observations using distance sampling methods. Journal of
  the Royal Statistical Society: Series C (Applied
  Statistics), 60(2):301-316.

  M.J. Cox, D.L. Borchers and N. Kelly. 2013. nupoint: An R
  package for density estimation from point transects in
  the presence of non-uniform animal density Methods in
  Ecology and Evolution 4(6):589-594

  Marques, T.A. , Buckland, S.T. , Borchers, D.L. , Tosh,
  D. and McDonald, R.A. 2010.  Point transect sampling
  along linear features.  Biometrics 66(4):1247-1255.
}

