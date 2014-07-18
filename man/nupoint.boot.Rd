\name{nupoint.boot}
\alias{nupoint.boot}
\title{Non-parametric bootstrap for depth preference}
\usage{
  nupoint.boot(observations, nboot, blockVar, initial.pars,
    grad.type, det.type, w, theta.max, n = NULL,
    attenuation = FALSE, grid.density = 100, lower.b,
    upper.b, angularDetect = FALSE, L = NULL, A = NULL,
    shoalBioMass = NULL)
}
\arguments{
  \item{observations}{\code{data.frame of shaol
  observations.  See details and }summary(krill)' object).}

  \item{nboot}{\code{number of bootstraps.}}

  \item{blockVar}{\code{string column name or integer
  column number in the observations data frame which is the
  bootstrap block variable.}}

  \item{initial.pars}{\code{intitial parameters for use by
  the nupoint.fit function (see the pars argument in
  nupoint.fit).}}

  \item{grad.type}{\code{parametric form of the
  environmental gradient.  See environ.fit.f.}}

  \item{det.type}{\code{parametric form of the detection
  function.  See detectF.}}

  \item{w}{\code{truncation distance in the y-dimension
  (depth for MBE)}}

  \item{theta.max}{\code{maximum observation angle rad,
  (MBE swath width). Maximum angle is pi/2.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the environmental gradient.  See
  details.}}

  \item{attenuation}{\code{Logical.  TRUE = seabed
  attenuation function is specified; FALSE = no seabed
  attenuation function.}}

  \item{grid.density}{\code{numerical integration grid
  density (number of elements in each (x and y)
  dimensions.}}

  \item{lower.b}{\code{lower parameter space bound for
  nupoint.fit function.}}

  \item{upper.b}{\code{upper parameter space bound for
  nupoint.fit function.}}

  \item{angularDetect}{\code{logical:FALSE. Should angular
  detection be considered (not yet implemented).}}

  \item{L}{\code{Total line transect length}}

  \item{A}{\code{Survey region area}}

  \item{shoalBioMass}{\code{either a numeric value giving
  the mean shoal biomass or a character string giving the
  regression equation, or NULL.  If NULL then survery
  region biomass is not calculated.}}
}
\value{
  2D array nrow = nboot; ncol = initial.pars, P*. Each row
  in the array contains parameter estimates for each
  bootstrap. Plus, if argument L is specified shoal
  volumetric and areal density. Plus, if argument A is
  specified, number of shoals, N. Plus, if shoal biomass
  (\code{shoalBioMass} argument) is specified, expected
  shoal biomass, E[b], and survey region biomass. See
  details and the \code{parallel.pdf} vignette.
}
\description{
  Calculate the variance of parameter and density estimates
  for a parallel density gradient using a non-parameteric
  bootstrap
}
\details{
  The \code{observations} argument is a data.frame object
  which must have the following structure:

  \describe{ \item{list("transect")}{Transect number}
  \item{list("x")}{Shoal cross-track distance}
  \item{list("y")}{Shoal depth} \item{list("r")}{Shoal
  radial distance} \item{list("theta")}{Shoal angle, rad}
  \item{list("bio.g")}{Shoal biomass}
  \item{list("z")}{Seabed depth under shoal} }

  For each bootstrap, the vertical distrubtion and
  detection function parameters are estimated along with P*
  (see \code{\link{nupoint.p.star.f}}). If total line
  transect length, L, is specified then shoal volumetric
  density, \code{shoal.vol.den.hat}, and shoal areal
  density \code{shoal.areal.den.hat} are also estimated. If
  survey region area, A, is specified then shoal abundance,
  \code{N.hat} is estimated. If the shoalBioMass argument
  is specified, then shoal expected biomass, \code{E.bio}
  is estimated as is survey region biomass, \code{E.bio
  area.biomass}.

  The \code{shoalBioMass} argument can be given as: NULL,
  is which case no shoal or area biomass is estimated; a
  fixed mean shoal biomass value (the same value is used in
  each bootstrap), or a linear model formula which will be
  calculated for each bootstrap e.g. \code{log(bio.g) ~
  theta}, see example.
}
\examples{
\dontrun{
boot.res <- nupoint.boot(observations=krill,
 nboot=50,
 blockVar='transect',
 initial.pars=c(81.736, 25.144, 38.194),
 grad.type='NORM',
 det.type='HNORM',
 w=100,
 theta.max=pi/3,
 n=NULL,
 attenuation=FALSE,
 grid.density=100,
 lower.b=c(1,1,1),
 upper.b=c(120,100,100),
 angularDetect=FALSE,
 L=11*2.5*1e3,
 A=3300000,
 shoalBioMass='log(bio.g) ~ theta')
}
}
\seealso{
  \code{\link{detectF}},
  \code{\link{nupoint.fit}},'nupoint.p.star',
  \code{\link{nupoint.vol.density}},'krill'
}
\keyword{misc}

