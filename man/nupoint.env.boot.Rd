\name{nupoint.env.boot}
\alias{nupoint.env.boot}
\title{Non-parametric bootstrap for environmental preference}
\usage{
  nupoint.env.boot(sightings, nboot, blockVar = NULL,
    initial.pars, grad.type, det.type, n = NULL, lower.b,
    upper.b)
}
\arguments{
  \item{sightings}{\code{list object of sightings and
  survey region variables (see details and the }env.pdf'
  vignette).'}

  \item{nboot}{\code{number of bootstraps.}}

  \item{blockVar}{\code{string column name or integer
  column number in the sightings$sighting.mat data frame
  which is the bootstrap block variable.}}

  \item{blockVar}{\code{string column name or integer
  column number in the sightings$sighting.mat data frame
  which is the bootstrap block variable.}}

  \item{initial.pars}{\code{intitial parameters for use by
  the nupoint.env.fit function (see the pars argument in
  nupoint.env.fit).}}

  \item{grad.type}{\code{environmental gradient type (See
  nupoint.env.fit).}}

  \item{det.type}{\code{detection function form. (See
  detecF and nupoint.env.fit).}}

  \item{n}{\code{number of distributions in a mixture. (See
  nupoint.env.fit).}}

  \item{lower.b}{\code{lower parameter space bound for
  nupoint.env.fit function.}}

  \item{upper.b}{\code{upper parameter space bound for
  nupoint.env.fit function.}}
}
\value{
  2D array nrow = nboot; ncol = initial.pars, each row
  contains parameter estimates for each bootstrap.
}
\description{
  Calculate the variance of parameter estimates for the
  environmental preference gradient using a non-parameteric
  bootstrap.
}
\details{
  The \code{\link{sightings}} argument is a list object
  which must have the following structure:

  \describe{ \item{list("sighting.mat")}{whale cue
  sightings;} \item{list("x.mat")}{survey region
  description: x-coordinates;} \item{list("y.mat")}{survey
  region description: y-coordinates;}
  \item{list("rd.mat")}{survey region description: radial
  distance;} \item{list("z.mat")}{survey region
  description: environment variable;}
  \item{list("zGradmat")}{survey region description: change
  in environment variable;} \item{list("x")}{survey region
  description: x-coordinate vector;}
  \item{list("y")}{survey region description: y-coordinate
  vector;} \item{list("obsx")}{observer x-coordinate;}
  \item{list("obsy")}{observer y-coordinate;}
  \item{list("wx")}{x-dimension truncation distance;}
  \item{list("wy")}{y-dimension truncation distance;}
  \item{list("wz")}{z-dimension truncation distance;}
  \item{list("minz")}{Minimum z dimension value;} }

  The \code{sighting.mat} object within
  \code{\link{sightings}} contains the following data for
  each observation

  \describe{ \item{list("x")}{x-coordinate of sighting;}
  \item{list("y")}{y-coordinate of sighting;}
  \item{list("r")}{radial distance from observer to
  sighting;} \item{list("z")}{seabed depth (or other
  environmental variable of interest) at each sighting;}
  \item{list("dzdy")}{seabed depth (or other environmental
  variable of interest) gradient with respect to the
  y-dimension at each sighting, and}
  \item{list("obs.period")}{observation period. This
  variable is used as the grouping parameter in
  \code{nupoint.env.boot} and is specified in the
  \code{blockVar} argument.} }

  Example sighting data \code{\link{sightings}} is included
  with the package.
}
\examples{
\dontrun{
boot.res=nupoint.env.boot(sightings=sightings,
              nboot=100,
              blockVar='obs.period',
              initial.pars=c(1171.479,503.87,3020.442),
              grad.type='NORM',
              det.type='HNORM',
              lower.b=c(-2000,1,1),
              upper.b=c(10000,10000,10000))
par(mfrow=c(1,3))
sapply(colnames(boot.res),function(x) hist(boot.res[,x],main=x,xlab='Estimate'))
}
}
\seealso{
  \code{\link{detectF}}, \code{\link{nupoint.env.fit}},
  \code{\link{sightings}}
}
\keyword{misc}

