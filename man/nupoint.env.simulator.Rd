\name{nupoint.env.simulator}
\alias{nupoint.env.simulator}
\title{Simulator for shore-based distance observations of animals distributed along
an environmental gradient}
\usage{
  nupoint.env.simulator(pars = c(60, 10, 50), z.mat = NULL,
    xlim = c(0, 200), ylim = c(0, 100),
    grid.resolution = 1, grad.type = "NORM",
    det.type = "HNORM", observer.coords = c(100, 0),
    nbr.targets = 350,
    environment.simulator.control = list(c(X = 50, Y = 10, sd = 60), c(X = 90, Y = 0, sd = 30)),
    slope.control = NULL, mask.mat = NULL, mask.ang = 0,
    plot = TRUE, perp.lines = NULL, n = NULL)
}
\arguments{
  \item{pars}{\code{vector of parameters to simulate from
  (see nupoint.env.fit).}}

  \item{z.mat}{\code{matrix:either a regular grid of the
  environmental feature throughout the survey region, or
  NULL if bathymetry is simulated using settings in the
  environment.simulator.control argument (see details).}}

  \item{xlim}{\code{x-coordinate range.}}

  \item{ylim}{\code{y-coordinate range.}}

  \item{grid.resolution}{\code{resolution of simulated
  survey region.}}

  \item{grad.type}{\code{parametric form of the
  environmental gradient.  See nupoint.env.fit.}}

  \item{det.type}{\code{parametric form of the detection
  function c("HNORM","HAZARD"). See detectF.}}

  \item{observer.coordinates}{\code{x,y observer
  position.}}

  \item{nbr.targets}{\code{Number of targets in survey
  region. When Pr(detect)=1 this is the number of detected
  targets.}}

  \item{environment.simulator.control}{\code{list object
  comprised of vectors. Each vector describes a radial
  basis function which are used to simulate the environment
  when z.mat=NULL (see details).}}

  \item{slope.control}{\code{default NULL. Otherwise a two
  element vector comprising of slope and intercept
  specifying a linearly increasing seabed function.  See
  details.}}

  \item{mask.mat}{\code{default NULL.  Otherwise a user
  defined matrix of dimensions equal to z.mat specifying
  areas within the survey region to be excluded from
  analysis (e.g. simulated coastline, see details).}}

  \item{mask.ang}{\code{angle (rad) at observer, measured
  from baseline (y=0) contains unsuitable habit (e.g. land;
  see details).  Either specify mask.ang or mask.ang}}

  \item{plot}{\code{logical. Default TRUE.  Returns a panel
  plot of detection function and envrionmental gradient
  with all targets and seen targets shown on both plots.}}

  \item{perp.lines}{\code{Default = NULL (see details).  Or
  a list object describing perpendicular transects lines
  run from baseline, list(nbr.transects,min.stop,max.stop).
  The min.stop and max.stop settings are used to specify
  the minimum and maximum ltransect length in the
  y-dimension.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the environmental gradient.  See
  nupoint.env.fit.}}
}
\value{
  list object of simulator settings, survey region
  descriptors and observations (see details).
}
\description{
  Creates a simulated survey region and a set of
  observations collected by a shore-based observer when
  animals are distributed along an environmental gradient.
}
\details{
  The user can either pass in a matrix specifying the
  environment variable (e.g. seabed) using the \code{z.mat}
  argument or can create an environment variable using the
  \code{environment.simulator.control} argument. If a
  matrix is passed into the function via the \code{z.mat}
  argument, then its dimensions must equal those of the
  matrix (if any) used in \code{mask.mat}. Further the
  range of \code{xlim} divided by \code{grid.resolution}
  must equal the number of rows in the matrix (\code{z.mat}
  and \code{mask.mat}) arguments, and the range of
  \code{ylim} divided by \code{grid.resolution} must equal
  the number of columns in \code{z.mat} and
  \code{mask.mat}.

  If \code{z.mat} is unspecified, the user can create a
  \code{z.mat} within the \code{nupoint.env.simulator}
  function by specifying radial basis functions using the
  \code{environment.simulator.control} argument, which is a
  list object comprised of vectors. Each vector describes a
  radial basis function \code{c(x-coordinate, y-coordinate,
  standard deviation)}.

  If \code{slope.control} is unspecified, the default is
  for the simulated seabed to equal the y-coordinate at a
  given position. This means the seabed increases linearly
  with distance from shore (the x-axis). The user can
  change the default seabed by specifying an intercept
  parameter \code{slope.control[1]} and slope parameter
  \code{slope.control[2]} so that a given seabed depth,
  \code{z} is \code{z}=slope.control[1] + slope.control[2]
  x y-coordinate.

  A mask matrix may be specified by the user using the
  \code{mask.mat} argument to exclude unsuitable
  environment from the survey region (e.g. land in the case
  of beaked whales). This matrix must have identical
  dimensions to \code{z.mat}. Alternatively, areas within
  the survey region may be excluded using the
  \code{mask.ang} argument, which removes areas that are
  less than \code{mask.ang} degrees either side of the
  observer. \code{mask.ang} is measured from the x-axis.

  Perpendicular line transects can be simulated using the
  \code{perp.lines} argument. The default is
  \code{prep.lines=NULL}, meaning no perpendicular survey
  lines are simulated, otherwise an even number of integer
  lines are created and populated with simulated
  observations.

  \code{nupoint.env.simulator} returns a list comprised of
  the following: \code{[[1]]} a data frame of simulated
  sightings \code{\link{sightings}}; \code{[[2]]}
  \code{rd.mat} a matrix of distances from observer to each
  point in the simulated survey region; \code{[[3]]}
  \code{z.mat} a regular matrix of the environmental
  feature throughout the survey region; \code{[[4]]}
  \code{zGradmat} regular matrix of the environmental
  feature gradient in the y-dimension; \code{[[5]]} a
  sub-list of simulator settings (\code{pars}, \code{xlim},
  \code{ylim}, \code{grid.resolution}, \code{grad.type},
  \code{det.type}, \code{observer.coords},
  \code{nbr.targets} arguments), and \code{[[6:10]]}:
  simulated perpendicular transect information if
  \code{perp.lines!=NULL}.

  The data frame of simulated sightings is made up of
  sighting x-coordinate, \code{x}, sighting y-coordinate,
  \code{y}, sighting environmental feature value, \code{z},
  environmental feature gradient, \code{dzdy}, and
  detection distance \code{d}.

  Perpendicular transect information is simulated the
  returned list contains additional information:
  \code{[[6]]}: \code{perpendicular.y} matrix of
  y-coordinates, when each column contains the y-values
  along a given simulated perpendicular transect;
  \code{[[7]]}: \code{perpendicular.y} matrix of
  environmental feature values, when each column contains
  values along a given simulated perpendicular transect;
  \code{[[8]]}: \code{perpendicular.dzdy} matrix of
  environmental feature gradient in the y-dimension;
  \code{[[9]]}: \code{perpendicular.transect.info}
  perpendicular transect information giving the start and
  end x-coordinates \code{startx} and \code{stopx}, and the
  start and end y-coordinates \code{starty} and
  \code{stopy} for each simulated transect, and
  \code{[[10]]} \code{perpendicular.transect.data}:
  observations collected on each simulated transect
  (identical format to element \code{[[1]]}
  \code{\link{sightings}} plus a column named
  \code{transect} that assigns an observation to a
  perpendicular transect.
}
\examples{
environ.sim.dat=nupoint.env.simulator(pars=c(60,10,50),z.mat=NULL,xlim=c(0,200),ylim=c(0,100),
                      grid.resolution=1,grad.type='NORM',det.type='HNORM',
                      observer.coords=c(100,0),nbr.targets=350,
                      environment.simulator.control=list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
                      mask.mat=NULL,mask.ang=0,plot=TRUE,
                      perp.lines=NULL,n=NULL)
}
\seealso{
  \code{\link{detectF}}, \code{\link{nupoint.env.fit}}
}
\keyword{misc}

