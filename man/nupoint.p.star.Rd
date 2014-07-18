\name{nupoint.p.star}
\alias{nupoint.p.star}
\title{Calculate the probability of detecting a target within the covered area given
a parallel density gradient.}
\usage{
  nupoint.p.star(pars.hat, grad.type, det.type, n = NULL,
    w, theta.max, grid.density = 100, attenuation = FALSE,
    angularDetect = FALSE)
}
\arguments{
  \item{pars.hat}{\code{a single vector of parameters used
  to estimate the probablity of detecting target within
  covered area.  Parameters must be ordered:environment
  density gradient first, followed by detection function.
  These parameters will typically be estimated using
  nupoint.fit}}

  \item{grad.type}{\code{parametric form of the target
  object density function. See nupoint.fit and
  nupoint.env.fit.}}

  \item{det.type}{\code{parametric form of the detection
  function.  See detectF.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the depth preference function.}}

  \item{w}{\code{truncation distance in the y-dimension
  (depth for a multi-beam echosounder)}}

  \item{theta.max}{\code{maximum observation angle radians,
  (multi-beam echosounder swath width).  Maximum angle is
  pi/2.}}

  \item{grid.density}{\code{numerical integration grid
  density (number of elements in each (x and y)
  dimensions).}}

  \item{attenuation}{\code{logical:default FALSE.  If TRUE,
  nupoint.p.star uses the attenuation function, atten.f, in
  the R workspace.  See make.atten.f.}}

  \item{angularDetect}{\code{logical:default FALSE. Should
  angular detection be considered (not yet implemented).}}
}
\value{
  Probability of detecting a target in the covered region.
}
\description{
  This function calculates the probability of detecting a
  target, e.g. krill swarm, within the sampled area.
  Parameters estimated using \code{\link{nupoint.fit}} are
  typically used in the \code{pars.hat} argument
}
\details{
  This function calculates the probability of detecting a
  target (e.g. krill swarm, fish shoal, or single fish) in
  the area covered by a multi-beam echosounder. Typically,
  this function will calculate the probability of detecting
  a target using the parameter estimates calculated in
  \code{\link{nupoint.fit}}.
}
\examples{
\dontrun{
nupoint.p.star(pars.hat=c(50,20,60),grad.type='NORM',det.type='HNORM',n=NULL,w=100,
 theta.max=pi/3,grid.density=100,
 attenuation=FALSE,angularDetect=FALSE)
}
}
\seealso{
  \code{\link{nupoint.fit}}
}
\keyword{misc}

