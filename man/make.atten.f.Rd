\name{make.atten.f}
\alias{make.atten.f}
\title{Fit a attenuation function}
\usage{
  make.atten.f(prep.dist, form, starting.pars,
    plot = FALSE, summary = FALSE)
}
\arguments{
  \item{prep.dist}{\code{vector of a representative sample
  of prependicular distances(e.g. in the multibeam context,
  seabed depth under a krill swarm).}}

  \item{form}{\code{any standard R distribution (see
  Distributions) that can be acommodated by the FAmle
  package e.g. norm for normal.}}

  \item{starting.parameters}{\code{vector of starting
  parameters describing the distribution specified in
  form.}}

  \item{plot}{\code{logical. produce a plot of the
  attenuation function}}

  \item{summary}{\code{default=FALSE. If TRUE, print a
  summary of the attenuation function fit to the console.}}
}
\value{
  1) creates a function in the R workspace called
  \code{atten.f} 2) returns a \code{FAmle::mle} object.
}
\description{
  Fit an attenuation function calls the \code{mle} function
  from the \pkg{FAmle} package. The seabed attenuation
  function accounts for variation in maximum transect
  observation distance. In terms of multi-beam echosounder
  surveys, the attenuation function accounts for varying
  seabed depth. See Cox et al. (2011).
}
\details{
  This function must be run prior to
  \code{\link{nupoint.fit}} if an attenuation function is
  required in the likelihood fitted in
  \code{\link{nupoint.fit}}. The name of the \code{atten.f}
  must not be changed.
}
\section{Reference}{
  Cox, M.J., Borchers, D.L., Demer, D.A., Cutter, G.R., and
  Brierley, A.S. (2011). Estimating the density of
  Antarctic krill (Euphausia superba) from multi-beam
  echo-sounder observations using distance sampling
  methods. Journal of the Royal Statistical Society: Series
  C (Applied Statistics).
}
\examples{
#use the seabed depth under each krill swarm 'krill$z' to create an attenuation function based on a normal distribution
\dontrun{
ls()
atten.fit=make.atten.f(prep.dist=krill$z,form='norm',starting.pars=c(90,10))
ls() #function named atten.f now exists in workspace
}
}
\seealso{
  \code{\link{nupoint.fit}}
}
\keyword{misc}

