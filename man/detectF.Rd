\name{detectF}
\alias{detectF}
\title{Probability of detection at a given range}
\usage{
  detectF(rd, det.type, det.par)
}
\arguments{
  \item{rd}{\code{radial distance from observer to detected
  target.}}

  \item{det.type}{\code{parametric form of the detection
  function c("HNORM","HAZARD"). See details.}}

  \item{det.par}{\code{vector of detection function
  parameters.}}
}
\value{
  Detection probability at range \code{rd}.
}
\description{
  This function calculates the probability of detecting a
  target (e.g. krill swarm or beaked whale) at a given
  distance from the observer.
}
\details{
  The probability of detecting a target at a given radial
  distance from an observer can be modelled using the
  following parametric forms (\code{det.type} argument) and
  parameters passed to the function using the
  \code{det.par} argument:

  1) half-normal "HNORM", \code{det.par=c(sigma)}

  2) Hazard rate "HAZARD", \code{det.par=c(sigma,b)}

  NB If constructing an intital value vector, the
  \code{pars} argument in the \code{\link{nupoint.env.fit}}
  function, \code{pars=c(...,detect function parameters)}
  where \code{...} denotes environmental gradient
  parameters.
}
\section{Reference}{
  Buckland, S.T., Anderson,D. R., Burnham, K. P., Laake, J.
  L., Borchers,D. L. and Thomas, L. (2001) Introduction to
  Distance Sampling. Oxford: Oxford University Press, page
  47.
}
\examples{
rVec=0:100
pVec=detectF(rd=rVec, det.type='HNORM', det.par=50)
plot(rVec,pVec,type='l',xlab='Distance',ylab='Pr(detect)')
}
\seealso{
  \code{\link{nupoint.env.fit}}, \code{\link{nupoint.fit}}
}
\keyword{misc}

