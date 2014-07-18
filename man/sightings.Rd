\docType{data}
\name{sightings}
\alias{sightings}
\title{Sightings of whale blows from shore vantage point}
\description{
  Blows seen from land by Arranz and colleagues.
}
\details{
  \itemize{ Sightings data consists of a list composed of
  14 elements.  This is information that can be used as
  input to \code{nupoint.env.fit}.  Consult example
  associated with \code{nupoint.env.fit}. Elements of the
  sightings data list: \item sightings.mat dataframe of
  x,y,r,z,dzdy and obs.period of sightings \item x.mat
  x-coordinates of grid \item y.mat y-coordinates of grid
  \item rd.mat radial distances of grid points \item z.mat
  covariate value of grid points \item zGradmat derivative
  of coordinate in y dimension for all grid points \item x
  cutpoints in x-dimension of grid points \item y cutpoints
  in y-dimension of grid points \item obsx observer
  position in x \item obsy observer position in y \item wx
  maximum distance of study area in x-dimension \item wy
  maximum distance of study area in y-dimension \item wz
  maximum value of measured covariate (z) \item minz
  minimum value of measured covariate (z) }
}
\note{
  The beaked whale data presented in this package were
  collected by the University of La Laguna (ULL) with
  funding from 2003 to 2010, of the Woods Hole
  Oceanographic Institution (WHOI). In 2004, the studies
  were co-funded by the Department of the Environment of
  the Canary Islands, through an agreement made with the
  Ministry of Defence and in 2010, by the Department of the
  Environment of the Government of Spain, through the
  national project CETOBAPH. These data can not be used for
  analysis, databases or to be published in scientific
  journals, without the prior consent of the University of
  La Laguna (contact \email{snaguilar@ull.es},
  \email{arranz@ull.es}).
}
\references{
  P. Arranz, D.L. Borchers, N. Aguilar Soto, M.P. Johnson,
  and M.J. Cox. A new method to study inshore whale cue
  distribution from land-based observations. Marine Mammal
  Science, submitted.
}
\keyword{datasets}

