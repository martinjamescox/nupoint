\docType{package}
\name{nupoint}
\alias{nupoint}
\alias{nupoint-package}
\title{'nupoint': An R package for density estimation from point transects in the
presence of non-uniform animal density}
\description{
  The R package \pkg{nupoint} provides tools for estimating
  animal density, using distance sampling methods, from
  point transect surveys in which the conventional point
  transect assumption of uniform animal distribution in the
  vicinity of the point is violated.
}
\details{
  The key features of the package are as follows:

  1. By using detection angles as well as distances it
  allows estimation of both a conventional point transect
  radial distance detection function and (unlike
  conventional methods) non-uniform animal density with
  respect to some observable environmental feature within
  the searched area (such as distance from road or from
  shore, or altitude, or depth).

  2. Key outputs are estimates of density itself, density
  function parameters, detection function parameters and
  associated coefficients of variation and confidence
  intervals.

  3. The package also provides some tools for plotting and
  for model diagnosis, AIC statistics for model selection
  and goodness-of-fit test statistics.

  4. Finally, a facility to investigate estimator
  properties by simulation is also provided.

  The \pkg{nupoint} package handles two types of
  non-uniform density: (1) when density contours are
  parallel, these functions are named as \code{nupoint.*},
  and (2) when density varies with some environmental
  variable such as seabed depth, these functions are named
  \code{nupoint.env.*}.

  The package includes tools for fitting models
  \code{\link{nupoint.fit}} and
  \code{\link{nupoint.env.fit}}, model selection,
  goodness-of-fit testing \code{\link{nupoint.gof}} and
  \code{\link{nupoint.env.gof}}. Variance estimation is by
  non-parametric bootstrap \code{\link{nupoint.boot}} and
  \code{\link{nupoint.env.boot}}.

  A survey region and sightings simulator
  \code{\link{nupoint.env.simulator}} is also included.

  Two datasets: \code{\link{krill}} as an example of
  parallel density contours, and \code{\link{sightings}}
  for varying density with an environmental variable, are
  included and have accompanying vignettes,
  \code{parallel.pdf} and \code{env.pdf}.
}
\section{Authors}{
  Martin j Cox, David l Borchers and Natalie Kelly
}

\section{Sightings dataset}{
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
  \email{arranz@ull.es}). Multi-beam instrument was loaned
  by J. Condiotty of Simrad USA. Multi-beam data were
  collected in association with an NSF-funded (grant
  #06-OPP-33939) investigation of the Livingston Island
  nearshore environment. Logistical support was provided by
  US Antarctic Marine Living Resources Program, and
  engineering design and support by Sea Technology
  Services. The multi-beam instrument was deployed using
  funds provided by the Royal Society. We are grateful to
  Francois Aucoin for allowing us to use the functions
  \code{mle} and \code{distr} from the \pkg{FAmle} package.
  MJC is funded by Australian Research Council grant
  FS110200057.
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
  the presence of non-uniform animal density. Methods in
  Ecology and Evolution 4(6):589-594. DOI:
  10.1111/2041-210X.12058
}
\keyword{misc}

