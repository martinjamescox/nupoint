\name{nupoint.gof}
\alias{nupoint.gof}
\title{Chi-squared goodness-of-fit test for parallel density gradients}
\usage{
  nupoint.gof(y.obs, pars, w, grad.type, det.type,
    n = NULL, intervals, verbose = TRUE, plot = FALSE,
    attenuation = FALSE, angularDetect = FALSE,
    xmax = NULL, theta.max, grid.resolution = 100)
}
\arguments{
  \item{y.obs}{\code{observation vector: target y
  coordinate. This is water depth at the centre of an
  aggregation for multi-beam echosounder observations.}}

  \item{pars}{\code{vector of parameters estimated using
  nupoint.fit.}}

  \item{w}{\code{truncation distance in the y-dimension
  (depth for multi-beam echosounder observations)}}

  \item{grad.type}{\code{parametric form of the target
  object density function. See nupoint.fit and
  nupoint.env.fit.}}

  \item{det.type}{\code{parametric form of the detection
  function.  See detectF.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the environmental gradient.  See
  nupoint.fit.}}

  \item{intervals}{\code{Either an integer number of bins
  for chi-squared test (which will be used to generate bins
  of equal width).  Alternatively, bins of unequal can be
  specified using vector break points of length 1+ the
  number of desired bins.}}

  \item{verbose}{\code{return summary text of
  goodness-of-fit to the R console.}}

  \item{plot}{\code{logical default=FALSE. If TRUE a two
  panel diagnostic plot is displayed.}}

  \item{attenuation}{\code{Logical.  TRUE = seabed
  attenuation function is specified; FALSE = no seabed
  attenuation function.}}

  \item{angularDetect}{\code{logical default=FALSE. Should
  angular detection be considered (not yet implemented).}}

  \item{xmax}{\code{maximum x-dimension distance. Default
  is NULL, xmax = w * sin(theta).}}

  \item{theta.max}{\code{maximum observation angle radians,
  (swath width in the multi-beam echosounder case).
  Maximum angle is pi/2.}}

  \item{grid.resolution}{\code{numerical integration grid
  density (number of elements in each (x and y)
  dimensions).}}
}
\value{
  List object of length two. 1) Chi-squared p-value, and 2)
  goodness-of-fit table with columns: bin minimum, bin
  maximum, expected sightings, observed sightings, Chisq.
}
\description{
  This function calculates the one-dimensional chi-squared
  goodness-of-fit statistic for environment preference
  data.
}
\details{
  The \code{pars} argument is the vector of parameter
  estimates obtained from \code{\link{nupoint.fit}}. The
  arguments \code{grad.type}, \code{det.type} and \code{n}
  must ,match those used in \code{\link{nupoint.fit}} to
  obtain estimates of \code{pars}.
}
\examples{
\dontrun{
# norm.gof=nupoint.gof(y.obs=krill$y,
#                pars=c(81.736 ,25.144 ,38.194),
#                w=100,
#                grad.type='NORM',
#                det.type='HNORM',
#                n=NULL,
#                intervals=11,
#                angularDetect=FALSE,
#                theta.max=pi/3,
#                grid.resolution=100)
# -------------------------------------------------------
#   nupoint: 1D Chi-squared Goodness-of-Fit results
# -------------------------------------------------------
#   bin.min bin.max  mids expected observed Chisq
# 1    10.74   18.69 14.71     0.95        2  1.15
# 2    18.69   26.63 22.66     3.96        4  0.00
# 3    26.63   34.58 30.60    10.67        9  0.26
# 4    34.58   42.52 38.55    21.39       26  0.99
# 5    42.52   50.47 46.50    28.88       28  0.03
# 6    50.47   58.42 54.44    41.75       34  1.44
# 7    58.42   66.36 62.39    44.67       41  0.30
# 8    66.36   74.31 70.34    40.51       46  0.74
# 9    74.31   82.26 78.28    30.78       41  3.39
# 10   82.26   90.20 86.23    19.08       16  0.50
# 11   90.20   98.15 94.17     8.35        4  2.27
# -------------------------------------------------------
#   Chi-squ. statistic   =  11.06844
# Number of parameters =  3
# Chi-squ. df          =  7
# Chi-squ. GoF p-value =  0.1356615
# -------------------------------------------------------
}
}
\seealso{
  \code{\link{nupoint.fit}}, \code{\link{detectF}},
  \code{\link{nupoint.env.fit}}
}
\keyword{misc}

