\name{nupoint.env.gof}
\alias{nupoint.env.gof}
\title{Chi-squared goodness-of-fit test for environmental data}
\usage{
  nupoint.env.gof(pars, r.mat, z.mat, minz, wz, z.obs,
    grad.type, det.type, intervals, plot = FALSE, n = NULL,
    dzdy.mat, breaks = NULL)
}
\arguments{
  \item{pars}{\code{vector of parameters estimated using
  nupoint.env.fit.}}

  \item{r.mat}{\code{matrix:radial distances from observer
  to each regular grid point in the survey region.}}

  \item{z.mat}{\code{matrix:regular grid of the
  environmental feature throughout the survey region.}}

  \item{minz}{\code{minimum environmental feature value.}}

  \item{wz}{\code{truncation distance in the z-dimension
  (environmental feature).}}

  \item{z.obs}{\code{observation vector:environmental
  feature.}}

  \item{grad.type}{\code{parametric form of the
  environmental gradient (this is the function describing
  habitat preference).  See nupoint.env.fit.}}

  \item{det.type}{\code{parametric form of the detection
  function.  See detectF.}}

  \item{intervals}{\code{Either an integer number of bins
  for chi-squared test (which will be used to generate bins
  of equal width).  Alternatively, bins of unequal can be
  specified using vector break points of length 1+ the
  number of desired bins.}}

  \item{plot}{\code{logical default=FALSE. If TRUE a four
  panel diagnostic plot is displayed.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the environmental gradient.  See
  nupoint.env.fit.}}

  \item{dzdy.mat}{\code{matrix:regular grid of the rate of
  change of environmental feature with respect to the
  y-dimension.}}

  \item{breaks}{\code{vector of break points to calculate
  chi-squared goodness-of-fit.  The argument
  nbr.intervals=NULL when breaks is used.}}
}
\value{
  List object of length two. 1) Chi-squared p-value, and 2)
  goodness-of-fit table with columns: bin minimum, bin
  maximum, expected sightings, observed sightings, Chisq
}
\description{
  This function calculates the one-dimensional chi-squared
  goodness-of-fit statistic for environment preference
  data.
}
\details{
  The \code{pars} argument is the vector of parameter
  estimates obtained from \code{\link{nupoint.env.fit}}.
  The arguments \code{grad.type}, \code{det.type} and
  \code{n} must ,match those used in
  \code{\link{nupoint.env.fit}} to obtain estimates of
  \code{pars}.
}
\examples{
##Don't run:
attach(sightings) # subset of beaked whale data from Arranz (submitted)
#fit a normal environmental feature preference and half-normal detection function:
environ.fit=nupoint.env.fit(pars=c(1000,200,3000),
             z=sighting.mat$z, rd=sighting.mat$r, dzdy=sighting.mat$dzdy,
       z.mat=z.mat, dzdy.mat=zGradmat, rd.mat=rd.mat,
             minz=minz, wx=wx, wy=wy, wz=wz,
             grad.type="NORM", det.type="HNORM",
             lower.b=c(-2000,1,1),upper.b=c(10000,10000,10000))
#after running nupoint.env.fit using the example sightings data with 'grad.type="NORM", det.type="HNORM"':
tt=nupoint.env.gof(pars=environ.fit$par, r.mat=rd.mat, z.mat=z.mat, minz=minz,wz=wz,
       z.obs=sighting.mat$z, grad.type="NORM", det.type="HNORM",
       intervals=13,plot=FALSE,dzdy.mat=zGradmat)
#-------------------------------------------------------
#1D Chi-squared Goodness-of-Fit results
#-------------------------------------------------------
#   bin.min bin.max    mids expected observed Chisq
#1   103.22  241.73  172.47     5.08        7  0.73
#2   241.73  380.25  310.99     8.81        6  0.89
#3   380.25  518.77  449.51    15.43       16  0.02
#4   518.77  657.29  588.03    22.63       17  1.40
#5   657.29  795.81  726.55    24.77       24  0.02
#6   795.81  934.32  865.07    28.88       29  0.00
#7   934.32 1072.84 1003.58    25.65       30  0.74
#8  1072.84 1211.36 1142.10    21.33       24  0.33
#9  1211.36 1349.88 1280.62    16.85       19  0.28
#10 1349.88 1488.40 1419.14    11.22        9  0.44
#11 1488.40 1626.91 1557.66     6.46        7  0.04
#12 1626.91 1765.43 1696.17     2.93        0  2.93
#13 1765.43 1903.95 1834.69     0.97        2  1.10
#-------------------------------------------------------
#Chi-squ. statistic   =  8.926069
#Number of parameters =  3
#Chi-squ. df          =  9
#Chi-squ. GoF p-value =  0.4441265
#-------------------------------------------------------
#
#example using user defined bin break points
binV=seq(sightings$minz,sightings$wz,length=14) #GoF bin break point vector.
binV=binV[-5] #remove break 5
tt=nupoint.env.gof(pars=environ.fit$par, r.mat=rd.mat, z.mat=z.mat, minz=minz,wz=wz,
       z.obs=sighting.mat$z, grad.type="NORM", det.type="HNORM",
       intervals=binV,plot=FALSE,dzdy.mat=zGradmat)
#-------------------------------------------------------
#1D Chi-squared Goodness-of-Fit results
#-------------------------------------------------------
#   bin.min bin.max    mids expected observed Chisq
#1   103.22  241.73  172.47     5.08        7  0.73
#2   241.73  380.25  310.99     8.81        6  0.89
#3   380.25  518.77  449.51    15.43       16  0.02
#4   518.77  795.81  657.29    47.40       41  0.86
#5   795.81  934.32  865.07    28.88       29  0.00
#6   934.32 1072.84 1003.58    25.65       30  0.74
#7  1072.84 1211.36 1142.10    21.33       24  0.33
#8  1211.36 1349.88 1280.62    16.85       19  0.28
#9  1349.88 1488.40 1419.14    11.22        9  0.44
#10 1488.40 1626.91 1557.66     6.46        7  0.04
#11 1626.91 1765.43 1696.17     2.93        0  2.93
#12 1765.43 1903.95 1834.69     0.97        2  1.10
#-------------------------------------------------------
#Chi-squ. statistic   =  8.366981
#Number of parameters =  3
#Chi-squ. df          =  8
#Chi-squ. GoF p-value =  0.3984676
#-------------------------------------------------------
detach(sightings)
##End don't run
}
\seealso{
  \code{\link{nupoint.env.fit}}
}
\keyword{misc}

