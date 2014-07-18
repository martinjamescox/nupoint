\name{nupoint.fit}
\alias{nupoint.fit}
\title{Maximum likelihood estimate of parameters describing parallel density
gradient and range detection functions}
\usage{
  nupoint.fit(pars, sight.x, sight.y, sight.z = NULL, w,
    theta.max, grad.type, det.type, n = NULL,
    grid.density = 100, angularDetect = FALSE,
    verbose = FALSE, lower.b, upper.b,
    optim.control = NULL)
}
\arguments{
  \item{pars}{\code{a single vector of starting parameters
  to be estimated. Parameters must be ordered:environment
  density gradient first, followed by detection function.
  See nupoint.env.fit.}}

  \item{sight.x}{\code{observation vector: target x
  coordinate. This is cross-track distance for multi-beam
  echosounder(MBE) observations.}}

  \item{sight.y}{\code{observation vector: target y
  coordinate. This is water depth at the centre of an
  aggregation for MBE observations.}}

  \item{sight.z}{\code{observation vector: default NULL if
  no attenuation function is required (all transects extend
  to w). For MBE observations, this is water depth to the
  centre of a swarm.}}

  \item{w}{\code{truncation distance in the y-dimension
  (depth for MBE)}}

  \item{theta.max}{\code{maximum observation angle radians,
  (MBE swath width). Maximum angle is pi/2.}}

  \item{grad.type}{\code{parametric form of the
  environmental gradient.  See nupoint.env.fit.}}

  \item{det.type}{\code{parametric form of the detection
  function.  See detectF.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the environmental gradient.  See
  details.}}

  \item{grid.density}{\code{numerical integration grid
  density (number of elements in each (x and y)
  dimensions).}}

  \item{angularDetect}{\code{logical:FALSE. Should angular
  detection be considered (not yet implemented).}}

  \item{verbose}{\code{logical. Default FALSE. If TRUE
  parameter estimates and log-likelihood are printed to the
  console for each optim search iteration.}}

  \item{lower.b}{\code{vector of lower parameter space
  bounds.  Length must equal that of the pars argument. }}

  \item{upper.b}{\code{vector of upper parameter space
  bounds.  Length must equal that of the pars argument. }}

  \item{optim.control}{\code{control list object to be
  passed to the optim function.  See optim help and
  nupoint.env.fit.}}
}
\value{
  output as per optim and AIC.
}
\description{
  This function calculates maximum likelihood estimates of
  parameters for both a parallel density gradient and
  detection function, for which there are a variety of
  options for parametric forms.
}
\details{
  This likelihood function is for use when determining a
  density gradient with respect to a parallel feature e.g.
  road (Marques et al., 2010 and sea surface (Cox et al.,
  2011).
}
\examples{
##Don't run:
#Example 1
norm.fit <-   nupoint.fit(c(50,20,50),
                      sight.x=krill$x,
                      sight.y=krill$y,
                      w=100,
                      theta.max=pi/3,
                      grad.type="NORM",
                      det.type="HNORM",
                      lower.b=c(1,1,1),
                      upper.b=c(120,100,100))
# ---------------------------
#   Parallel density gradient likelihood settings
# ---------------------------
#   Depth preference parametric form: NORM
# range detection function, g(r), parametric form: HNORM
# Parameter starting values = 50 20 50
# Angular detection function applied: FALSE
# Seabed attenuation function applied: FALSE
# Truncation distance = 100
# Swath (max) angle, rad = 1.047198
# Numerical integration grid density= 100
# ---------------------------
#   Estimating parameters
# ---------------------------
#   Maximum likelihood results
# ---------------------------
#   parameter estimates
# 81.736 25.144 38.194
# AIC = -3.02
# ---------------------------
#Example 2
#create an attenuation function:
at.fit=make.atten.f(depths=krill$z,form='norm',starting.pars=c(90,10),plot=TRUE)
#fit a Beta vertical distribution:
mbe.beta=nupoint.fit(pars=c(3,1.5,40), sight.x=krill$x, sight.y=krill$y,sight.z=krill$z,
                     w=100,theta.max=pi/3,grad.type='BETA',det.type='HNORM',grid.density=100,lower.b=c(0.5,0.5,10),upper.b=c(20,20,1000))
#---------------------------
#Parallel density gradient likelihood settings
#---------------------------
#Depth preference parametric form: BETA
#range detection function, g(r), parametric form: HNORM
#Parameter starting values = 3 1.5 40
#Angular detection function applied: FALSE
#Seabed attenuation function applied: TRUE
#function(seabedDepth) 1-pnorm(seabedDepth,95.1141535839473,10.5625798419803)
#<environment: 0x06e6beb0>
#Truncation distance = 100
#Swath (max) angle, rad = 1.047198
#Numerical integration grid density= 100
#---------------------------
#Estimating parameters
#---------------------------
#Maximum likelihood results
#---------------------------
#parameter estimates
#3.52 1.281 38.914
#AIC = 484.63
#---------------------------
##End don't run
}
\references{
  Cox, M.J., Borchers, D.L., Demer, D.A., Cutter, G.R., and
  Brierley, A.S. (2011). Estimating the density of
  Antarctic krill (Euphausia superba) from multi-beam
  echo-sounder observations using distance sampling
  methods. Journal of the Royal Statistical Society: Series
  C (Applied Statistics), vol 60, part 2, pp. 301–316.

  Marques, T.A. , Buckland, S.T. , Borchers, D.L. , Tosh,
  D. and McDonald, R.A. (2010). \code{ Point Transect
  Sampling Along Linear Features } Biometrics , vol 66, no.
  4, pp. 1247-1255.
}
\seealso{
  \code{\link{nupoint.env.fit}}, \code{\link{detectF}},
  \code{optim}
}
\keyword{misc}

