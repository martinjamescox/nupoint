\name{nupoint.env.fit}
\alias{nupoint.env.fit}
\title{Maximum likelihood estimate of parameters describing environmental density
gradient and range detection functions}
\usage{
  nupoint.env.fit(pars, z, rd, dzdy, z.mat, dzdy.mat,
    rd.mat, minz, wx, wy, wz, grad.type, det.type,
    verbose = FALSE, n = NULL, diagnostic = FALSE,
    print.summary = TRUE, lower.b, upper.b,
    optim.control = NULL)
}
\arguments{
  \item{pars}{\code{a single vector of starting parameters
  to be estimated. Parameters must be ordered:environment
  density gradient first, followed by detection function.
  See details.}}

  \item{z}{\code{observation vector: environmental
  feature.}}

  \item{rd}{\code{observation vector: radial distance from
  observer to sighting.}}

  \item{dzdy}{\code{observation vector: rate of change of
  environmental feature with respect to the y-dimension.}}

  \item{z.mat}{\code{matrix:regular grid of the
  environmental feature throughout the survey region.}}

  \item{dzdy.mat}{\code{matrix:regular grid of the rate of
  change of environmental feature with respect to the
  y-dimension.}}

  \item{rd.mat}{\code{matrix:radial distances from observer
  to each regular grid point in the survey region.}}

  \item{minz}{\code{minimum environmental feature value.}}

  \item{wx}{\code{truncation distance in the x-dimension.}}

  \item{wy}{\code{truncation distance in the y-dimension.}}

  \item{wz}{\code{truncation distance in the z-dimension
  (environmental feature).}}

  \item{grad.type}{\code{parametric form of the
  environmental gradient.  See details.}}

  \item{det.type}{\code{parametric form of the detection
  function.  See details and detectF.}}

  \item{verbose}{\code{Logical. Default FALSE. If TRUE
  parameter estimates and log-likelihood are printed to the
  console for each optim search iteration.}}

  \item{n}{\code{Default NULL. Use only when
  grad.type=MNORM to specify an integer number of
  distributions in the normal mixture option for the
  parametric form of the environmental gradient.  See
  details.}}

  \item{diagnostic}{\code{Logical. Default FALSE. If TRUE
  diagnostic information for each call to the likelihood
  numerator and denominator functions is printed to the
  console.}}

  \item{lower.b}{\code{Vector of lower parameter space
  bounds.  Length must equal that of the pars argument. }}

  \item{upper.b}{\code{Vector of upper parameter space
  bounds.  Length must equal that of the pars argument. }}

  \item{optim.control}{\code{control list object to be
  passed to the optim function.  See optim help and
  details.}}
}
\value{
  output as per \code{optim} and AIC.
}
\description{
  This function calculates maximum likelihood estimates of
  parameters for both an environment density gradient and
  detection function, for which there are a variety of
  options for parametric forms.
}
\details{
  The length of the parameter starting value vector
  (\code{pars} argument) is dependent on the parametric
  forms of the environmental gradient (\code{grad.type}
  argument) and detection function (\code{det.type}
  argument). Parameters for the environmental gradient are
  placed first in the \code{pars} argument, followed by the
  detection function parameter(s).

  The following environmental gradient parametric forms are
  available of which the following have two parameters
  (plus detection function parameters):

  Parametric form name, \code{grad.type=}, parameter vector
  to be passed into the \code{pars} argument)

  Normal, "NORM" , \code{pars=c(mu,sd,....)}

  Log-normal, "LNORM", \code{pars=c(log mu, log sd, ....)}

  Beta, "BETA", \code{pars=c(shape_1,shape_2,....)}

  Uniform, "UNIFORM", \code{pars=c(....)}

  NB "...." denotes starting parameter values for the
  detection function (see detectF).

  The rate of change of the environmental feature with
  respect to the y-dimension for each observation
  \code{dzdy}, and the survey region \code{dzdy.mat} must
  be calculated by the user.

  Matrices describing the survey region (arguments
  \code{z.mat} , \code{dzdy.mat} , \code{rd.mat}) must have
  identical dimensions.

  Normal-mixture \code{grad.type=} "MNORM" may also be
  specified with the number of component normal
  distributions specified using the \code{n} argument. The
  starting parameter vector has length \code{(3*n-1)} plus
  detection function parameters and has the form:

  \code{pars =
  c(mu_1,sd_1,alpha_1,mu_2,sd_2,alpha_2,.....mu_n,sd_n,....)}

  where "..." are detection function parameters (see
  \code{\link{detectF}}).

  NB the weight, \code{alpha}, for the final distribution
  is not estimated so must not be included in \code{pars}.
  The final weight is \code{1-sum(alpha_1,
  ...,alpha_(n-1))}.

  Vectors of lower and upper parameter space bounds,
  \code{lower.b} and \code{upper.b} must both be the same
  length as pars and include the starting parameter value
  within their range.

  Options for the \code{optim} call within
  \code{nupoint.env.fit} may be controlled via the
  \code{optim.control} argument (see the \code{control}
  argument in \code{optim} for details). Because this is a
  maximisation problem, care must be taken not to overwrite
  the \code{fnscale} variable in \code{optim.control} with
  a positive number.
}
\examples{
attach(sightings) # subset of beaked whale data from Arranz (submitted)
#fit a normal environmental feature preference and half-normal detection function:
environ.fit=nupoint.env.fit(pars=c(1000,200,3000),
             z=sighting.mat$z, rd=sighting.mat$r, dzdy=sighting.mat$dzdy,
       z.mat=z.mat, dzdy.mat=zGradmat, rd.mat=rd.mat,
             minz=minz, wx=wx, wy=wy, wz=wz,
             grad.type="NORM", det.type="HNORM",
             lower.b=c(-2000,1,1),upper.b=c(10000,10000,10000))
\dontrun{
# ------------------------------------------------------
# Environmental gradient likelihood settings
# ------------------------------------------------------
# Environment preference parametric form: NORM
# range detection function, g(r), parametric form: HNORM
# Parameter starting values = 1000 200 3000
# Truncation distances
# x= 6450  ; y= 8850  ; z= 1903.951
# ------------------------------------------------------
# Estimating parameters
# ------------------------------------------------------
# Maximum likelihood results
# ------------------------------------------------------
# parameter point estimates = 1171.479 503.87 3020.442
# AIC = 3893.23
# ------------------------------------------------------
detach(sightings)
}
##example of a 2 component mixture normal:
attach(sightings)
det.type="HNORM"
sigma.r=3000 #half-normal sigma detection function parameter
grad.type="MNORM" #seabed depth cue distribution shape.
nDist=2 #number of normal distributions.
wt=rep(1/nDist,nDist) #distribution weights (final element removed later)
mu=seq(minz,wz,length=nDist) #distribute means along the environmental gradient
sigma=rep(300,nDist)
pars=as.vector(matrix(c(mu,sigma,wt),ncol=nDist,byrow=TRUE))
pars=c(pars[-length(pars)],sigma.r)
#bounds
mumin=rep(-2000,nDist)
sigmamin=rep(1,nDist)
alphamin=rep(-15,nDist)
sigma.rmin=1
lower.b=as.vector(matrix(c(mumin,sigmamin,alphamin),ncol=nDist,byrow=TRUE))
lower.b=c(lower.b[-length(lower.b)],sigma.rmin)
mumax=rep(5e+3,nDist)
sigmamax=rep(1e+8,nDist)
alphamax=rep(15,nDist)
sigma.rmax=1e+5
upper.b=as.vector(matrix(c(mumax,sigmamax,alphamax),ncol=nDist,byrow=TRUE))
upper.b=c(upper.b[-length(upper.b)],sigma.rmax)
#5) fit using optim with the  nupoint.env.fit function
environ.fit=nupoint.env.fit(pars=pars,
             z=sighting.mat$z,
                      rd=sighting.mat$r,
                       dzdy=sighting.mat$dzdy,
                      z.mat=z.mat,
                      dzdy.mat=zGradmat,
                      rd.mat=rd.mat,
                      minz=minz,
                      wx=wx,
                      wy=wy,
                      wz=wz,
                      grad.type=grad.type,
                      det.type=det.type,
                      n=nDist,lower.b=lower.b,upper.b=upper.b)
\dontrun{
#------------------------------------------------------
#Environmental gradient likelihood settings
#------------------------------------------------------
#Environment preference parametric form: MNORM
#range detection function, g(r), parametric form: HNORM
#Mixture of normal distributions with starting values =
#   mu sigma weight
# mixture-1  103.2158   300    0.5
# mixture-2 1903.9510   300    0.5
# Detection function starting values  = 3000
# Truncation distances
# x= 6450  ; y= 8850  ; z= 1903.951
#------------------------------------------------------
#Estimating parameters
#------------------------------------------------------
#Maximum likelihood results
# ------------------------------------------------------
#parameter point estimates = -349.812 226.903 0.29 1145.626 429.287 3023.626
#AIC = 3894.54
#------------------------------------------------------
}
detach(sightings)
}
\seealso{
  \code{\link{nupoint.env.gof}}, \code{optim},
  \code{\link{detectF}}
}
\keyword{misc}

