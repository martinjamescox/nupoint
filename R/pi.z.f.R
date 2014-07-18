#' Scaled density for a variety of distributions
#'@param g.type  gradient type (Normal etc). c('NORM','BETA','LOGNORM','UNIFORM','TWEEDIE','MNORM')
#'@param pars  parameter vector containing only the gradient distribution parameters. Within the likelihood this is provided by par.unpack.F 
#'@param  z  depth
#'@param z.lim  vector (min,max) seabed or dimension of interest.
#'@param seabedAtObs  for multibeam (krill) only. this allows density to be rescaled on an
#                       observation-by-observation basis.  must be the same length as z.
#'@param n  number of multivariate normal distributions
#'@param  mnormPlotFlag  argument to be passed mnorm.f (diagnostic) - should mvnorm be plotted
#'@param mnormVerboseFlag  argument to be passed mnorm.f (diagnostic)
#'@return scaled density for a given depth (dimension of interest).
#'CALLS: tweedief; mnorm.f
#'@details the scaling of density could be removed for the by observation calculation of the integration
#' grid.  that way the integration grid could be calculated for each parameter estimate, then trimmed to 
#' seabed depth and then scaled.
#'variable maxz can be overwritten by seabed at each observation in multibeam case
pi.z.f <- function(g.type,pars,z,z.lim,seabedAtObs=NULL,n=NULL,mnormPlotFlag=FALSE,
                mnormVerboseFlag=FALSE){ 
    maxz=z.lim[2] 
  if(length(seabedAtObs)!=0){
    maxz=seabedAtObs
  #catch length mismatch in   
  if(length(z)!=length(seabedAtObs))
    warning('pi.z.f: length mismatch in object depth and seabed depth vectors')
    }
  switch(g.type,
         NORM = dnorm(z,pars[1],pars[2])/diff(pnorm(c(z.lim[1],maxz),pars[1],pars[2])),
         BETA = dbeta(z/z.lim[2],pars[1],pars[2])/diff(pbeta(c(z.lim[1],maxz)/z.lim[2],pars[1],pars[2])),
         LOGNORM = dlnorm(z,pars[1],pars[2])/diff(plnorm(c(z.lim[1],maxz),pars[1],pars[2])),
         UNIFORM = dunif(z,z.lim[1],z.lim[2]),
         TWEEDIE = tweedief(z,pars,zlim=c(z.lim[1],maxz)),
         MNORM = mnorm.f(z,pars,z.lim,n,seabedAtObs=seabedAtObs,
                         plot=mnormPlotFlag,verbose=mnormVerboseFlag))}