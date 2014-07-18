#'Log-likelihood for the density gradient with respect to a linear feasture.
#'
#'@param pars a single vector of starting parameters to be estimated.  Parameters must be ordered: environment density gradient first, followed by detection function. See environ.fit.f'
#'@param sight.x observation vector:  target x coordinate. This is cross-track distance for multi-beam echosounder(MBE) observations.'
#'@param swarm.z observation vector:  target y coordinate. This is water depth to the centre of a swarm for MBE observations.'
#'@param grid.density numerical integration grid density (number of elements in each (x and y) dimensions.'
#'@param td truncation distance in the y-dimension (depth for MBE)'
#'@param sb.z observation vector:  default NULL if no attenuation function is required (all transects extend to w). FOr MBE observations, This is water depth to the centre of a swarm.'
#'@param g.type parametric form of the environmental gradient.  See environ.fit.f.'
#'@param det.type parametric form of the detection function.  See detectF.'
#'@param n  number of normal mixture distributions in the MNORM function.
#'@param angularDetect logical: FALSE. Should angular detection be considered (not yet implemented).'
#'@param verbose Logical. Default FALSE. If TRUE parameter estimates and log-likelihood are printed to the console for each optim search iteration.'
#'@param mask area within the point transect where observation is possible e.g. MBE swath. created using swathInOutF in mbe.fit.f
#'@return -log lik.
shore.only <- function(pars,sight.x, swarm.z, grid.density,td,sb.z,g.type,det.type,n=NULL,verbose=TRUE,mask,xV,yV,angularDetect)
{

  pars=par.unpack.F(g.type,pars,n=n) #  #unpack parameters
  attenuationFlag=FALSE
  attenuationFlag[length(sb.z)>0]=TRUE
  top <- shore.get.numerator.mbe(pars, sight.x=sight.x,swarm.depth=swarm.z,
                                 w=td,seabed.z=sb.z,grad.type=g.type,det.type=det.type,n=n,
                                 angularDetect=angularDetect,x.max=max(xV))
  bottom <- shore.get.denominator.mbe(pars, num.detects=length(swarm.z),
                                      w=td,grad.type=g.type,det.type=det.type,n=n,
                                      mask=mask,xV=xV,yV=yV,attenuation=attenuationFlag,angularDetect=angularDetect,xmax=max(xV))
  shore.likhood <- top + bottom
  
  if(verbose){
    cat('current density gradient parameter estimates',pars[[1]],'\n')
    cat('current range detection function parameter estimate(s)',pars[[2]],'\n')
    cat('numerator',top,'\n')
    cat('denominator',bottom,'\n')
    cat('log-likelihood',shore.likhood,'\n')
  }
  shore.likhood[is.nan(shore.likhood)]=-9e99
  shore.likhood[is.infinite(shore.likhood)]=-9e99
  return(shore.likhood)
}