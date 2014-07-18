#'Numerical integration for likelihood denominator
#'
#'@param pars   model parameters, pars[1]=depth gradient par 1; pars[2]=depth gradient par 2; pars[3]=half normal rng detection function
#'@param num.detects  number of animals detected
#'@param w   truncation distance (y-dimension)
#'@param grad.type  parametric form of the environmental gradient.  See environ.fit.f.
#'@param det.type  parametric form of the detection function.  See detectF.
#'@param n  number of normal mixture distributions in the MNORM function.
#'@param mask  area within the point transect where observation is possible e.g. MBE swath. created using swathInOutF in mbe.fit.f
#'@param xV  vector of x dimension coordinates for integration grid.
#'@param yV  vector of y dimension coordinates for integration grid.
#'@param attenuation  Logical TRUE - use attenuation function (NB attenuation function must be called atten.f); FALSE - no attenuation.
#'@param angularDetect  logical: FALSE. Should angular detection be considered (not yet implemented).
#'@return likelihood denominator.
shore.get.denominator.mbe <- function(pars,num.detects,w,grad.type,det.type,n=NULL,mask,xV,yV,attenuation,angularDetect,xmax)   
{
  temp <- outer(xV, yV, shore.calc.integral.grid, pars=pars,trunc.dist=w,grad.type=grad.type,det.type=det.type,n=n,
                attenuation=attenuation,angularDetect=angularDetect)
  temp = temp *(1/xmax) * mask
  denominator <- mean(temp,na.rm=TRUE)*w*(1/xmax)
  denominator <- num.detects * log(1/denominator)
  return(denominator)
}