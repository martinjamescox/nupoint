  #20121107: 
#' Calculate value at a given numerical integration cell in the likelihood denominator.
#' @param x and y coordinates
#' @param pars model parameters, pars[1]=depth gradient par 1; pars[2]=depth gradient par 2; pars[3]=half normal rng detection function
#' @param trunc.dist  truncation distance (y-dimension)
#' @param grad.type parametric form of the environmental gradient.  See environ.fit.f.
#' @param det.type  parametric form of the detection function.  See detectF.
#' @param  n  number of normal mixture distributions in the MNORM function.
#' @param attenuation Logical TRUE - use attenuation function (NB attenuation function must be called atten.f); FALSE - no attenuation.
#' @param angularDetect  logical: FALSE. Should angular detection be considered (not yet implemented).
#' @return value at denominator numerical integration grid cell.
shore.calc.integral.grid <- function(x,y,pars, trunc.dist,grad.type,det.type,n=NULL,attenuation,angularDetect)
{
  r=sqrt(x**2 + y**2)
  seabedAttenuation=1
  if(attenuation)
    seabedAttenuation=atten.f(y)
  angDet=1
  if(angularDetect)
    angDet=angDetF(x=atan(x/y))
  f=detectF(rd=r,det.type=det.type,det.par=pars[[2]]) * r * pi.z.f(g.type=grad.type,pars=pars[[1]],z=y,z.lim=c(0,trunc.dist),n=n) * angDet * seabedAttenuation
  return(f)
}
