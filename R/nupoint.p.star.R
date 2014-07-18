#'Calculate the probability of detecting a target within the covered area given
#'a parallel density gradient.
#'
#'This function calculates the probability of detecting a target, e.g. krill
#'swarm, within the sampled area. Parameters estimated using
#'\code{\link{nupoint.fit}} are typically used in the \code{pars.hat} argument
#'
#'This function calculates the probability of detecting a target (e.g. krill
#'swarm, fish shoal, or single fish) in the area covered by a multi-beam
#'echosounder. Typically, this function will calculate the probability of
#'detecting a target using the parameter estimates calculated in
#'\code{\link{nupoint.fit}}.
#'
#'@param pars.hat \code{a single vector of parameters used to estimate the
#'probablity of detecting target within covered area.  Parameters must be
#'ordered:environment density gradient first, followed by detection function.
#'These parameters will typically be estimated using nupoint.fit}
#'@param grad.type \code{parametric form of the target object density function.
#'See nupoint.fit and nupoint.env.fit.}
#'@param det.type \code{parametric form of the detection function.  See
#'detectF.}
#'@param n \code{Default NULL. Use only when grad.type=MNORM to specify an
#'integer number of distributions in the normal mixture option for the
#'parametric form of the depth preference function.}
#'@param w \code{truncation distance in the y-dimension (depth for a multi-beam
#'echosounder)}
#'@param theta.max \code{maximum observation angle radians, (multi-beam
#'echosounder swath width).  Maximum angle is pi/2.}
#'@param grid.density \code{numerical integration grid density (number of
#'elements in each (x and y) dimensions).}
#'@param attenuation \code{logical:default FALSE.  If TRUE, nupoint.p.star uses
#'the attenuation function, atten.f, in the R workspace.  See make.atten.f.}
#'@param angularDetect \code{logical:default FALSE. Should angular detection be
#'considered (not yet implemented).}
#'@return Probability of detecting a target in the covered region.
#'@seealso \code{\link{nupoint.fit}}
#'@export
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'nupoint.p.star(pars.hat=c(50,20,60),grad.type='NORM',det.type='HNORM',n=NULL,w=100,
#'  theta.max=pi/3,grid.density=100,
#'  attenuation=FALSE,angularDetect=FALSE)
#'}
#'
nupoint.p.star <- function(pars.hat,grad.type,det.type,n=NULL,w,theta.max,grid.density=100,
                       attenuation=FALSE,angularDetect=FALSE)
{
  #20130108 nupoint.p.star
  #we start by calling the denominator function and passing
  #the fitted parameters and assumed parameteric forms of the
  #detection and density gradient functions. 
  #NB angular detection is not currently implemented, so the angularDetect argument
  #must remain angularDetect=FALSE
  browser()
  shore.calc.integral.grid.P=
  function(x,y,pars, trunc.dist,grad.type,det.type,attenuation,angularDetect,x.max)
  {
    r=sqrt(x**2 + y**2)
    seabedAttenuation=1
    if(attenuation)
      seabedAttenuation=atten.f(y)
    angDet=1
    if(angularDetect)
      angDet=angDetF(x=atan(x/y))
    f=detectF(rd=r,det.type=det.type,det.par=pars[[2]]) *  
      pi.z.f(g.type=grad.type,pars=pars[[1]],z=y,z.lim=c(0,trunc.dist)) * 
      angDet * seabedAttenuation 
    return(f)
  }
  
  #checks
  #no seabed attenuation function in workspace
  if(attenuation){
    if(!('atten.f'%in% ls(envir=.GlobalEnv))){
      warning('Seabed attenuation function, atten.f, missing from workspace.')
      break
    }
  }
  #end checks.
  parL=par.unpack.F(grad.type=grad.type,pars=pars.hat,n=n)
  yV= seq(0,w,length=grid.density)+(0.5*w/grid.density)  
  xV=seq(0,w*sin(theta.max),length=grid.density)+(0.5*w/grid.density)  
  xV=xV[-length(xV)]; yV=yV[-length(yV)]
  mask = outer(xV,yV,swathInOutF,w,theta=theta.max)
  temp <- outer(xV, yV, shore.calc.integral.grid.P, pars=parL,trunc.dist=w,grad.type=grad.type,det.type=det.type,
                attenuation=attenuation,angularDetect=angularDetect)
  temp = temp * mask
  temp = temp * (1/w*sin(theta.max))
  PStar <- mean(temp,na.rm=TRUE)*w*w*sin(theta.max)
  return(PStar)
}
