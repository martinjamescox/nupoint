#' Calculate log-likelihood numerator
#' @param pars    model parameters, pars[1]=depth gradient par 1; pars[2]=depth gradient par 2; pars[3]=half normal rng detection function
#' @param sight.x  observation vector:  target x coordinate. This is cross-track distance for multi-beam echosounder(MBE) observations.'
#' @param swarm.depth depth at geometric swarm centre (m)
#' @param w     truncation distance (m)
#' @param seabed.z  seabed depth under geometric swarm centre (m)
#' @param grad.type   parametric form of the environmental gradient.  See environ.fit.f.
#' @param det.type   parametric form of the detection function.  See detectF.
#' @param n    number of normal mixture distributions in the MNORM function.
#' @param angularDetect logical: FALSE. Should angular detection be considered (not yet implemented).
#' @param x.max    maximum distance in the x-dimension (w * sin(theta.max))
#' @return numerator for the log-likelihood
shore.get.numerator.mbe <- function(pars,sight.x,swarm.depth, w,seabed.z,grad.type,det.type,n=NULL,angularDetect,x.max)
{
  radial.dist=sqrt(sight.x**2+swarm.depth**2)
  seabedAttenuation=0
  if(length(seabed.z)>0)
    seabedAttenuation=log(atten.f(seabed.z) )
  angDetV=0
  if(angularDetect)
    angDetV=log(angDetF(x=atan(sight.x/swarm.depth)))
  numerator=sum(log(detectF(rd=radial.dist,det.type=det.type,det.par=pars[[2]])) + log(pi.z.f(g.type=grad.type,pars=pars[[1]],n=n,z=swarm.depth,z.lim=c(0,w))) + 
    angDetV + seabedAttenuation  + log(1/x.max) + log(radial.dist))
  return(numerator)
}
