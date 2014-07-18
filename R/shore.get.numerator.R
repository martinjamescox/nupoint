#'Numerator of the environmental variable (depth in case of whale) log-likelihood. NB this is an internal function.
#'@param    pars   model parameters, paramter vector optimised using optim. e.g.: 
#'                           pars[1]=depth gradient par 1;
#'                           pars[2]=depth gradient par 2;
#'                           pars[3]=half normal rng detection function. NB pars unpacked using
#'                           par.unpack.F        
#'@param           z   seabed depth under cue sighting.
#'@param           dzdy rate of change of seabed depth at cue sighting.
#'@param           rd  radial distance from observer to cue.
#'@param           minz minimum z (depth) NB assumes z is positive.
#'@param           wx x-dimension truncation distance (for x-direction uniform sampling).
#'@param           wz   seabed depth truncation distance. NB assumes z is positive
#'@param           grad.type  cue depth disribution type e.g. "NORM"  - this is currently coded for
#'@param                           a radial detection function only
#'@param           det.type  radial detection function form.
#'@param           n     number of distributions in a multinomial distribution (default NULL).
#'@param           diagnostic  print diagnostic data.
#'@return   single value of the log-likelihood numerator.

shore.get.numerator <- function(pars,z,dzdy,rd, minz,wx,wz,grad.type,det.type,n=NULL, diagnostic=FALSE)
{
  #CALLS:     par.unpack.F; pi.z.f; 
  #unpack parameters
  #parList=par.unpack.F(grad.type,pars,n)
  #grad.pars=parList[[1]]; det.par=parList[[2]]
  grad.pars=pars[[1]]; det.par=pars[[2]]
  pi.z=pi.z.f(g.type=grad.type,pars=grad.pars,z=z,z.lim=c(minz,wz),n=n)
  
  if(diagnostic){
    print('Numerator')
    cat('grad pars=',grad.pars,'\n')
    cat('detect pars =',det.par,'\n')
    cat('pi_z = ',summary(pi.z),'\n')} #end of if diagnostic.
  detV=log(detectF(rd=rd,det.type=det.type,det.par=det.par))
  OUT = sum(detV + log(abs(dzdy)) + log(pi.z) + log(1/(2*wx)))
  return(OUT)
}