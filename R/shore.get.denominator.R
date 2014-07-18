#'Denominator of the environmental variable (depth in case of whale) log-likelihood. NB this is an internal function.
#'@param    pars    model parameters, paramter vector optimised using optim. e.g.: 
#'                           pars[1]=depth gradient par 1;
#'                           pars[2]=depth gradient par 2;
#'                           pars[3]=half normal rng detection function. NB pars unpacked using
#'                           par.unpack.F 
#'@param           num.detects  number of detections.
#'@param           z.mat         gridded seabed depth throughout survey area (regular grid; positive depths).
#'@param           dzdy.mat      rate of change of seabed depth throughout survey area.
#'@param           rd.mat      matrix of radial distances from observer to each grid cell in survey area.
#'@param           minz         minimum z (depth) NB assumes z is positive.
#'@param           wx           x-dimension truncation distance (for x-direction uniform sampling).
#'@param           wy           y-dimension truncation distance.
#'@param           wz        seabed depth truncation distance. NB assumes z is positive
#'@param           grad.type  cue depth disribution type e.g. "NORM"  - this is currently coded for
#'                           a radial detection function only
#'@param           n         number of distributions in a multinomial distribution (default NULL).
#'@param           diagnostic  print diagnostic data.
#'@return   single value of the log-likelihood denominator.
shore.get.denominator <- function(pars,num.detects,wx,wy,minz,wz,grad.type,
         z.mat,dzdy.mat,rd.mat,det.type,n=NULL, diagnostic=FALSE)
{

  #CALLS:     par.unpack.F; detectF; pi.z.f; 
  
  #unpack parameters
  # parList=par.unpack.F(grad.type,pars,n)
  #grad.pars=parList[[1]]; det.par=parList[[2]]
  grad.pars=pars[[1]]; det.par=pars[[2]]
  det.mat=detectF(rd=rd.mat,det.type=det.type,det.par=det.par)
  pi.z=pi.z.f(g.type=grad.type,pars=grad.pars,z=z.mat,z.lim=c(minz,wz),n=n,
              mnormVerboseFlag=diagnostic)
  
  if(diagnostic){
    print('Denominator')
    cat('gradient type=',grad.type,'\n')
    cat('z lim =',c(minz,wz),'\n')
    cat('z matrix =',summary(as.vector(z.mat)),'\n')  
    cat('grad pars=',grad.pars,'\n')
    cat('detect pars =',det.par,'\n')
    cat('det.mat mean =',mean(det.mat,na.rm=TRUE),'\n')
    cat('pi.z summary =',summary(as.vector(pi.z)),'\n')}
  denominator <- num.detects * -log(sum( pi.z * det.mat*abs(dzdy.mat),na.rm=TRUE)/(2*wx))####
  return(denominator)
}