#'Log-likelihood forenvironmental variable (depth in case of whale) NB this function is optimised by optim.
#'@param    pars     model parameters, paramter vector optimised using optim. e.g.: 
#'                           pars[1]=depth gradient par 1;
#'                           pars[2]=depth gradient par 2;
#'                           pars[3]=half normal rng detection function. NB (i) pars unpacked using
#'                           par.unpack.F, (ii) structure for MNORM \eqn{\mu_1,\sigma_1,\alpha_1,...,\mu_n,\sd_n},
#'                               there is no \eqn{alpha_n}
#'@param           z       seabed depth under cue sighting.
#'@param           dzdy     rate of change of seabed depth at cue sighting.
#'@param           rd      radial distance from observer to cue.
#'@param           z.mat   gridded seabed depth throughout survey area (regular grid; positive depths).
#'@param           dzdy.mat   rate of change of seabed depth throughout survey area.
#'@param           rd.mat     matrix of radial distances from observer to each grid cell in survey area.
#'@param           minz     minimum z (depth) NB assumes z is positive.
#'@param           wx     x-dimension truncation distance (for x-direction uniform sampling).
#'@param           wy     y-dimension truncation distance.
#'@param           wz     seabed depth truncation distance. NB assumes z is positive
#'@param           grad.type  cue depth disribution type e.g. "NORM"  - this is currently coded for
#'                           a radial detection function only
#'@param           n    number of distributions in a multinomial distribution (default NULL).
#'@param           diagnostic  print diagnostic data from the hore.get.numerator and shore.get.denominator.
#'@param           verbose   prints the log-likelihood and parameter estimates for s,                         
#'@return   single value of the log-likelihood denominator.
environ.mle.f <- function(pars,z, rd,dzdy,z.mat,dzdy.mat,rd.mat,minz,wx,wy,wz,g.type,det.type,
         verbose=FALSE,n=NULL,diagnostic=TRUE)
{
  #CALLS:     shore.get.numerator; shore.get.denominator
  pars=par.unpack.F(grad.type=g.type,pars,n)
  top <- shore.get.numerator(pars=pars,z=z,dzdy=dzdy,rd=rd, minz=minz,wz=wz,
                             wx=wx,grad.type=g.type,det.type=det.type,n=n,
                             diagnostic=diagnostic)
  bottom <- shore.get.denominator(pars=pars,num.detects=length(z),
                                  wx=wx,wy=wy,minz=minz,wz=wz,grad.type=g.type,
                                  z.mat=z.mat,dzdy.mat=dzdy.mat,
                                  rd.mat=rd.mat,det.type=det.type,n=n,
                                  diagnostic=diagnostic)
  shore.likhood <- top + bottom
  
  if(verbose){
    print(unlist(pars))
    print(shore.likhood)}
  if(is.infinite(shore.likhood)) shore.likhood= -9e99
  if(is.nan(shore.likhood)) shore.likhood= -9e99
  return(shore.likhood)
}