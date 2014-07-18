#'Density of mixture of normal distributions
#'@param  z observation (seabed depth in the case of the whales; depth of swarm in case of multibeam)
#'@param  pars parameter vector:\eqn{\mu_1,\sigma_1,\alpha_1,...,\mu_n,\sd_n}, there is no alpha_n.
#'@param  z.lim z-dimension (whale) or y-dimension (multi-beam) truncation distance - used to scale the density.
#'@param  n number of multivariate normal distributions
#'@param  seabedAtObs for multibeam (krill) only. this allows density to be rescaled on an
#'                     observation-by-observation basis.  must be the same length as z.
#'@param  plot TRUE/FALSE plot density of mixture normal distributions
#'@param  verbose print parameter vectors split into mean, variance and weights.
#'@return scaled density

mnorm.f <- function(z,pars,z.lim,n,seabedAtObs=NULL,plot=FALSE,verbose=FALSE){
  require(nor1mix,quietly = TRUE)
  #basic error checking for incorrect parameter vector length
  if(((length(pars)+1)/3)!=n)
    warning('mnormf: number of parameters to number of distribution mismatch')
  
  #set w_z truncation distance to maximum observed seabed depth z max.
  zMax=z.lim[2]
  if(length(seabedAtObs)>0){
    #error checking for seabed at obs. vector (seabedAtObs) and depth of observation (z) vector
    if(length(seabedAtObs)!=length(z))
      warning('mnorm.f: mismatch in seabed depth and obs. vector lengths.')
    zMax = seabedAtObs
  } #end if length(seabedAtObs)>0
  
  alphas=pars[seq(3,3*(n-1),3)]   #extract the alpha values from the parameter vector
  
  #implements the weight estimation idea detailed in Dave Miller's thesis:
  FalphaV=pgamma(cumsum(exp(alphas)),3,2)
  wts=c(FalphaV[1],diff(FalphaV))
  wts[n]=1-sum(wts[1:(n-1)]) #calculate the weight for the nth distribution
  
  norMixObj=norMix(mu=pars[seq(1,(n*3)-2,3)], 
                   sig2 = pars[seq(2,(n*3)-1,3)]**2, w = wts) #creates a normal mixture object.
  
  den=dnorMix(x=z, obj=norMixObj)/diff(pnorMix(q=c(z.lim[1],zMax), obj=norMixObj))
  if(verbose){cat('mnorm.f pars=',pars,'\n')
              cat('mnorm.f alphas=',alphas,'\n')
              cat('mnorm.f mu=',pars[seq(1,(n*2)-1,2)],'\n')
              cat('mnorm.f sd=',pars[seq(2,(n*2),2)],'\n')
              
              cat('mnorm.f mixing parameters=',wts,'\n')
              print(norMixObj)
  }
  if(plot) plot(norMixObj,xlim=c(0,zMax))
  return(den)
}