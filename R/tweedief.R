#'Density of a scaled Tweedie distribution.
#'@param  x  observation (seabed depth in the case of the whales; depth of swarm in case of multibeam)
#'@param         pars  parameter vector for the distribution (does not include the detection function)
#'                 pars[1] mu; pars[2] phi, and pars[3] power.
#'@param         zlim  z-dimension (whale) or y-dimension (multi-beam) truncation distance - used to scale
#'                 the density.
#'@return scaled density of the tweedie distribution for parameters par
#'@note this function is required because passing an NA into dtweedie causes an error.  This function
#'         returns an NA when it is passed an NA in the x argument.
tweedief <- function(x,pars,zlim){
  #DEPENDS: package tweedie
    d=rep(NA,length(x))
  d[is.na(x)==FALSE]=dtweedie(x[is.na(x)==FALSE],mu=pars[1],phi=pars[2],power=pars[3])/ #mu,phi,power
    diff(ptweedie(q=zlim, mu=pars[1],phi=pars[2],power=pars[3])) #
  return(d)
}