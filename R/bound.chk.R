#'Check starting values are within parameter space bounds when using ,method="L-BFGS-B" in optim()
#'@param pars  initial parameter estimates
#'@param        lower  lower bound
#'@param        upper  upper bound
#'@param        f.n  character string of function within which the bounds are being checked e.g. 'mbe.fit.f for multibeam.
#'@return warning flag  logical FALSE = estimates are fine. TRUE = estimates outside upper or lower bounds.

bound.chk.f <- function(pars,lower,upper,f.n)
{
warning.flag=FALSE
  LC=which(pars<=lower)
  UC=which(pars>=upper)
  if(length(LC)>0)
  {
    warning(cat('Parameter starting values at or less than lower bound in',f.n,'\n'))
    warning.flag=TRUE
  }
  if(length(UC)>0)
  {
    warning(cat('Parameter starting values at or greater than upper bound',f.n,'\n'))
    warning.flag=TRUE
  }
  return(warning.flag)
}