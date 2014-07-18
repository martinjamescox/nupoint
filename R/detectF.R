#'Probability of detection at a given range
#'
#'This function calculates the probability of detecting a target (e.g. krill
#'swarm or beaked whale) at a given distance from the observer.
#'
#'The probability of detecting a target at a given radial distance from an
#'observer can be modelled using the following parametric forms
#'(\code{det.type} argument) and parameters passed to the function using the
#'\code{det.par} argument:
#'
#'1) half-normal "HNORM", \code{det.par=c(sigma)}
#'
#'2) Hazard rate "HAZARD", \code{det.par=c(sigma,b)}
#'
#'NB If constructing an intital value vector, the \code{pars} argument in the
#'\code{\link{nupoint.env.fit}} function, \code{pars=c(...,detect function
#'parameters)} where \code{...} denotes environmental gradient parameters.
#'
#'@param rd \code{radial distance from observer to detected target.}
#'@param det.type \code{parametric form of the detection function
#'c("HNORM","HAZARD"). See details.}
#'@param det.par \code{vector of detection function parameters.}
#'@return Detection probability at range \code{rd}.
#'@section Reference: Buckland, S.T., Anderson,D. R., Burnham, K. P., Laake, J.
#'L., Borchers,D. L. and Thomas, L. (2001) Introduction to Distance Sampling.
#'Oxford: Oxford University Press, page 47.
#'@seealso \code{\link{nupoint.env.fit}}, \code{\link{nupoint.fit}}
#'@keywords misc
#'@export
#'@examples
#'
#'rVec=0:100
#'pVec=detectF(rd=rVec, det.type='HNORM', det.par=50)
#'plot(rVec,pVec,type='l',xlab='Distance',ylab='Pr(detect)')
#'
detectF <- function(rd,det.type,det.par)
{
  #20120903: detection function
  #INPUTS: det.type = detection function type: c('HNORM','HAZARD',
  #                   'EXP','EXPQ','LOGISTIC')
  #         det.par = detection function parameter vector
  #         rd = radial distance, but could also be angle.
  #RETURNS: probability of detection for rd | det.par
  k=1
  if(length(det.par==3))
    k=det.par[3]
  out=switch(det.type,
       HNORM =exp(-rd^2/(2*det.par[1]^2)),
       HAZARD = k* (1-exp(-(rd/det.par[1])^-det.par[2])),
       EXP = k* exp(-(rd/det.par[1])**det.par[2]),
       EXPQ = k* exp(-det.par[1]*rd-det.par[2]*rd**2),
       LOGISTIC = 1/(1+exp(-(rd-det.par[1])/det.par[2])))
  out[out>1]=1
  return(out)
}
