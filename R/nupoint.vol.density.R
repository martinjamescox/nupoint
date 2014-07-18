#'Calculate volumetric density of targets for a parallel density survey
#'
#'This function calculates volumetric density of shoals/targets using the
#'probability of detecting a target (e.g. hare, fish shoal, krill swarm) in a
#'rectangular region.
#'
#'
#'@param n.seen \code{number of targets detected.}
#'@param L \code{total line transect length.}
#'@param w \code{vertical dimension truncation distance.}
#'@param theta.max \code{maximum observation angle rad, (for multi-beam
#'echosounders, this is swath width).  Maximum angle is pi/2.}
#'@param P.star \code{probability of detecting a target, e.g. krill swarm, in a
#'rectangular area (dimensions:w * w * sin(theta.max)).}
#'@return Target volumetric density (e.g. number of krill swarms per cubic
#'metre of seawater) within a survey region.
#'@seealso \code{\link{nupoint.p.star}}, and Figure 1 and Section 5.1 of the
#'\code{parallel.pdf} vignette
#'@export
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'nupoint.vol.density(n.seen=251, L=11*2.5*1e3, w=100, theta.max=pi/3, P.star=0.079)
#'}
#'
nupoint.vol.density <- function(n.seen,L,w,theta.max,P.star) n.seen / L * 1/(2*w*sin(theta.max)*w*P.star)


