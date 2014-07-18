#'Unpack parameters so that optim can run using a single parameter vector.
#'@param  grad.type  gradient type (Normal etc). c('NORM','BETA','LOGNORM','UNIFORM','TWEEDIE','MNORM')
#'@param  pars  parameter vector containing all parameters to be estimated by optim.  pars vector
#'               NB structure for MNORM \eqn{\mu_1,\sigma_1,\alpha_1},...,\eqn{\mu_n,\sd_n},then detection 
#'                 function parameters.  there is no alpha_n.
#'@param  n  number of distributions in a multinomial distribution (default NULL). 
#'@return list([1]=density gradient parameters; [2]=detection function parameters).

par.unpack.F <- function(grad.type,pars,n){
  if(grad.type=="MNORM") {
    if(n<2)
      warning('shore.get.numerator: insufficient number of distributions specified in mixture normal perference distribution')
    
    grad.pars=pars[1:((3*n)-1)] #asummes Mx distribution has two parameters
    det.par <- pars[(3*n):length(pars)]
    
    }else{ #single distribution
    grad.pars <- pars[1:2]
    det.par <- pars[3:length(pars)]
    if(grad.type=="UNIFORM"){
      grad.pars=NA
      det.par=pars}
    if(grad.type =='TWEEDIE'){
      grad.pars=pars[1:3] #mu,phi,power
      det.par <- pars[4:length(pars)]}
    }
  return(list(grad.pars,det.par))
}