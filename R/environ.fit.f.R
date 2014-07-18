environ.fit.f <- function(pars,z, rd,dzdy,z.mat,dzdy.mat,rd.mat,minz,wx,wy,wz,grad.type,det.type,
                       verbose=FALSE,n=NULL,diagnostic=FALSE,lower.b,upper.b,optim.control=NULL)
{
   bound.chk=bound.chk.f(pars,lower=lowerB,upper=upperB,f.n='mbe.fit.f')
  if(bound.chk) return('chk_parameter_boundary')
 
  #print model fit
  cat('------------------------------------------------------\n')
  cat('Environmental gradient likelihood settings\n')
  cat('------------------------------------------------------\n')
  cat('Environment preference parametric form:',grad.type,'\n')
  cat('range detection function, g(r), parametric form:',det.type,'\n')
  if(grad.type=='MNORM') {
    cat('Mixture of normal distributions with starting values =\n')
    print(matrix(c(pars[-((n*3):length(pars))],1-sum(pars[(1:(n-1)*3)])),nrow=nDist,byrow=TRUE,
                 dimnames=list(paste('mixture',1:n,sep='-'),c('mu','sigma','weight'))))
    cat('Detection function starting values  =',pars[((n*3):length(pars))],'\n')
    pars[(1:(n-1)*3)]=log(qgamma(pars[(1:(n-1)*3)],3,2)) #scale wts onto cumuluative gamma dist'n
  } else {
    cat('Parameter starting values =',pars,'\n')}
  cat('Truncation distances\n')
  cat('x=',wx,' ; y=',wy,' ; z=',wz,'\n')
  cat('------------------------ ------------------------------\n')
  
  parscales=10^round(log10(abs(pars))) #parameter scales
  
  #optim control
  control.list=list(fnscale=-1700,parscale=parscales,maxit=200)
  if(length(optim.control)>0){ #instance where user has defined optim control settings
    control.list=optim.control
    if(any(names(control.list) %in% 'fnscale') ==FALSE) #check user has specified fnscale and parscales
      control.list$fnscale=-1700
    if(any(names(control.list) %in% 'parscales') ==FALSE)
      control.list$parscales=parscales
    if(any(names(control.list) %in% 'maxit') ==FALSE)
      control.list$maxit=200
    
  } #end of user defined optim control list handling
  
  grad.likhood=NULL
  cat('Estimating parameters\n')
  grad.likhood <- optim(pars,environ.mle.f,
                        z=z, 
                        rd=rd,
                        dzdy=dzdy,
                        z.mat=z.mat,
                        dzdy.mat=dzdy.mat,
                        rd.mat=rd.mat,
                        minz=minz,
                        wx=wx,
                        wy=wy,
                        wz=wz,
                        g.type=grad.type,
                        det.type=det.type,
                        verbose=verbose,
                        n=n,
                        diagnostic=diagnostic,
                        hessian=TRUE,method="L-BFGS-B", lower=lower.b,
                        upper=upper.b,
                        control=control.list)
  if(grad.likhood$convergence!=0)
    warning('environ.fit.f: convergence error code returned by optim')
  if(any(pars==grad.likhood$par))
    warning('environ.fit.f: one or more parameter estimates equals starting value.  Check convergence.')
  if(any(pars==lower.b))
    warning('environ.fit.f: one or more parameter estimates at lower bound.  Check convergence.')
  if(any(pars==upper.b))
    warning('environ.fit.f: one or more parameter estimates at upper bound.  Check convergence.')
  
  
  cat('------------------------------------------------------ \n')
  cat('Maximum likelihood results \n')
  cat('------------------------------------------------------ \n')
  cat('parameter point estimates =',round(grad.likhood$par,3),'\n')
  grad.likhood$AIC=-2*grad.likhood$value+2*length(grad.likhood$par)
  cat('AIC =',round(grad.likhood$AIC,2),'\n')
  cat('------------------------------------------------------ \n')
  #additional model metrics
  return(grad.likhood)
}