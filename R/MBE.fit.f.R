MBE.fit.f <- function(pars,sight.x,sight.y,sight.sb.z=NULL, w, theta.max,grad.type, det.type, 
                   grid.density=100, angularDetect=FALSE,verbose=FALSE,lower.b,upper.b,optim.control=NULL){
  #### THIS FUNCTION IS DOCUMENTED IN NUPOINT #####
  #20121107: uses optim to obtain MLE for multibeam type problems.
  #INPUTS:
  #pars: 'a single vector of starting parameters to be estimated.  Parameters must be ordered: environment density gradient first, followed by detection function. See environ.fit.f'
  #sight.x: 'observation vector:  target x coordinate. This is cross-track distance for multi-beam echosounder(MBE) observations.'
  #sight.y: 'observation vector:  target y coordinate. This is water depth to the centre of a swarm for MBE observations.'
  #sight.sb.z: 'observation vector:  default NULL if no attenuation function is required (all transects extend to w). FOr MBE observations, This is water depth to the centre of a swarm.'
  #w: 'truncation distance in the y-dimension (depth for MBE)'
  #theta.max: 'maximum observation angle rad, (MBE swath width).  Maximum angle is pi/2.'
  #grad.type: 'parametric form of the environmental gradient.  See environ.fit.f.'
  #det.type: 'parametric form of the detection function.  See detectF.'
  #grid.density: 'numerical integration grid density (number of elements in each (x and y) dimensions.'
  #angularDetect: 'logical: FALSE. Should angular detection be considered (not yet implemented).'
  #verbose: 'Logical. Default FALSE. If TRUE parameter estimates and log-likelihood are printed to the console for each optim search iteration.'
  #lower.b: 'Vector of lower parameter space bounds.  Length must equal that of the pars argument. '
  #upper.b: 'Vector of upper parameter space bounds.  Length must equal that of the pars argument. '
  #optim.control: 'control list object to be passed to the optim function.  See optim help and details.'
  #RETURNS:
  #optim output with AIC.
  
  #checks
  #no seabed attenuation function in workspace
  if(length(sight.sb.z)>0){
    if(!('atten.f'%in% ls(envir=.GlobalEnv))){
      warning('Seabed attenuation function missing when seabed at observation data passed to mbe.fit.f in the sight.sb.z ARG.')
      cat('WARNING: seabed attenuation function atten.f not in workspace.\n')
      cat('Either run make.atten.f to create atten.f or remove seabed data from mbe.fit.f (sight.sb.z=NULL)\n')
      break
    }
  }
  #no angular detection function in workspace
  if(angularDetect){
    if(!('angDetF'%in% ls(envir=.GlobalEnv))){
      warning('Angular detection function, angDetF, missing when angularDetect=TRUE in mbe.fit.f')
      cat('WARNING: Angular detection function, angDetF, not in workspace.')
      cat('Either create angDetF or set angularDetect=FALSE in mbe.fit.f','\n')
      break
    }
  }
  if(any(sight.y>w)){
    warning('One or one observations exceed the truncation distance, w.')
    return('Observation(s)_exceed_w')
  }
  bound.chk=bound.chk.f(pars,lower=lower.b,upper=upper.b,f.n='mbe.fit.f')
  if(bound.chk) return('chk_parameter_boundary')
  
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
  
  #print model fit
  cat('---------------------------','\n')
  cat('Perpendicular density gradient likelihood settings','\n')
  cat('---------------------------','\n')
  cat('Depth preference parametric form:',grad.type,'\n')
  cat('range detection function, g(r), parametric form:',det.type,'\n')
  cat('Parameter starting values =',pars,'\n')
  cat('Angular detection function applied:',angularDetect,'\n')
  if(angularDetect) print(angDetF)
  cat('Seabed attenuation function applied:',ifelse(length(sight.sb.z)>0,'TRUE','FALSE'),'\n')
  if(length(sight.sb.z)>0) print(atten.f)
  cat('Truncation distance =',w,'\n')
  cat('Swath (max) angle, rad =',theta.max,'\n')
  cat('Numerical integration grid density=',grid.density,'\n')
  cat('---------------------------','\n')
  
  xV=yV= seq(0,w,length=grid.density)+(0.5*w/grid.density)  
  xV=xV[-length(xV)]; yV=yV[-length(yV)]
  
  mask = outer(xV,yV,swathInOutF,w,theta=theta.max)
  grad.likhood=NULL
  cat('Estimating parameters \n')
  grad.likhood <- optim(pars, shore.only,sight.x=sight.x,swarm.z=sight.y, grid.density=grid.density,td=w,sb.z=sight.sb.z,
                        g.type=grad.type,det.type=det.type,verbose=verbose,
                        mask=mask,xV=xV,yV=yV,angularDetect=angularDetect,
                        hessian=TRUE,method="L-BFGS-B", lower=lower.b, upper=upper.b,
                        control=control.list)
  grad.likhood$hessian=-grad.likhood$hessian
  if(grad.likhood$convergence!=0)
    warning('environ.fit.f: convergence error code returned by optim')
  if(any(grad.likhood$par < 1.01*pars  & grad.likhood$par > 0.99*pars ))
    warning('environ.fit.f: one or more parameter estimates within 1 percent of starting value.  Check convergence.')
  if(any(floor(grad.likhood$par)==lower.b))
    warning('environ.fit.f: one or more parameter estimates at lower bound.  Check convergence.')
  if(any(floor(grad.likhood$par)==upper.b))
    warning('environ.fit.f: one or more parameter estimates at upper bound.  Check convergence.')
  
  cat('--------------------------- \n')
  cat('Maximum likelihood results \n')
  cat('---------------------------\n')
  cat('parameter estimates \n')
  cat(round(grad.likhood$par,3),'\n')
  grad.likhood$AIC=-2*grad.likhood$value+2*length(grad.likhood$par)
  cat('AIC =',round(grad.likhood$AIC,2),'\n')
  cat('---------------------------\n')
  return(grad.likhood)
  
 
				   }
