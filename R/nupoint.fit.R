#'Maximum likelihood estimate of parameters describing parallel density
#'gradient and range detection functions
#'
#'This function calculates maximum likelihood estimates of parameters for both
#'a parallel density gradient and detection function, for which there are a
#'variety of options for parametric forms.
#'
#'This likelihood function is for use when determining a density gradient with
#'respect to a parallel feature e.g. road (Marques et al., 2010 and sea surface
#'(Cox et al., 2011).
#'
#'@param pars \code{a single vector of starting parameters to be estimated.
#'Parameters must be ordered:environment density gradient first, followed by
#'detection function. See nupoint.env.fit.}
#'@param sight.x \code{observation vector: target x coordinate. This is
#'cross-track distance for multi-beam echosounder(MBE) observations.}
#'@param sight.y \code{observation vector: target y coordinate. This is water
#'depth at the centre of an aggregation for MBE observations.}
#'@param sight.z \code{observation vector: default NULL if no attenuation
#'function is required (all transects extend to w). For MBE observations, this
#'is water depth to the centre of a swarm.}
#'@param w \code{truncation distance in the y-dimension (depth for MBE)}
#'@param theta.max \code{maximum observation angle radians, (MBE swath width).
#'Maximum angle is pi/2.}
#'@param grad.type \code{parametric form of the environmental gradient.  See
#'nupoint.env.fit.}
#'@param det.type \code{parametric form of the detection function.  See
#'detectF.}
#'@param n \code{Default NULL. Use only when grad.type=MNORM to specify an
#'integer number of distributions in the normal mixture option for the
#'parametric form of the environmental gradient.  See details.}
#'@param grid.density \code{numerical integration grid density (number of
#'elements in each (x and y) dimensions).}
#'@param angularDetect \code{logical:FALSE. Should angular detection be
#'considered (not yet implemented).}
#'@param verbose \code{logical. Default FALSE. If TRUE parameter estimates and
#'log-likelihood are printed to the console for each optim search iteration.}
#'@param lower.b \code{vector of lower parameter space bounds.  Length must
#'equal that of the pars argument. }
#'@param upper.b \code{vector of upper parameter space bounds.  Length must
#'equal that of the pars argument. }
#'@param optim.control \code{control list object to be passed to the optim
#'function.  See optim help and nupoint.env.fit.}
#'@return output as per optim and AIC.
#'@seealso \code{\link{nupoint.env.fit}}, \code{\link{detectF}}, \code{optim}
#'@references Cox, M.J., Borchers, D.L., Demer, D.A., Cutter, G.R., and
#'Brierley, A.S. (2011). Estimating the density of Antarctic krill (Euphausia
#'superba) from multi-beam echo-sounder observations using distance sampling
#'methods. Journal of the Royal Statistical Society: Series C (Applied
#'Statistics), vol 60, part 2, pp. 301–316.
#'
#'Marques, T.A. , Buckland, S.T. , Borchers, D.L. , Tosh, D. and McDonald, R.A.
#'(2010). \code{ Point Transect Sampling Along Linear Features } Biometrics ,
#'vol 66, no. 4, pp. 1247-1255.
#'@export
#'@keywords misc
#'@examples
#'
#'##Don't run:
#'#Example 1
#'norm.fit <-   nupoint.fit(c(50,20,50),
#'                       sight.x=krill$x,
#'                       sight.y=krill$y,
#'                       w=100,
#'                       theta.max=pi/3,
#'                       grad.type="NORM",
#'                       det.type="HNORM",
#'                       lower.b=c(1,1,1),
#'                       upper.b=c(120,100,100))
#'# ---------------------------
#'#   Parallel density gradient likelihood settings
#'# ---------------------------
#'#   Depth preference parametric form: NORM
#'# range detection function, g(r), parametric form: HNORM
#'# Parameter starting values = 50 20 50
#'# Angular detection function applied: FALSE
#'# Seabed attenuation function applied: FALSE
#'# Truncation distance = 100
#'# Swath (max) angle, rad = 1.047198
#'# Numerical integration grid density= 100
#'# ---------------------------
#'#   Estimating parameters
#'# ---------------------------
#'#   Maximum likelihood results
#'# ---------------------------
#'#   parameter estimates
#'# 81.736 25.144 38.194
#'# AIC = -3.02
#'# ---------------------------
#'#Example 2
#'#create an attenuation function:
#'at.fit=make.atten.f(depths=krill$z,form='norm',starting.pars=c(90,10),plot=TRUE)
#'#fit a Beta vertical distribution:
#'mbe.beta=nupoint.fit(pars=c(3,1.5,40), sight.x=krill$x, sight.y=krill$y,sight.z=krill$z,
#'                      w=100,theta.max=pi/3,grad.type='BETA',det.type='HNORM',grid.density=100,lower.b=c(0.5,0.5,10),upper.b=c(20,20,1000))
#'#---------------------------
#'#Parallel density gradient likelihood settings
#'#---------------------------
#'#Depth preference parametric form: BETA
#'#range detection function, g(r), parametric form: HNORM
#'#Parameter starting values = 3 1.5 40
#'#Angular detection function applied: FALSE
#'#Seabed attenuation function applied: TRUE
#'#function(seabedDepth) 1-pnorm(seabedDepth,95.1141535839473,10.5625798419803)
#'#<environment: 0x06e6beb0>
#'#Truncation distance = 100
#'#Swath (max) angle, rad = 1.047198
#'#Numerical integration grid density= 100
#'#---------------------------
#'#Estimating parameters
#'#---------------------------
#'#Maximum likelihood results
#'#---------------------------
#'#parameter estimates
#'#3.52 1.281 38.914
#'#AIC = 484.63
#'#---------------------------
#'##End don't run
#'
nupoint.fit <- function(pars,sight.x,sight.y,sight.z=NULL, w, theta.max,grad.type, det.type, n=NULL,
                   grid.density=100, angularDetect=FALSE,verbose=FALSE,lower.b,upper.b,optim.control=NULL){
  #### THIS FUNCTION IS DOCUMENTED IN NUPOINT #####
  #20121107: uses optim to obtain MLE for multibeam type problems.
  #INPUTS:
  #pars: 'a single vector of starting parameters to be estimated.  Parameters must be ordered: environment density gradient first, followed by detection function. See environ.fit.f'
  #sight.x: 'observation vector:  target x coordinate. This is cross-track distance for multi-beam echosounder(MBE) observations.'
  #sight.y: 'observation vector:  target y coordinate. This is water depth to the centre of a shoal for MBE observations.'
  #sight.z: 'observation vector:  default NULL if no attenuation function is required (all transects extend to w). For MBE observations, this is water depth to the centre of a shoal.'
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
  if(length(sight.z)>0){
    if(!('atten.f'%in% ls(envir=.GlobalEnv))){
      warning('perpendicular distance attenuation function missing when seabed at observation data passed to nupoint.fit in the sight.z ARG.')
      cat('WARNING: perpendicular distance attenuation function atten.f not in workspace.\n')
      cat('Either run make.atten.f to create atten.f or remove perpendicular distance (set ARG sight.z=NULL) from nupoint.fit (sight.z=NULL)\n')
      return('no_attenuation_function_atten.f_in_workspace')
    }
  if(length(sight.z)!=length(sight.y) ){
	warning('Observation vector (ARGS sight.z and sight.y) length mismatch')
	return('Observation_vector_length_mismatch')}
  }
  #observation vector length mistmatch
  if(length(sight.x)!=length(sight.y)) {
	warning('Observation vector (ARGS sight.x and sight.y) length mismatch')
	return('Observation_vector_length_mismatch')}
  
  
  #no angular detection function in workspace
  if(angularDetect){
    if(!('angDetF'%in% ls(envir=.GlobalEnv))){
      warning('Angular detection function, angDetF, missing when angularDetect=TRUE in nupoint.fit')
      cat('WARNING: Angular detection function, angDetF, not in workspace.')
      cat('Either create angDetF or set angularDetect=FALSE in nupoint.fit','\n')
      return('angular_detection_function_missing')
    }
  }
  if(any(sight.y>w)){
    warning('One or one observations exceed the truncation distance, w.')
    return('Observation(s)_exceed_w')
  }
  bound.chk=bound.chk.f(pars,lower=lower.b,upper=upper.b,f.n='nupoint.fit')
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
  cat('Parallel density gradient likelihood settings','\n')
  cat('---------------------------','\n')
  cat('Depth preference parametric form:',grad.type,'\n')
  cat('range detection function, g(r), parametric form:',det.type,'\n')
  cat('Parameter starting values =',pars,'\n')
  cat('Angular detection function applied:',angularDetect,'\n')
  if(angularDetect) print(angDetF)
  cat('Seabed attenuation function applied:',ifelse(length(sight.z)>0,'TRUE','FALSE'),'\n')
  if(length(sight.z)>0) print(atten.f)
  cat('Truncation distance =',w,'\n')
  cat('Swath (max) angle, rad =',theta.max,'\n')
  cat('Numerical integration grid density=',grid.density,'\n')
  cat('---------------------------','\n')
  
  xV= seq(0,w*sin(theta.max),length=grid.density)+(0.5*w/grid.density)  
  yV= seq(0,w,length=grid.density)+(0.5*w/grid.density)  
  xV=xV[-length(xV)]; yV=yV[-length(yV)]
  
  mask = outer(xV,yV,swathInOutF,w,theta=theta.max)
  grad.likhood=NULL
  cat('Estimating parameters \n')
  grad.likhood <- optim(pars, shore.only,sight.x=sight.x,swarm.z=sight.y, grid.density=grid.density,td=w,sb.z=sight.z,
                        g.type=grad.type,det.type=det.type,n=n,verbose=verbose,
                        mask=mask,xV=xV,yV=yV,angularDetect=angularDetect,
                        hessian=TRUE,method="L-BFGS-B", lower=lower.b, upper=upper.b,
                        control=control.list)
  grad.likhood$hessian=-grad.likhood$hessian
  if(grad.likhood$convergence!=0)
    warning('nupoint.fit: convergence error code returned by optim')
  if(any(grad.likhood$par < 1.01*pars  & grad.likhood$par > 0.99*pars ))
    warning('nupoint.fit: one or more parameter estimates within 1 percent of starting value.  Check convergence.')
  if(any(floor(grad.likhood$par)==lower.b))
    warning('nupoint.fit: one or more parameter estimates at lower bound.  Check convergence.')
  if(any(floor(grad.likhood$par)==upper.b))
    warning('nupoint.fit: one or more parameter estimates at upper bound.  Check convergence.')
  
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
