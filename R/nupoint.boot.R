#'Non-parametric bootstrap for depth preference
#'
#'Calculate the variance of parameter and density estimates for a parallel
#'density gradient using a non-parameteric bootstrap
#'
#'The \code{observations} argument is a data.frame object which must have the
#'following structure:
#'
#'\describe{ \item{list("transect")}{Transect number} \item{list("x")}{Shoal
#'cross-track distance} \item{list("y")}{Shoal depth} \item{list("r")}{Shoal
#'radial distance} \item{list("theta")}{Shoal angle, rad}
#'\item{list("bio.g")}{Shoal biomass} \item{list("z")}{Seabed depth under
#'shoal} }
#'
#'For each bootstrap, the vertical distrubtion and detection function
#'parameters are estimated along with P* (see \code{\link{nupoint.p.star.f}}).
#'If total line transect length, L, is specified then shoal volumetric density,
#'\code{shoal.vol.den.hat}, and shoal areal density \code{shoal.areal.den.hat}
#'are also estimated. If survey region area, A, is specified then shoal
#'abundance, \code{N.hat} is estimated. If the shoalBioMass argument is
#'specified, then shoal expected biomass, \code{E.bio} is estimated as is
#'survey region biomass, \code{E.bio area.biomass}.
#'
#'The \code{shoalBioMass} argument can be given as: NULL, is which case no
#'shoal or area biomass is estimated; a fixed mean shoal biomass value (the
#'same value is used in each bootstrap), or a linear model formula which will
#'be calculated for each bootstrap e.g. \code{log(bio.g) ~ theta}, see example.
#'
#'@param observations \code{data.frame of shaol observations.  See details and
#'}summary(krill)' object).
#'@param nboot \code{number of bootstraps.}
#'@param blockVar \code{string column name or integer column number in the
#'observations data frame which is the bootstrap block variable.}
#'@param initial.pars \code{intitial parameters for use by the nupoint.fit
#'function (see the pars argument in nupoint.fit).}
#'@param grad.type \code{parametric form of the environmental gradient.  See
#'environ.fit.f.}
#'@param det.type \code{parametric form of the detection function.  See
#'detectF.}
#'@param w \code{truncation distance in the y-dimension (depth for MBE)}
#'@param theta.max \code{maximum observation angle rad, (MBE swath width).
#'Maximum angle is pi/2.}
#'@param n \code{Default NULL. Use only when grad.type=MNORM to specify an
#'integer number of distributions in the normal mixture option for the
#'parametric form of the environmental gradient.  See details.}
#'@param attenuation \code{Logical.  TRUE = seabed attenuation function is
#'specified; FALSE = no seabed attenuation function.}
#'@param grid.density \code{numerical integration grid density (number of
#'elements in each (x and y) dimensions.}
#'@param lower.b \code{lower parameter space bound for nupoint.fit function.}
#'@param upper.b \code{upper parameter space bound for nupoint.fit function.}
#'@param angularDetect \code{logical:FALSE. Should angular detection be
#'considered (not yet implemented).}
#'@param L \code{Total line transect length}
#'@param A \code{Survey region area}
#'@param shoalBioMass \code{either a numeric value giving the mean shoal
#'biomass or a character string giving the regression equation, or NULL.  If
#'NULL then survery region biomass is not calculated.}
#'@return 2D array nrow = nboot; ncol = initial.pars, P*. Each row in the array
#'contains parameter estimates for each bootstrap. Plus, if argument L is
#'specified shoal volumetric and areal density. Plus, if argument A is
#'specified, number of shoals, N. Plus, if shoal biomass (\code{shoalBioMass}
#'argument) is specified, expected shoal biomass, E[b], and survey region
#'biomass. See details and the \code{parallel.pdf} vignette.
#'@seealso \code{\link{detectF}}, \code{\link{nupoint.fit}},'nupoint.p.star',
#'\code{\link{nupoint.vol.density}},'krill'
#'@export
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'boot.res <- nupoint.boot(observations=krill,
#'  nboot=50,
#'  blockVar='transect',
#'  initial.pars=c(81.736, 25.144, 38.194),
#'  grad.type='NORM',
#'  det.type='HNORM',
#'  w=100,
#'  theta.max=pi/3,
#'  n=NULL,
#'  attenuation=FALSE,
#'  grid.density=100,
#'  lower.b=c(1,1,1),
#'  upper.b=c(120,100,100),
#'  angularDetect=FALSE,
#'  L=11*2.5*1e3,
#'  A=3300000,
#'  shoalBioMass='log(bio.g) ~ theta')
#'}
#'
nupoint.boot <- function (observations, nboot, blockVar, 
                       initial.pars, grad.type, det.type,
                      w, theta.max, n=NULL, attenuation=FALSE,
                      grid.density=100, lower.b, upper.b,
                      angularDetect=FALSE,
                      L=NULL, A=NULL, shoalBioMass=NULL) 
{
  #20130201: function created.
  #Initialise variables:
  #find grouping variable vector:
  
  if(!is.integer(blockVar) & !is.character(blockVar))  {warning('blockVar argument is not of type integer or character')
                                                        return('error:nupoint.boot')
  }
  if(is.character(blockVar)){
    if(!blockVar %in% names(observations)) { warning('blockVar ARG column name not found. Exiting.')
      return('error:nupoint.boot')
    }
    blockVar=which(names(observations)==blockVar)
    }
  if(is.integer(blockVar))
  { 
    if(blockVar>ncol(observations)){warning('Integer blockVar argument exceeds observations matrix ncol. Exiting')
                                    return('error:nupoint.boot')
    } #end error trap for integer blockVar
    groupVar=unlist(observations[,blockVar])
  }
  #end of blockVar handling and error check
  unique.grp=unique(groupVar)
  cat('---------------------------------------\n')
  cat('Non-parameteric bootstrap \n')
  cat('Blocking variable levels and number of observations per level \n')
  print(table(groupVar))
  cat('---------------------------------------\n')
  #check if transect length is available
  
  additionalCols=0
  do.den=do.N=do.bio=FALSE #flags for additional calculations
  colNames=c(paste('par_hat',1:length(initial.pars),sep='_'),'P.star.hat')
  if(length(L)==0) {message('No total transect length, L ARG, in nupoint.boot FUN specified.  Shoal densities will not be estimated')
                    if(length(A)>0 & length(L)==0) {warning('Survey area specified, but not total transect length. Exiting') 
                                                    return('error:nupoint.boot: A ARG specified, but not L ARG')}
                    } else {
    additionalCols=2
    do.den=TRUE  
    colNames=c(colNames,'shoal.vol.den.hat','shoal.areal.den.hat')
    
    #check if the survery area is available
    if(length(A)==0) message('No survey area (survey area, A, argument in nupoint.boot FUN) specified.  Nbr of shoals will not be estimated')
    if(length(A)>0 & !is.numeric(A)) {warning('Non-numeric survey area specified. Exiting.') 
                                                        return('error:nupoint.boot: non-numeric survey area')}
    #check if biomass is available
    if(length(A)==0 & length(shoalBioMass)>0) {warning('shoalBioMass specified, but survey area, A, not specified. Exiting') 
                                                         return('error:nupoint.boot: shoalBioMass specified, but survey area not specified')}
    #end of error trapping
  } #end of L check else
 
  
  if(length(A)>0){
    additionalCols=additionalCols+1
    colNames=c(colNames,'N.hat')
    do.N=TRUE
  }
  if(length(shoalBioMass)>0)
  {
    bioFlag=ifelse(is.numeric(shoalBioMass),paste('A fixed value of',shoalBioMass,'biomass units'),paste('a linear regression of ',shoalBioMass))
    message(paste('nupoint.boot will esimate biomass using',bioFlag))
    additionalCols=additionalCols+2
    colNames=c(colNames,'E.bio','area.biomass')
    do.bio=TRUE
  }
  required.nbr.obs=nrow(observations)  #minimum required number of observations in each bootstrap
  res.mat=matrix(NA,nrow=nboot,ncol=length(initial.pars)+1+additionalCols) #fill in later
  dimnames(res.mat)=list(1:nboot,colNames)
  
  pb <- txtProgressBar(min = 0, max = nboot, style = 3) #progress bar
  #end of variable initialisation
  
  #start boostrap:
  for(i in 1:nboot) {
    setTxtProgressBar(pb, i) #set progress bar at the start of each iteration 
    ##create current bootstrap data
    nbr.c.obs=0 #number of observations in the current bootstrap data set.
    while(nbr.c.obs<required.nbr.obs)
    {
      TMP.obs=observations[which(observations[,blockVar]==unique.grp[sample(1:length(unique.grp),1)]),]
      if(nbr.c.obs==0){ sight.mat=TMP.obs } else {sight.mat=rbind(sight.mat,TMP.obs)}
      nbr.c.obs=nrow(sight.mat)
    }
    attenuation.Obj=NULL
    if(attenuation) attenuation.Obj=sight.mat$z
    ##fit model
    c.fit=0
    c.fit=try(nupoint.fit(pars=initial.pars,
                        sight.x=sight.mat$x,
                        sight.y=sight.mat$y,
                        sight.z=attenuation.Obj, 
                        w=w, 
                        theta.max=theta.max,
                        grad.type=grad.type, 
                        det.type=det.type, 
                        n=n,
                        grid.density=grid.density, 
                        angularDetect=angularDetect,
                        verbose=FALSE,
                        lower.b=lower.b,
                        upper.b=upper.b,
                        optim.control=NULL))

  if(length(c.fit)==7) { #check if optim converged
      #add estimated parameters, estimated P*, vol den and areal.den
      P.Star=nupoint.p.star(pars.hat=c.fit$par,grad.type=grad.type,det.type=det.type,n=n,w=w,theta.max=theta.max,grid.density=grid.density,
                          attenuation=attenuation,angularDetect=angularDetect)
      c.res=c(c.fit$par,P.Star)
      if(do.den){D3=nupoint.vol.density(n.seen=nrow(sight.mat),L=L,w=w,theta.max=theta.max,P.star=P.Star)
                D2=D3*w
                  c.res=c(c.res,D3,D2)}
      if(do.N) {N=A*D2
        c.res=c(c.res,N)}
      if(do.bio){
        if(is.character(shoalBioMass)){
        lmString=paste('b.lm <-lm(',shoalBioMass,',data=sight.mat)',sep='')
        eval(parse(text=lmString))
        e.bio=as.vector(exp(coefficients(b.lm)[1] + summary(b.lm)$coefficients[1,2]**2/2))} else
        {e.bio=shoalBioMass}
        svyBio=e.bio*N
        c.res=c(c.res,e.bio,svyBio)
        
      } #end do.bio
      res.mat[i,]=c.res  
    } else {
      warning(paste('nupoint.fit failed to converge during bootstrap',i))} #end of if loop (optim convergence check )
  } #end of boostrap
  close(pb) #close progress bar
  return(res.mat)
}
