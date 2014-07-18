#'Non-parametric bootstrap for environmental preference
#'
#'Calculate the variance of parameter estimates for the environmental
#'preference gradient using a non-parameteric bootstrap.
#'
#'The \code{\link{sightings}} argument is a list object which must have the
#'following structure:
#'
#'\describe{ \item{list("sighting.mat")}{whale cue sightings;}
#'\item{list("x.mat")}{survey region description: x-coordinates;}
#'\item{list("y.mat")}{survey region description: y-coordinates;}
#'\item{list("rd.mat")}{survey region description: radial distance;}
#'\item{list("z.mat")}{survey region description: environment variable;}
#'\item{list("zGradmat")}{survey region description: change in environment
#'variable;} \item{list("x")}{survey region description: x-coordinate vector;}
#'\item{list("y")}{survey region description: y-coordinate vector;}
#'\item{list("obsx")}{observer x-coordinate;} \item{list("obsy")}{observer
#'y-coordinate;} \item{list("wx")}{x-dimension truncation distance;}
#'\item{list("wy")}{y-dimension truncation distance;}
#'\item{list("wz")}{z-dimension truncation distance;}
#'\item{list("minz")}{Minimum z dimension value;} }
#'
#'The \code{sighting.mat} object within \code{\link{sightings}} contains the
#'following data for each observation
#'
#'\describe{ \item{list("x")}{x-coordinate of sighting;}
#'\item{list("y")}{y-coordinate of sighting;} \item{list("r")}{radial distance
#'from observer to sighting;} \item{list("z")}{seabed depth (or other
#'environmental variable of interest) at each sighting;}
#'\item{list("dzdy")}{seabed depth (or other environmental variable of
#'interest) gradient with respect to the y-dimension at each sighting, and}
#'\item{list("obs.period")}{observation period. This variable is used as the
#'grouping parameter in \code{nupoint.env.boot} and is specified in the
#'\code{blockVar} argument.} }
#'
#'Example sighting data \code{\link{sightings}} is included with the package.
#'
#'@param sightings \code{list object of sightings and survey region variables
#'(see details and the }env.pdf' vignette).'
#'@param nboot \code{number of bootstraps.}
#'@param blockVar \code{string column name or integer column number in the
#'sightings$sighting.mat data frame which is the bootstrap block variable.}
#'@param blockVar \code{string column name or integer column number in the
#'sightings$sighting.mat data frame which is the bootstrap block variable.}
#'@param initial.pars \code{intitial parameters for use by the nupoint.env.fit
#'function (see the pars argument in nupoint.env.fit).}
#'@param grad.type \code{environmental gradient type (See nupoint.env.fit).}
#'@param det.type \code{detection function form. (See detecF and
#'nupoint.env.fit).}
#'@param n \code{number of distributions in a mixture. (See nupoint.env.fit).}
#'@param lower.b \code{lower parameter space bound for nupoint.env.fit
#'function.}
#'@param upper.b \code{upper parameter space bound for nupoint.env.fit
#'function.}
#'@return 2D array nrow = nboot; ncol = initial.pars, each row contains
#'parameter estimates for each bootstrap.
#'@seealso \code{\link{detectF}}, \code{\link{nupoint.env.fit}},
#'\code{\link{sightings}}
#'@export
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'boot.res=nupoint.env.boot(sightings=sightings,
#'               nboot=100,
#'               blockVar='obs.period',
#'               initial.pars=c(1171.479,503.87,3020.442),
#'               grad.type='NORM',
#'               det.type='HNORM',
#'               lower.b=c(-2000,1,1),
#'               upper.b=c(10000,10000,10000))
#'par(mfrow=c(1,3))
#'sapply(colnames(boot.res),function(x) hist(boot.res[,x],main=x,xlab='Estimate'))
#'}
#'
nupoint.env.boot <- function (sightings,nboot,blockVar=NULL,initial.pars,grad.type,det.type,n=NULL,lower.b,upper.b) 
{
  #find grouping variable vector:
  if(!is.integer(blockVar) & !is.character(blockVar))  {warning('blockVar argument is not of type integer or character')
      return('error:nupoint.env.boot')
  } 
  if(is.character(blockVar)){
	blockVar=which(names(sightings$sighting.mat)==blockVar)
    groupVar=sightings$sighting.mat[,blockVar]}
  if(is.integer(blockVar))
    { 
      if(blockVar>ncol(sightings$sighting.mat)){warning('Integer blockVar argument exceeds sightings matrix ncol')
        return('error:nupoint.env.boot')
      } #end error trap for integer blockVar
    groupVar=unlist(sightings$sighting.mat[,blockVar])
    }
                 
  #end of groupVar handling and error check
  unique.grp=unique(groupVar)
  cat('---------------------------------------\n')
  cat('Non-parameteric bootstrap \n')
  cat('Grouping variable levels and number of observations per level \n')
  print(table(groupVar))
  cat('---------------------------------------\n')
  
  required.nbr.obs=nrow(sightings$sighting.mat) #minimum required number of observations in each bootstrap
  res.mat=matrix(NA,nrow=nboot,ncol=length(initial.pars))
  dimnames(res.mat)=list(1:nboot,paste('par_hat',1:length(initial.pars),sep='_'))
  
  pb <- txtProgressBar(min = 0, max = nboot, style = 3) #progress bar
  #end of variable initialisation

  #start boostrap:
  for(i in 1:nboot) {
    setTxtProgressBar(pb, i) #set progress bar at the start of each iteration 
    ##create current bootstrap data
    nbr.c.obs=0 #number of observations in the current bootstrap data set.
    while(nbr.c.obs<required.nbr.obs)
    {
      TMP.obs=sightings$sighting.mat[which(sightings$sighting.mat[,blockVar]==unique.grp[sample(1:length(unique.grp),1)]),]
      if(nbr.c.obs==0){ sight.mat=TMP.obs } else {sight.mat=rbind(sight.mat,TMP.obs)}
      nbr.c.obs=nrow(sight.mat)
    }
    ##fit model
    c.fit=0
    c.fit=try(nupoint.env.fit(pars=initial.pars,
                            z=sight.mat$z, #use the resampled obs
                            rd=sight.mat$r,#use the resampled obs
                            dzdy=sight.mat$dzdy,#use the resampled obs
                            z.mat=sightings$z.mat,
                            dzdy.mat=sightings$zGradmat,
                            rd.mat=sightings$rd.mat,
                            minz=sightings$minz,
                            wx=sightings$wx,
                            wy=sightings$wy,
                            wz=sightings$wz,
                            grad.type=grad.type,
                            det.type=det.type,
                            n=n,
                            lower.b=lower.b,upper.b=upper.b) )
    
  if(length(c.fit)==7) { #check if optim converged
  res.mat[i,]=c.fit$par  
  } else {
    warning(paste('nupoint.env.fit failed to converge during bootstrap',i))} #end of if loop (optim convergence check )
  } #end of boostrap
  close(pb) #close progress bar
  return(res.mat)
}
