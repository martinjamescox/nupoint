#'Fit a attenuation function
#'
#'Fit an attenuation function calls the \code{mle} function from the
#'\pkg{FAmle} package. The seabed attenuation function accounts for variation
#'in maximum transect observation distance. In terms of multi-beam echosounder
#'surveys, the attenuation function accounts for varying seabed depth. See Cox
#'et al. (2011).
#'
#'This function must be run prior to \code{\link{nupoint.fit}} if an
#'attenuation function is required in the likelihood fitted in
#'\code{\link{nupoint.fit}}. The name of the \code{atten.f} must not be
#'changed.
#'
#'@param prep.dist \code{vector of a representative sample of prependicular
#'distances(e.g. in the multibeam context, seabed depth under a krill swarm).}
#'@param form \code{any standard R distribution (see Distributions) that can be
#'acommodated by the FAmle package e.g. norm for normal.}
#'@param starting.parameters \code{vector of starting parameters describing the
#'distribution specified in form.}
#'@param plot \code{logical. produce a plot of the attenuation function}
#'@param summary \code{default=FALSE. If TRUE, print a summary of the
#'attenuation function fit to the console.}
#'@return 1) creates a function in the R workspace called \code{atten.f} 2)
#'returns a \code{FAmle::mle} object.
#'@section Reference: Cox, M.J., Borchers, D.L., Demer, D.A., Cutter, G.R., and
#'Brierley, A.S. (2011). Estimating the density of Antarctic krill (Euphausia
#'superba) from multi-beam echo-sounder observations using distance sampling
#'methods. Journal of the Royal Statistical Society: Series C (Applied
#'Statistics).
#'@seealso \code{\link{nupoint.fit}}
#'@export
#'@keywords misc
#'@examples
#'
#'#use the seabed depth under each krill swarm 'krill$z' to create an attenuation function based on a normal distribution
#'\dontrun{
#'ls()
#'atten.fit=make.atten.f(prep.dist=krill$z,form='norm',starting.pars=c(90,10))
#'ls() #function named atten.f now exists in workspace
#'}
#'
make.atten.f <- function(prep.dist,form,starting.pars,plot=FALSE,summary=FALSE){ 
  #### THIS FUNCTION IS DOCUMENTED IN NUPOINT #####
  #20121107: make attenuation function for gradient with respect to linear feature
#   INPUTS:
#   prep.dist: 'vector of a representative sample of prep.dist (or other observation describing attenuation.'
#   form: 'any standard R distribution (see Distributions) that can be acoomodated by the FAmle package e.g. norm for  normal.'
#   starting.parameters: 'vector of starting parameters describing the distribution specified in form.'
#   plot: 'logical. produce a plot of the attenuation function'
#   summary: 'print a summary of the attenuation function fit to the console.'
#RETURNS: a FAmle summary object and an attenuation function is returned to the R workspace called atten.f
  #require(FAmle,quietly=TRUE)
  #mleObj=FAmle::mle(x=prep.dist,dist=form,start=starting.pars)
  mleObj=mle(x=prep.dist,dist=form,start=starting.pars)
  if(summary) print(mleObj$fit)
  functionString=paste('atten.f', "<<- function(seabedDepth) 1-p",form,'(seabedDepth,',
                       paste(mleObj$par.hat,collapse=','),')',sep='')
  eval(parse(text=functionString))
  cat('the following seabed attenuation function has been returned to the R workspace:\n')
  cat(functionString,'\n')
  if(plot)
  {
    plot(rev(sort(prep.dist)),1:length(prep.dist)/length(prep.dist),xlab='seabed depth',
         ylab='Attenuation effect',type='l',ylim=c(0,1))
    x=seq(min(prep.dist),max(prep.dist),1)
    lines(x,atten.f(x),col='grey',lwd=2)
    legend('topright',c('ECDF','Model fit'),col=c('black','grey'),lty=1,lwd=c(1,2))  
  }  
  return(mleObj)
}
