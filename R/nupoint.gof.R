#'Chi-squared goodness-of-fit test for parallel density gradients
#'
#'This function calculates the one-dimensional chi-squared goodness-of-fit
#'statistic for environment preference data.
#'
#'The \code{pars} argument is the vector of parameter estimates obtained from
#'\code{\link{nupoint.fit}}. The arguments \code{grad.type}, \code{det.type}
#'and \code{n} must ,match those used in \code{\link{nupoint.fit}} to obtain
#'estimates of \code{pars}.
#'
#'@param y.obs \code{observation vector: target y coordinate. This is water
#'depth at the centre of an aggregation for multi-beam echosounder
#'observations.}
#'@param pars \code{vector of parameters estimated using nupoint.fit.}
#'@param w \code{truncation distance in the y-dimension (depth for multi-beam
#'echosounder observations)}
#'@param grad.type \code{parametric form of the target object density function.
#'See nupoint.fit and nupoint.env.fit.}
#'@param det.type \code{parametric form of the detection function.  See
#'detectF.}
#'@param n \code{Default NULL. Use only when grad.type=MNORM to specify an
#'integer number of distributions in the normal mixture option for the
#'parametric form of the environmental gradient.  See nupoint.fit.}
#'@param intervals \code{Either an integer number of bins for chi-squared test
#'(which will be used to generate bins of equal width).  Alternatively, bins of
#'unequal can be specified using vector break points of length 1+ the number of
#'desired bins.}
#'@param verbose \code{return summary text of goodness-of-fit to the R
#'console.}
#'@param plot \code{logical default=FALSE. If TRUE a two panel diagnostic plot
#'is displayed.}
#'@param attenuation \code{Logical.  TRUE = seabed attenuation function is
#'specified; FALSE = no seabed attenuation function.}
#'@param angularDetect \code{logical default=FALSE. Should angular detection be
#'considered (not yet implemented).}
#'@param xmax \code{maximum x-dimension distance. Default is NULL, xmax = w *
#'sin(theta).}
#'@param theta.max \code{maximum observation angle radians, (swath width in the
#'multi-beam echosounder case).  Maximum angle is pi/2.}
#'@param grid.resolution \code{numerical integration grid density (number of
#'elements in each (x and y) dimensions).}
#'@return List object of length two. 1) Chi-squared p-value, and 2)
#'goodness-of-fit table with columns: bin minimum, bin maximum, expected
#'sightings, observed sightings, Chisq.
#'@seealso \code{\link{nupoint.fit}}, \code{\link{detectF}},
#'\code{\link{nupoint.env.fit}}
#'@export
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'# norm.gof=nupoint.gof(y.obs=krill$y,
#'#                pars=c(81.736 ,25.144 ,38.194),
#'#                w=100,
#'#                grad.type='NORM',
#'#                det.type='HNORM',
#'#                n=NULL,
#'#                intervals=11,
#'#                angularDetect=FALSE,
#'#                theta.max=pi/3,
#'#                grid.resolution=100)
#'# -------------------------------------------------------
#'#   nupoint: 1D Chi-squared Goodness-of-Fit results
#'# -------------------------------------------------------
#'#   bin.min bin.max  mids expected observed Chisq
#'# 1    10.74   18.69 14.71     0.95        2  1.15
#'# 2    18.69   26.63 22.66     3.96        4  0.00
#'# 3    26.63   34.58 30.60    10.67        9  0.26
#'# 4    34.58   42.52 38.55    21.39       26  0.99
#'# 5    42.52   50.47 46.50    28.88       28  0.03
#'# 6    50.47   58.42 54.44    41.75       34  1.44
#'# 7    58.42   66.36 62.39    44.67       41  0.30
#'# 8    66.36   74.31 70.34    40.51       46  0.74
#'# 9    74.31   82.26 78.28    30.78       41  3.39
#'# 10   82.26   90.20 86.23    19.08       16  0.50
#'# 11   90.20   98.15 94.17     8.35        4  2.27
#'# -------------------------------------------------------
#'#   Chi-squ. statistic   =  11.06844
#'# Number of parameters =  3
#'# Chi-squ. df          =  7
#'# Chi-squ. GoF p-value =  0.1356615
#'# -------------------------------------------------------
#'}
#'
nupoint.gof <- function(y.obs,pars,w,grad.type,det.type,n=NULL,intervals,verbose=TRUE,plot=FALSE,attenuation=FALSE,
                        angularDetect=FALSE,xmax=NULL,theta.max,grid.resolution=100)
{
  if(length(xmax)==0) xmax=w*sin(theta.max)
  #20130321: modified to include fixed intervals or breaks for the GoF bins:
  if(!is.numeric(intervals)) stop('intervals Arg must be numeric.')
  if(length(intervals)==1){
	#bit of error trapping:
	if((intervals-floor(intervals))!=0)
		stop('intervals ARG has length = 1, so must be an integer specifying the number of intervals.')
	#instance where user specifies number of intervals:
	breaks=seq(min(y.obs),max(y.obs),length=intervals+1) #vector of lower and upper bounds in y-dimension
  } else {
	breaks=intervals
  }
  #check interval range
  if(min(breaks)>min(y.obs))
	stop('Minimum interval specified in intervals ARG > minimum value in y.obs ARG.')
  if(max(breaks)<max(y.obs))
	stop('Maximum interval specified in intervals ARG < maximum value in y.obs ARG.')

  require(fields,quietly=TRUE)
  pars=par.unpack.F(grad.type,pars,n)
  xV=seq(0,xmax,length.out=grid.resolution)
  yV=seq(0,w,length.out=grid.resolution)
  mask= outer(xV,yV,swathInOutF,w=w,theta=theta.max)
  temp=outer(xV,yV,shore.calc.integral.grid,pars, trunc.dist=w,grad.type,det.type,n,attenuation=attenuation,angularDetect) *(1/xmax) *mask
  z.mat=matrix(rep(yV,length(xV)),nrow(temp),byrow=T) #matrix of depths aligned with temp matrix
  cut.mat=matrix(cut(z.mat,breaks,labels=FALSE),nrow=nrow(temp)) #determine a matrix of depth intervals in y-dimension.
  ppn=tapply(as.vector(temp),as.vector(cut.mat),sum,na.rm=TRUE) #integrate each depth
  ppn=ppn/sum(ppn)
  eV=length(y.obs) * ppn
  histObj=hist(y.obs,breaks,plot=F)
  obsV=histObj$counts
  
  #graphics:
  if(plot){
    par(mfrow=c(2,1))
    image.plot(x=xV,y=yV,z=100*temp/sum(temp,na.rm=T),
               xlab='x-coordinate',ylab='y-coordinate',xlim=c(0,100),
               main='Expected % of detections',
               ylim=c(100,0))
    histline(height=eV,breaks,ylim=c(0,max(c(obsV,eV))),xlab='y-coordinate',ylab='Frequency')
    histline(height=obsV,breaks,lineonly=TRUE,outline=TRUE,fill=FALSE,col="red")
  }
  
  #GoF calcs:
  out=data.frame(bin.min=breaks[-length(breaks)],bin.max=breaks[-1],mids=histObj$mids,
                 expected=eV,observed=obsV,Chisq=(obsV-eV)^2/eV)
  Chisq=sum((out$observed-out$expected)**2/out$expected,na.rm=TRUE)
  npar=length(unlist(pars))
  Chisq.df=nrow(out)-npar-1
  chiSqp=1-pchisq(Chisq,Chisq.df)
  if(Chisq.df>0 & plot==T) mtext(paste('chi-squ. p-value=',  round(chiSqp,2)),cex=1)
  if(verbose){
    cat('-------------------------------------------------------\n')
    cat('nupoint: 1D Chi-squared Goodness-of-Fit results\n')
    cat('-------------------------------------------------------\n')
    print(round(out,2))
    Chisq=sum((out$observed-out$expected)**2/out$expected,na.rm=TRUE)
    Chisq.df=nrow(out)-npar-1
    cat('-------------------------------------------------------\n')
    cat('Chi-squ. statistic   = ',Chisq,'\n')
    cat('Number of parameters = ',npar,'\n')
    cat('Chi-squ. df          = ',Chisq.df,'\n')
    chiSqp=NA
    if(Chisq.df>0) {
      chiSqp=1-pchisq(Chisq,Chisq.df)
      cat('Chi-squ. GoF p-value = ',chiSqp,'\n')
    }else {
      cat("Degrees of freedom <1 so can't get p-value\n")
    }
    cat('-------------------------------------------------------\n')
  }  
  return(list(chiTable=out,ChiSqup=chiSqp))
}
