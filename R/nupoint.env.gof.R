#'Chi-squared goodness-of-fit test for environmental data
#'
#'This function calculates the one-dimensional chi-squared goodness-of-fit
#'statistic for environment preference data.
#'
#'The \code{pars} argument is the vector of parameter estimates obtained from
#'\code{\link{nupoint.env.fit}}. The arguments \code{grad.type},
#'\code{det.type} and \code{n} must ,match those used in
#'\code{\link{nupoint.env.fit}} to obtain estimates of \code{pars}.
#'
#'@param pars \code{vector of parameters estimated using nupoint.env.fit.}
#'@param r.mat \code{matrix:radial distances from observer to each regular grid
#'point in the survey region.}
#'@param z.mat \code{matrix:regular grid of the environmental feature
#'throughout the survey region.}
#'@param minz \code{minimum environmental feature value.}
#'@param wz \code{truncation distance in the z-dimension (environmental
#'feature).}
#'@param z.obs \code{observation vector:environmental feature.}
#'@param grad.type \code{parametric form of the environmental gradient (this is
#'the function describing habitat preference).  See nupoint.env.fit.}
#'@param det.type \code{parametric form of the detection function.  See
#'detectF.}
#'@param intervals \code{Either an integer number of bins for chi-squared test
#'(which will be used to generate bins of equal width).  Alternatively, bins of
#'unequal can be specified using vector break points of length 1+ the number of
#'desired bins.}
#'@param plot \code{logical default=FALSE. If TRUE a four panel diagnostic plot
#'is displayed.}
#'@param n \code{Default NULL. Use only when grad.type=MNORM to specify an
#'integer number of distributions in the normal mixture option for the
#'parametric form of the environmental gradient.  See nupoint.env.fit.}
#'@param dzdy.mat \code{matrix:regular grid of the rate of change of
#'environmental feature with respect to the y-dimension.}
#'@param breaks \code{vector of break points to calculate chi-squared
#'goodness-of-fit.  The argument nbr.intervals=NULL when breaks is used.}
#'@return List object of length two. 1) Chi-squared p-value, and 2)
#'goodness-of-fit table with columns: bin minimum, bin maximum, expected
#'sightings, observed sightings, Chisq
#'@seealso \code{\link{nupoint.env.fit}}
#'@export
#'@keywords misc
#'@examples
#'
#'##Don't run:
#'attach(sightings) # subset of beaked whale data from Arranz (submitted)
#'#fit a normal environmental feature preference and half-normal detection function:
#'environ.fit=nupoint.env.fit(pars=c(1000,200,3000),
#'              z=sighting.mat$z, rd=sighting.mat$r, dzdy=sighting.mat$dzdy,
#'        z.mat=z.mat, dzdy.mat=zGradmat, rd.mat=rd.mat,
#'              minz=minz, wx=wx, wy=wy, wz=wz,
#'              grad.type="NORM", det.type="HNORM",
#'              lower.b=c(-2000,1,1),upper.b=c(10000,10000,10000))
#'#after running nupoint.env.fit using the example sightings data with 'grad.type="NORM", det.type="HNORM"':
#'tt=nupoint.env.gof(pars=environ.fit$par, r.mat=rd.mat, z.mat=z.mat, minz=minz,wz=wz,
#'        z.obs=sighting.mat$z, grad.type="NORM", det.type="HNORM",
#'        intervals=13,plot=FALSE,dzdy.mat=zGradmat)
#'#-------------------------------------------------------
#'#1D Chi-squared Goodness-of-Fit results
#'#-------------------------------------------------------
#'#   bin.min bin.max    mids expected observed Chisq
#'#1   103.22  241.73  172.47     5.08        7  0.73
#'#2   241.73  380.25  310.99     8.81        6  0.89
#'#3   380.25  518.77  449.51    15.43       16  0.02
#'#4   518.77  657.29  588.03    22.63       17  1.40
#'#5   657.29  795.81  726.55    24.77       24  0.02
#'#6   795.81  934.32  865.07    28.88       29  0.00
#'#7   934.32 1072.84 1003.58    25.65       30  0.74
#'#8  1072.84 1211.36 1142.10    21.33       24  0.33
#'#9  1211.36 1349.88 1280.62    16.85       19  0.28
#'#10 1349.88 1488.40 1419.14    11.22        9  0.44
#'#11 1488.40 1626.91 1557.66     6.46        7  0.04
#'#12 1626.91 1765.43 1696.17     2.93        0  2.93
#'#13 1765.43 1903.95 1834.69     0.97        2  1.10
#'#-------------------------------------------------------
#'#Chi-squ. statistic   =  8.926069
#'#Number of parameters =  3
#'#Chi-squ. df          =  9
#'#Chi-squ. GoF p-value =  0.4441265
#'#-------------------------------------------------------
#'#
#'#example using user defined bin break points
#'binV=seq(sightings$minz,sightings$wz,length=14) #GoF bin break point vector.
#'binV=binV[-5] #remove break 5
#'tt=nupoint.env.gof(pars=environ.fit$par, r.mat=rd.mat, z.mat=z.mat, minz=minz,wz=wz,
#'        z.obs=sighting.mat$z, grad.type="NORM", det.type="HNORM",
#'        intervals=binV,plot=FALSE,dzdy.mat=zGradmat)
#'#-------------------------------------------------------
#'#1D Chi-squared Goodness-of-Fit results
#'#-------------------------------------------------------
#'#   bin.min bin.max    mids expected observed Chisq
#'#1   103.22  241.73  172.47     5.08        7  0.73
#'#2   241.73  380.25  310.99     8.81        6  0.89
#'#3   380.25  518.77  449.51    15.43       16  0.02
#'#4   518.77  795.81  657.29    47.40       41  0.86
#'#5   795.81  934.32  865.07    28.88       29  0.00
#'#6   934.32 1072.84 1003.58    25.65       30  0.74
#'#7  1072.84 1211.36 1142.10    21.33       24  0.33
#'#8  1211.36 1349.88 1280.62    16.85       19  0.28
#'#9  1349.88 1488.40 1419.14    11.22        9  0.44
#'#10 1488.40 1626.91 1557.66     6.46        7  0.04
#'#11 1626.91 1765.43 1696.17     2.93        0  2.93
#'#12 1765.43 1903.95 1834.69     0.97        2  1.10
#'#-------------------------------------------------------
#'#Chi-squ. statistic   =  8.366981
#'#Number of parameters =  3
#'#Chi-squ. df          =  8
#'#Chi-squ. GoF p-value =  0.3984676
#'#-------------------------------------------------------
#'detach(sightings)
#'##End don't run
#'
nupoint.env.gof <- function(pars,r.mat,z.mat,minz,wz,z.obs,grad.type,det.type,intervals,plot=FALSE,
              n=NULL,dzdy.mat,breaks=NULL)
{
  #20121031 nupoint.env.gof: chi-sq GoF for whale data (this is the 1D version)
  #INPUTS:    pars        =   model parameters, paramter vector optimised using optim. e.g.: 
  #                           pars[1]=depth gradient par 1;
  #                           pars[2]=depth gradient par 2;
  #                           pars[3]=half normal rng detection function. NB pars unpacked using
  #                           par.unpack.F        
  #           r.mat       =   radial distance from observer to centre of each integration grid cell.
  #           z.mat       =   seabed depth at each integration grid.
  #           minz        =   minimum z (depth) NB assumes z is positive.
  #           wz          =   seabed depth truncation distance. NB assumes z is positive
  #           z.obs       =   vector of depths at each observed cue,
  #           grad.type   =  cue depth disribution type e.g. "NORM"  - this is currently coded for
  #                           a radial detection function only
  #           det.type  = detection function form.
  #           intervals = number of intervals or vector of break points for use in GoF.
  #           n           = number of distributions in a multinomial distribution (default NULL).
  #           plot        = TRUE/FALSE 4 panel plot of GoF
  #           dzdy.mat    = matrix of derivatives of z with respct to y
  #REQUIRES:  fields package, functions: histline; detectF; pi.z.f
  #RETURNS:   single value of the log-likelihood numerator.
  #CALLS:     par.unpack.F; pi.z.f; package fields for plotting.
  parList=par.unpack.F(grad.type,pars,n)
  grad.pars=parList[[1]]; det.par=parList[[2]]
  
  detV=detectF(det.type=det.type,det.par=det.par,rd=r.mat)
  
  pi.z=pi.z.f(g.type=grad.type,pars=grad.pars,z=z.mat,z.lim=c(minz,wz),n=n)
  p.pi.z=detV*pi.z*as.vector(dzdy.mat)
  
  #20130321: modified to include fixed intervals or breaks for the GoF bins:
  if(!is.numeric(intervals)) stop('intervals Arg must be numeric.')
  if(length(intervals)==1){
	#bit of error trapping:
	if((intervals-floor(intervals))!=0)
		stop('intervals ARG has length = 1, so must be an integer specifying the number of intervals.')
	#instance where user specifies number of intervals:
	cutV=seq(minz,wz,length=intervals+1) #vector of lower and upper bounds in y-dimension
  } else {
	cutV=intervals
  }
  #check interval range
  if(min(cutV)>minz)
	stop('Minimum interval specified in intervals ARG > minz ARG.')
  if(max(cutV)<wz)
	stop('Maximum interval specified in intervals ARG < wz ARG.')

  
  
  #if(is.null(breaks)) cutV=seq(floor(minz),ceiling(wz),length.out=(nbr.intervals+1)) #determine seabed intervals
  #else {
  #  cutV=breaks;ncut=length(cutV)
  #  if(breaks[1]>floor(minz)) cutV[1]=floor(minz)
  #  if(breaks[ncut]<ceiling(wz)) cutV[ncut]=ceiling(wz)
  #}
  
  
  mids=cutV[-length(cutV)]+diff(cutV)/2 #mid-point of seabed intervals.
  cut.mat=matrix(cut(z.mat,cutV,labels=FALSE),nrow=nrow(z.mat)) #determine a matrix of seabed intervals.
  expectedV=rep(0,length(mids)) #create an empty matrix of expected number of sightings.
  ppn=tapply(as.vector(p.pi.z),as.vector(cut.mat),sum,na.rm=TRUE)
  ppn=ppn/sum(ppn)
  expected=ppn*length(z.obs)
  expectedV[as.numeric(names(expected))]=expected
  observedV=rep(0,length(mids))
  observedTAB=tapply(z.obs,cut(z.obs,cutV,labels=FALSE),length)
  observedV[as.numeric(names(observedTAB))]=observedTAB
  out=data.frame(bin.min=cutV[-length(cutV)],bin.max=cutV[-1],mids=mids,
                 expected=expectedV,observed=observedV,Chisq=(observedV-expectedV)^2/expectedV)
  cat('-------------------------------------------------------\n')
  cat('1D Chi-squared Goodness-of-Fit results\n')
  cat('-------------------------------------------------------\n')
  print(round(out,2))
  Chisq=sum((out$observed-out$expected)**2/out$expected,na.rm=TRUE)
  npar=length(pars)
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
  
  if(plot){
    require(fields,quietly = TRUE)
    par(mfrow=c(2,2),mar=c(3,2,2,2))
    image.plot(z.mat,main='z(x,y)',cex.main=0.75)
    image.plot(cut.mat,main='interval',cex.main=0.75)
    image.plot(100*p.pi.z/sum(p.pi.z,na.rm=TRUE),main='expected % of detections',cex.main=0.75)
    
    ylim=c(0,max(out$expected,out$observed))
    histline(height=out$expected,breaks=cutV,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=ylim,xlab="Seabed depth, m",ylab="Frequency",main="Observed (red) & Expected (black)",cex.main=0.75)
    histline(height=out$observed,breaks=cutV,lineonly=TRUE,outline=TRUE,fill=FALSE,col="red")
    if(Chisq.df>0) mtext(paste('chi-squ. p-value=',  round(chiSqp,2)),cex=0.5)
  }
  return(list(p.value=chiSqp,Chisq.table=out))
}  #end GoFenvironF function
