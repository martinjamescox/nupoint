#' Chi-sq GoF for whale data (this is the 1D version)
#' @param    pars        =   model parameters, paramter vector optimised using optim. e.g.: 
  #                           pars[1]=depth gradient par 1;
  #                           pars[2]=depth gradient par 2;
  #                           pars[3]=half normal rng detection function. NB pars unpacked using
  #                           par.unpack.F        
#' @param r.mat  radial distance from observer to centre of each integration grid cell.
#' @param z.mat  seabed depth at each integration grid.
#' @param minz  minimum z (depth) NB assumes z is positive.
#' @param wz   seabed depth truncation distance. NB assumes z is positive
#' @param z.obs  vector of depths at each observed cue,
#' @param grad.type cue depth disribution type e.g. "NORM"  - this is currently coded for a radial detection function only
#' @param det.type detection function form.
#' @param nbr.intervals  number of intervals for use in GoF.
#' @param n  number of distributions in a multinomial distribution (default NULL).
#' @param plot TRUE/FALSE 4 panel plot of GoF
#' @param dzdy.mat matrix of derivatives of z with respct to y
#' @note  Requires fields package, functions: histline; detectF; pi.z.f
#' @return   single value of the log-likelihood numerator.
#' @details     Calls par.unpack.F; pi.z.f; package fields for plotting.
GoF.environ.fit.f <- function(pars,r.mat,z.mat,minz,wz,z.obs,grad.type,det.type,nbr.intervals,plot=FALSE,
                              n=NULL,dzdy.mat,breaks=NULL)
{
  parList=par.unpack.F(grad.type,pars,n)
  grad.pars=parList[[1]]; det.par=parList[[2]]
  
  detV=detectF(det.type=det.type,det.par=det.par,rd=r.mat)
  
  pi.z=pi.z.f(g.type=grad.type,pars=grad.pars,z=z.mat,z.lim=c(minz,wz),n=n)
  p.pi.z=detV*pi.z*as.vector(dzdy.mat)
  if(is.null(breaks)) cutV=seq(floor(minz),ceiling(wz),length.out=(nbr.intervals+1)) #determine seabed intervals
  else {
    cutV=breaks;ncut=length(cutV)
    if(breaks[1]>floor(minz)) cutV[1]=floor(minz)
    if(breaks[ncut]<ceiling(wz)) cutV[ncut]=ceiling(wz)
  }
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
