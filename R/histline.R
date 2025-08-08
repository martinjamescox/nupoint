#'@title makes histograms in line form
#'@description Takes bar heights (height) and cutbpoints (breaks), and constructs a line-only 
#'@param height vector of bar heights, and breaks vector of cutpoints
#'@param breaks vector of cutpoints
#'@param lineonly if TRUE, uses lines() to plot histogram, else uses plot()
#'@param outline if TRUE, only outline of histogram is plotted
#'@param fill if TRUE, uses polygon() to fill bars
#'@param ylim y-axis limits
#'@param xlab x-axis label
#'@param ylab y-axis label
#'@param ... additional arguments passed to plot() or lines() if lineonly=TRUE, or polygon() if fill=TRUE
histline <- function(height,
                     breaks,
                     lineonly=FALSE,
                     outline=FALSE,
                     fill=FALSE,
                     ylim=range(height),
                     xlab="x",
                     ylab="y",
                     ...)
#' histogram from them using the function plot() (if lineonly==FALSE) or lines()
#' (if lineonly==TRUE). 
#' If fill==TRUE, uses polygon() to fill bars
#' If fill==TRUE, valid arguments to plot() or lines() are passed via argument(s) "..."
#' If outline==TRUE, only outline of histogram is plotted
#' If fill!=TRUE, valid arguments to polygon() are passed via argument(s) "..."

{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}