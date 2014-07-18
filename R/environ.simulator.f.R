#' Simulate data with non-parallel gradient
#' @param pars Parameters of the gradient and detection function (respectively)
#' @param z.mat Matrix of depth (or other covariate) values
#' @param xlim Study area limits in x-direction
#' @param ylim Study area limits in y-direction
#' @param grid.resolution Separation distance of covariate values
#' @param grad.type Distribution used for animal gradient
#' @param det.type Distribution used for detection function (HNORM or HAZARD??)
#' @param observer.coords Location of observer
#' @param nbr.targets Number of animals to simulate
#' @param environment.simulator.control List defining location of habitat patches being simulated
#' @param mask.mat Habitat mask (not used)
#' @param mask.ang Habitat maks angles (not used)
#' @param plot Flag to indicating plotting
#' @param perp.lines Structure list(nbr.transects=2,min.stop=70,max.stop=100)
#' @param n Number of parameters passed in associated with gradient and detection function
#' @seealso detectF
#' @return list consisting of 
#' \describe{
#'    \item{sightings}{the simulated sightings}
#'    \item{rd.mat}{matrix of radial distances to grid points}
#'    \item{z.mat}{covariate value at grid points}
#'    \item{zGradMat}{gradient of covariate values at grid points}
#'    \item{settings}{characteristics describing the simulation such as grid description, detection and distribution models}
#'    }
environ.simulator.f <- function(pars=c(60,10,50),z.mat=NULL,xlim=c(0,200),ylim=c(0,100),
                       grid.resolution=1,grad.type='NORM',det.type='HNORM',
                       observer.coords=c(100,0),nbr.targets=350,
                       environment.simulator.control=list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
                       mask.mat=NULL,mask.ang=0,plot=TRUE,
                       perp.lines=NULL,n=NULL){
  if(length(mask.mat)!=0 & mask.ang>0)
	{
	
		warning('Please specify either a mask matrix (mask.mat) OR mask angle (mask.ang)')
		return('error: both mask.mat and mask.ang specified')
	}
  message('Commencing data simulation')
  
  #1) unpack variables
  min.x.coord=xlim[1]
  max.x.coord=xlim[2]
  min.y.coord=ylim[1]
  max.y.coord=ylim[2]
  shore.x=observer.coords[1]
  shore.y=observer.coords[2]
    
  #2) construct survey area
  #2-1)coordinate lists
  x.coord=seq((min.x.coord+(grid.resolution/2)),(max.x.coord-(grid.resolution/2)),grid.resolution)
  y.coord=seq((min.y.coord+(grid.resolution/2)),(max.y.coord-(grid.resolution/2)),grid.resolution)
  

  #some error checking
  if(length(z.mat)!=0)
  {
	if(nrow(z.mat)!=length(y.coord))
		{
			warning('Mismatch between z.mat dimensions and ylim and/or grid.resolution')
			return('z.mat or y-coordinate vector size mismatch')
		}
		if(ncol(z.mat)!=length(x.coord))
		{
			warning('Mismatch between z.mat dimensions and xlim and/or grid.resolution')
			return('z.mat or x-coordinate vector size mismatch')
		}
  }

  #3) distance matrix
  rd.f=function(x,y) sqrt((x-shore.y)**2+(y-shore.x)**2)
  dist.mat=outer(y.coord,x.coord,rd.f)
  #3-1) detection function
  parList=par.unpack.F(grad.type,pars,n)
  grad.pars=parList[[1]]; det.par=parList[[2]]
  
  det.mat=detectF(rd=dist.mat,det.type=det.type,det.par=det.par)
  #4) Seabed depth matrix.  The states of the of the z.mat are:
  #z.mat = TRUE - create a seabed matrix 
  #z.mat = a matrix (2D array of seabed values) these are used in the function. see also
  #dzdy.mat
  #z.mat = FALSE - no seabed is used and the simulator will create a multibeam type dataset.
  if(length(z.mat)==0){
    #seabed.simulator.control form = list[[one list location for each knot]]
    #c(difference knot X coordinate, difference knot Y coordinate, knot sigma parameter)
    #example: seabed.simulator.control=list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30))
    #4-1) depth matrix
    z.mat=matrix(rep(y.coord,length(x.coord)),nrow=length(y.coord))
    #4-2) make irregular seabed
    coords=expand.grid(y.coord,x.coord)
    coords$depth=(as.vector(z.mat)) #jitter
    names(coords)[1:2]=c("Y","X")
    
    depthAdjustMat=matrix(0,nrow=nrow(coords),ncol=length(environment.simulator.control))
    for(i in 1:length(environment.simulator.control)){ #start simulator control loop.
      dist = sqrt((coords$X-environment.simulator.control[[i]][1])**2+ (coords$Y-environment.simulator.control[[i]][2])**2)
      prob = exp(-dist^2/(2*environment.simulator.control[[i]][3]^2))
      depthAdjustMat[,i] = coords$depth * (1+prob)} #end of simulator control loop.
    coords$z=apply(depthAdjustMat,1,mean)
    
    z.mat=matrix(coords$z,nrow=length(y.coord))
    z.max=max(z.mat,na.rm=TRUE)
  } #end of if(z.mat==NULL)
  
  #bit more error checking
  if(length(mask.mat)!=0){
  if(nrow(mask.mat)!=nrow(z.mat) | ncol(mask.mat)!=ncol(z.mat))
	warning('Array size mismatch in mask.mat and z.mat')
	return('Array size mismatch in mask.mat and z.mat')
  }
  
  y.mat=matrix(rep(y.coord,length(x.coord)),nrow=nrow(z.mat))
 
  #4-3) change in seabed depth matrix dz/dy
  #spline.call=function(y){
  #  f= splinefun(y.coord,y,method = "natural")
  #}
  #sp.list=apply(z.mat,2,spline.call)
  #
  #dzdy.mat=matrix(NA,nrow=nrow(z.mat),ncol=ncol(z.mat))
  #for(i in 1:ncol(z.mat))
   # dzdy.mat[,i]=abs(sp.list[[i]](1:nrow(z.mat),deriv=1))
  
  #rate of change in seabed depth matrix
	spline.call=function(y,z) f= splinefun(y,z)

	sp.list=list()

	for(i in 1:ncol(z.mat))
		sp.list[[i]]=spline.call(y.mat[,i],z.mat[,i])

	zGradmat=matrix(NA,nrow=nrow(z.mat),ncol=ncol(z.mat))
	for(i in 1:ncol(zGradmat))
		zGradmat[,i]=abs(sp.list[[i]](y.mat[,i],deriv=1))

  
  
  #5) mask matrix - not coded into simulator - probably replace depth matrix with
  if(mask.ang>0){
    rel.x.mat=Mod(matrix(rep(x.coord,length(y.coord)),nrow=nrow(z.mat),byrow=TRUE)-shore.x)
    #y.mat=matrix(rep(y.coord,length(x.coord)),nrow=nrow(z.mat))
    mask.dist=rel.x.mat[1,]*tan(pi*(mask.ang/180))
    
    for(i in 1:ncol(y.mat))
      z.mat[which(y.mat[,i]<mask.dist[i]),i]=NA
  } #end of mask creation
  
  #6) seabed depth preference $\pi_y$ matrix
  pi.y=pi.z.f(g.type=grad.type,pars=grad.pars,z=z.mat,
              z.lim=c(min(z.mat,na.rm=TRUE),max(z.mat,na.rm=TRUE)),n=n)

			  #7) simulate detections:
  x.coord.vec=y.coord.vec=z.coord.vec=dzdy.vec=
    cue.dist.vec=detect.vec=vector(mode="numeric",length=nbr.targets)
  
  for(i in 1:nbr.targets)
  {
    x.loc=sample(1:ncol(pi.y),size=1)
    x.coord.vec[i]=x.coord[x.loc]
    y.locations=1:nrow(pi.y)
    rm.ele=which(is.na(z.mat[,x.loc]))
    if(length(rm.ele)>0){
      y.loc=sample(y.locations[-rm.ele],size=1,prob=pi.y[-rm.ele,x.loc])}else{
        y.loc=sample(y.locations,size=1,prob=pi.y[,x.loc])}
    y.coord.vec[i]=y.coord[y.loc]
    z.coord.vec[i]=z.mat[y.loc,x.loc]
    cue.dist.vec[i]=dist.mat[y.loc,x.loc]
    dzdy.vec[i]=abs(sp.list[[x.loc]](y.coord.vec[i],deriv=1))
    detect.vec[i]=rbinom(n=1, size=1, prob=detectF(rd=cue.dist.vec[i],det.type=det.type,det.par=det.par))# exp(-cue.dist.vec[i]^2/(2*s2^2)))
    }
    
  sim=data.frame(x=x.coord.vec,y=y.coord.vec,z=z.coord.vec,dzdy=dzdy.vec,
                 d=cue.dist.vec,detect=detect.vec)
  
  #8) perpendicular lines
  if(length(perp.lines)>0){ #simulate perpendicular lines
    #data structure: perp.lines=list(nbr.transects=2,min.stop=70,max.stop=100)
   
    #error checking
    if(perp.lines$min.stop>max(y.coord) | perp.lines$max.stop>max(y.coord))
      warning('simulatorF call error. Simulated transect stop y-coordinate misalignment: 
              transect max coordinate potentially exceeds simulated area maximum y-coordinate.')
    if(length(perp.lines)!=3)
      warning('simulatorF call error. Simulator data structure input incorrect length.')
    #by transect
    transect.x.coord.start = transect.x.coord.end = transect.y.coord.start=
      transect.y.coord.end = vector(mode='numeric',length = perp.lines$nbr.transects)
    
    transect.x.coord.start=transect.x.coord.end=
      x.coord[round(seq(1,length(x.coord),length.out=perp.lines$nbr.transects+2),0)[-c(1,(perp.lines$nbr.transects+2))]]
    transect.y.coord.start=rep(0,perp.lines$nbr.transects)
    transect.y.coord.end=round(runif(perp.lines$nbr.transects,
                                     perp.lines$min.stop,perp.lines$max.stop),0)
    #transect.y.coord.start[i]=mask.dist[which.min(Mod(transect.x.coord.start[i]-x.coord))]
    if(mask.ang>0)
      transect.y.coord.start=Mod(transect.x.coord.start-shore.x)*tan(pi*(mask.ang/180))
    transect.x.coord.start<<-transect.x.coord.start
    perp.out.dat=sim[which(sim$x %in% transect.x.coord.start),]
	#perp.out.dat=perp.out.dat[,-which(names(perp.out.dat)=='detect')] #remove detect column from the perpendicular transect data.
    perp.out.dat$transect=as.numeric(as.factor(perp.out.dat$x))#sort(rep(tran.nbr[tran.table>0],tran.table[tran.table>0]))
  
    perp.out.info=data.frame(startx=transect.x.coord.start,
                                         stopx=transect.x.coord.end,
                                         starty=transect.y.coord.start,
                                         stopy=transect.y.coord.end)
    perp.out.y=y.mat[,which(x.coord %in% transect.x.coord.start)]
    perp.out.z=z.mat[,which(x.coord %in% transect.x.coord.start)]
    perp.out.dzdy=zGradmat[,which(x.coord %in% transect.x.coord.start)]
    } #end perpendicular lines.
  
  #9) plots
  if(plot)
  {
    polyF=function(mask.ang)
    {
      if(mask.ang>0)  
      {  
          polygon(c(shore.y,0,max(mask.dist),shore.y),c(shore.x,max(x.coord)+grid.resolution,
                    max(x.coord)+grid.resolution,shore.x),col='grey')
          polygon(c(shore.y,0,max(mask.dist),shore.y),c(shore.x,min(x.coord)-grid.resolution,
                    min(x.coord)-grid.resolution,shore.x),col='grey')
      }
    } #end polyF
    
    require(fields,quietly = TRUE)
    par(mfrow=c(2,1),mar=c(3,2,2,2))
    #image.plot(y=x.coord,x=y.coord,z=dist.mat,main="Radial distance")
    #polyF(mask.ang)
    #points(shore.y,shore.x,cex=3,pch=19)
    
    image.plot(x=x.coord,y=y.coord,z=t(det.mat),xlab='X',ylab='Y',main="Detection function, g(r)")
    polyF(mask.ang)
    if(length(perp.lines)>0){
      for(i in 1:nrow(perp.out.info)) lines(c(perp.out.info$startx[i],perp.out.info$stopx[i]),
											c(perp.out.info$starty[i],perp.out.info$stopy[i]))}
    
	points(shore.x,shore.y,cex=3,pch=19)
    points(x.coord.vec,y.coord.vec,pch=19,cex=0.3)
    points(sim$x[sim$detect==1],sim$y[sim$detect==1],col="white",pch=19,cex=0.5)
    
    image.plot(x=x.coord,y=y.coord,z=t(z.mat),xlab='X',ylab='Y',main='z(x,y)')#expression(paste("Environmental variable, ",pi[xy],'(z(x,y))')))
    polyF(mask.ang)
    points(shore.x,shore.y,cex=3,pch=19)
    points(x.coord.vec,y.coord.vec,pch=19,cex=0.3)
    
    points(sim$x[sim$detect==1],sim$y[sim$detect==1],col="white",pch=19,cex=0.5)
     if(length(perp.lines)>0){
      for(i in 1:nrow(perp.out.info)) lines(c(perp.out.info$startx[i],perp.out.info$stopx[i]),
											c(perp.out.info$starty[i],perp.out.info$stopy[i]))}
    
    #image.plot(y=x.coord,x=y.coord,z=zGradmat,main="Environmental gradient")
    #polyF(mask.ang)
    #points(shore.y,shore.x,cex=3,pch=19)
  }
  sim=sim[sim$detect==1,] #retain only sighted cues.
  settings=list(pars=pars,xlim=xlim,ylim=ylim,
                grid.resolution=grid.resolution,grad.type=grad.type,det.type=det.type,
                observer.coords=observer.coords,nbr.targets=nbr.targets)
  out=list(sightings=sim,rd.mat=dist.mat,z.mat=z.mat,zGradmat=zGradmat,settings=settings)
  if(length(perp.lines)>0) {out[[length(out)+1]]=perp.out.y
                            out[[length(out)+1]]=perp.out.z
                            out[[length(out)+1]]=perp.out.dzdy
                            out[[length(out)+1]]=perp.out.info
                            out[[length(out)+1]]=perp.out.dat
                            names(out)[(length(out)-4):length(out)]=c('perpendicular.y',
                                                                      'perpendicular.z',
                                                                      'perpendicular.dzdy',
                                                      'perpendicular.transect.info',
                                                      'perpendicular.transect.data')}
  return(out)

}
