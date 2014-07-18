#'Simulator for shore-based distance observations of animals distributed along
#'an environmental gradient
#'
#'Creates a simulated survey region and a set of observations collected by a
#'shore-based observer when animals are distributed along an environmental
#'gradient.
#'
#'The user can either pass in a matrix specifying the environment variable
#'(e.g. seabed) using the \code{z.mat} argument or can create an environment
#'variable using the \code{environment.simulator.control} argument. If a matrix
#'is passed into the function via the \code{z.mat} argument, then its
#'dimensions must equal those of the matrix (if any) used in \code{mask.mat}.
#'Further the range of \code{xlim} divided by \code{grid.resolution} must equal
#'the number of rows in the matrix (\code{z.mat} and \code{mask.mat})
#'arguments, and the range of \code{ylim} divided by \code{grid.resolution}
#'must equal the number of columns in \code{z.mat} and \code{mask.mat}.
#'
#'If \code{z.mat} is unspecified, the user can create a \code{z.mat} within the
#'\code{nupoint.env.simulator} function by specifying radial basis functions
#'using the \code{environment.simulator.control} argument, which is a list
#'object comprised of vectors. Each vector describes a radial basis function
#'\code{c(x-coordinate, y-coordinate, standard deviation)}.
#'
#'If \code{slope.control} is unspecified, the default is for the simulated
#'seabed to equal the y-coordinate at a given position. This means the seabed
#'increases linearly with distance from shore (the x-axis). The user can change
#'the default seabed by specifying an intercept parameter
#'\code{slope.control[1]} and slope parameter \code{slope.control[2]} so that a
#'given seabed depth, \code{z} is \code{z}=slope.control[1] + slope.control[2]
#'x y-coordinate.
#'
#'A mask matrix may be specified by the user using the \code{mask.mat} argument
#'to exclude unsuitable environment from the survey region (e.g. land in the
#'case of beaked whales). This matrix must have identical dimensions to
#'\code{z.mat}. Alternatively, areas within the survey region may be excluded
#'using the \code{mask.ang} argument, which removes areas that are less than
#'\code{mask.ang} degrees either side of the observer. \code{mask.ang} is
#'measured from the x-axis.
#'
#'Perpendicular line transects can be simulated using the \code{perp.lines}
#'argument. The default is \code{prep.lines=NULL}, meaning no perpendicular
#'survey lines are simulated, otherwise an even number of integer lines are
#'created and populated with simulated observations.
#'
#'\code{nupoint.env.simulator} returns a list comprised of the following:
#'\code{[[1]]} a data frame of simulated sightings \code{\link{sightings}};
#'\code{[[2]]} \code{rd.mat} a matrix of distances from observer to each point
#'in the simulated survey region; \code{[[3]]} \code{z.mat} a regular matrix of
#'the environmental feature throughout the survey region; \code{[[4]]}
#'\code{zGradmat} regular matrix of the environmental feature gradient in the
#'y-dimension; \code{[[5]]} a sub-list of simulator settings (\code{pars},
#'\code{xlim}, \code{ylim}, \code{grid.resolution}, \code{grad.type},
#'\code{det.type}, \code{observer.coords}, \code{nbr.targets} arguments), and
#'\code{[[6:10]]}: simulated perpendicular transect information if
#'\code{perp.lines!=NULL}.
#'
#'The data frame of simulated sightings is made up of sighting x-coordinate,
#'\code{x}, sighting y-coordinate, \code{y}, sighting environmental feature
#'value, \code{z}, environmental feature gradient, \code{dzdy}, and detection
#'distance \code{d}.
#'
#'Perpendicular transect information is simulated the returned list contains
#'additional information: \code{[[6]]}: \code{perpendicular.y} matrix of
#'y-coordinates, when each column contains the y-values along a given simulated
#'perpendicular transect; \code{[[7]]}: \code{perpendicular.y} matrix of
#'environmental feature values, when each column contains values along a given
#'simulated perpendicular transect; \code{[[8]]}: \code{perpendicular.dzdy}
#'matrix of environmental feature gradient in the y-dimension; \code{[[9]]}:
#'\code{perpendicular.transect.info} perpendicular transect information giving
#'the start and end x-coordinates \code{startx} and \code{stopx}, and the start
#'and end y-coordinates \code{starty} and \code{stopy} for each simulated
#'transect, and \code{[[10]]} \code{perpendicular.transect.data}: observations
#'collected on each simulated transect (identical format to element
#'\code{[[1]]} \code{\link{sightings}} plus a column named \code{transect} that
#'assigns an observation to a perpendicular transect.
#'
#'@param pars \code{vector of parameters to simulate from (see
#'nupoint.env.fit).}
#'@param z.mat \code{matrix:either a regular grid of the environmental feature
#'throughout the survey region, or NULL if bathymetry is simulated using
#'settings in the environment.simulator.control argument (see details).}
#'@param xlim \code{x-coordinate range.}
#'@param ylim \code{y-coordinate range.}
#'@param grid.resolution \code{resolution of simulated survey region.}
#'@param grad.type \code{parametric form of the environmental gradient.  See
#'nupoint.env.fit.}
#'@param det.type \code{parametric form of the detection function
#'c("HNORM","HAZARD"). See detectF.}
#'@param observer.coordinates \code{x,y observer position.}
#'@param nbr.targets \code{Number of targets in survey region. When
#'Pr(detect)=1 this is the number of detected targets.}
#'@param environment.simulator.control \code{list object comprised of vectors.
#'Each vector describes a radial basis function which are used to simulate the
#'environment when z.mat=NULL (see details).}
#'@param slope.control \code{default NULL. Otherwise a two element vector
#'comprising of slope and intercept specifying a linearly increasing seabed
#'function.  See details.}
#'@param mask.mat \code{default NULL.  Otherwise a user defined matrix of
#'dimensions equal to z.mat specifying areas within the survey region to be
#'excluded from analysis (e.g. simulated coastline, see details).}
#'@param mask.ang \code{angle (rad) at observer, measured from baseline (y=0)
#'contains unsuitable habit (e.g. land; see details).  Either specify mask.ang
#'or mask.ang}
#'@param plot \code{logical. Default TRUE.  Returns a panel plot of detection
#'function and envrionmental gradient with all targets and seen targets shown
#'on both plots.}
#'@param perp.lines \code{Default = NULL (see details).  Or a list object
#'describing perpendicular transects lines run from baseline,
#'list(nbr.transects,min.stop,max.stop).  The min.stop and max.stop settings
#'are used to specify the minimum and maximum ltransect length in the
#'y-dimension.}
#'@param n \code{Default NULL. Use only when grad.type=MNORM to specify an
#'integer number of distributions in the normal mixture option for the
#'parametric form of the environmental gradient.  See nupoint.env.fit.}
#'@return list object of simulator settings, survey region descriptors and
#'observations (see details).
#'@seealso \code{\link{detectF}}, \code{\link{nupoint.env.fit}}
#'@export
#'@keywords misc
#'@examples
#'
#'environ.sim.dat=nupoint.env.simulator(pars=c(60,10,50),z.mat=NULL,xlim=c(0,200),ylim=c(0,100),
#'                       grid.resolution=1,grad.type='NORM',det.type='HNORM',
#'                       observer.coords=c(100,0),nbr.targets=350,
#'                       environment.simulator.control=list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
#'                       mask.mat=NULL,mask.ang=0,plot=TRUE,
#'                       perp.lines=NULL,n=NULL)
#'
nupoint.env.simulator <- function(pars=c(60,10,50),z.mat=NULL,xlim=c(0,200),ylim=c(0,100),
                       grid.resolution=1,grad.type='NORM',det.type='HNORM',
                       observer.coords=c(100,0),nbr.targets=350,
                       environment.simulator.control=list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
					   slope.control=NULL,
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
	linear.seabed=y.coord #20130218 modification based on David's email.
	if(length(slope.control)>0)
	{message(paste('Using linear seabed parameters; intercept=',slope.control[1],'slope=',slope.control[2]))
	if(!is.numeric(slope.control)) {warning('Non-numeric slope.control parameters')
		return('nupoint.env.simulator_non_numeric_seabed_slope')} #end error trap
		linear.seabed=slope.control[1]+linear.seabed*slope.control[2]
   browser()
	} else {message('Using y-coordinates to generate a linearly increasing seabed') }
    z.mat=matrix(rep(linear.seabed,length(x.coord)),nrow=length(y.coord)) 
    #20130219: check for environemtnal feature ==0
	if(all(z.mat==0)) {warning('All elements of the environmental feature == 0.  Unable to calculate expected sighting locations.')
		warning('Check z.mat ARG or in the case of environmental preference simulation, check simulation settings: environment.simulator.control, slope.control ARGS')
		return('environmental_feature_values_all_zero')}
	#4-2) make irregular seabed
    coords=expand.grid(y.coord,x.coord)
    coords$depth=(as.vector(z.mat)) 
    names(coords)[1:2]=c("Y","X")
    
    depthAdjustMat=matrix(0,nrow=nrow(coords),ncol=length(environment.simulator.control))
    for(i in 1:length(environment.simulator.control)){ #start simulator control loop.
      dist = sqrt((coords$X-environment.simulator.control[[i]][1])**2+ (coords$Y-environment.simulator.control[[i]][2])**2)
      prob = exp(-dist^2/(2*environment.simulator.control[[i]][3]^2))
      depthAdjustMat[,i] = coords$depth * (1+prob)} #end of simulator control loop.
    coords$z=apply(depthAdjustMat,1,mean)
    coords$z[coords$z<0]=0
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
	#20130219: error check for infinite expected environmental preference. Typically occurs when z.mat == 0
	if(all(is.infinite(pi.y))) {warning('All expected location probabilities infinite')
	warning('Check environmental preference function (grad.type and pars ARGS') 
	warning('Check z.mat ARG or environmental preference simulation settings,environment.simulator.control, slope.control ARGS')
		return('expected_probabilities_all_infinite')}
			  #7) simulate detections:
  x.coord.vec=y.coord.vec=z.coord.vec=dzdy.vec=
    cue.dist.vec=detect.vec=vector(mode="numeric",length=nbr.targets)
  
  #20130219: rewritten 
  xV=rep(x.coord,each=length(y.coord))
  yV=rep(y.coord,each=length(x.coord))
  xM=matrix(xV,ncol=ncol(z.mat),byrow=FALSE)
  yM=matrix(yV,ncol=ncol(z.mat),byrow=TRUE)

  plocV=as.vector(pi.y) #weight vector for sightings locations.
  rows.zmat <- nrow(z.mat)  ###ER 20130724
  for(i in 1:nbr.targets)
  {
    sampleElement=sample(x=1:length(plocV) ,size=1,prob=plocV) #sample a location
	colSamp= 1 + sampleElement %/% rows.zmat
	rowSamp= sampleElement - (rows.zmat * sampleElement %/% rows.zmat)
#     rowSamp <- 1 + sampleElement %/% rows.zmat ###ER 20130724  integer divide
#     colCandidate <- sampleElement %% rows.zmat ###ER 20130724  modulus
#     if (colCandidate == 0) {
#       colSamp <- cols.zmat
#       rowSamp <- rowSamp - 1
#     } else {
#       colSamp <- colCandidate
#     }
# 	if(rowSamp==0 || colSamp==0) browser()
#   if(rowSamp>nrow(xM) || colSamp>ncol(xM)) browser()
    x.coord.vec[i]=xM[rowSamp,colSamp]
	y.coord.vec[i]=yM[rowSamp,colSamp]
	z.coord.vec[i]=z.mat[rowSamp,colSamp]
	cue.dist.vec[i]=dist.mat[rowSamp,colSamp]
	dzdy.vec[i]=abs(sp.list[[rowSamp]](y.coord.vec[i],deriv=1))
	detect.vec[i]=rbinom(n=1, size=1, prob=detectF(rd=cue.dist.vec[i],det.type=det.type,det.par=det.par))
	plocV[sampleElement]=0 #remove sampled location
	
	#x.loc=sample(1:ncol(pi.y),size=1)
    #x.coord.vec[i]=x.coord[x.loc]
    #y.locations=1:nrow(pi.y)
    #rm.ele=which(is.na(z.mat[,x.loc]))
    #if(length(rm.ele)>0){
    #  y.loc=sample(y.locations[-rm.ele],size=1,prob=pi.y[-rm.ele,x.loc])}else{
    #    y.loc=sample(y.locations,size=1,prob=pi.y[,x.loc])}
    #y.coord.vec[i]=y.coord[y.loc]
    #z.coord.vec[i]=z.mat[y.loc,x.loc]
    #cue.dist.vec[i]=dist.mat[y.loc,x.loc]
    #dzdy.vec[i]=abs(sp.list[[x.loc]](y.coord.vec[i],deriv=1))
    #detect.vec[i]=rbinom(n=1, size=1, prob=detectF(rd=cue.dist.vec[i],det.type=det.type,det.par=det.par))# exp(-cue.dist.vec[i]^2/(2*s2^2)))
  }
   #20130219: end of rewrite 
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
