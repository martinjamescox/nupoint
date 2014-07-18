#'Estimate animal abundance from actual data set
#'
#'Code as of 2 August 2013 is incomplete and not in a functional state.
#'
#'This function estimates animal abundance within the study area (grid)
#'by calculating density \eqn{\pi (z(x,y))} as a function of covariate for each grid cell.
#'
#'Calls to appropriate distribution (normal, lognormal, beta, uniform, 
#'mixture of normals) in association with the parameters estimated by
#'the likelihood routine (\code{nupoint.env.fit}) are summed to produce estimate.
#'
#'@param fit.obj fitted object
#'@param truncation distance proportion (default 0.9) such that sightings beyond 0.9*max.r are deleted
#'
#'@return list containing abundance estimate within covered region and
#'abundance estimate for entire study area (assuming grid cells are unit square in area)
#'
#'@details Should your grid cell sizes not be unit square, then multiply the
#'value returned by this function by the grid cell size to produce
#'abundance estimate in the units appropriate for your study.
#'@author Eric Rexstad
#'
#'@references M.J. Cox, D.L. Borchers, D.A. Demer, G.R. Cutter, and A.S. Brierley. 2011. Estimating the density of Antarctic krill (Euphausia superba) from multi-beam echo-sounder observations using distance sampling methods. Journal of the Royal Statistical Society: Series C (Applied Statistics), 60(2):301-316.
#'
#'M.J. Cox, D.L. Borchers and N. Kelly. 2013. nupoint: An R package for density estimation from point transects in the presence of non-uniform animal density Methods in Ecology and Evolution 4(6):589-594
#'
#'Marques, T.A. , Buckland, S.T. , Borchers, D.L. , Tosh, D. and McDonald, R.A.
#'2010.  Point transect sampling along linear features.  Biometrics 66(4):1247-1255.
#'
#'@export

est.abundance.whales <- function(environ.sim.dat, trunc.prop=0.9) {
# nsim <-200
# popn <- numeric(nsim)
# for (k in 1:nsim) {
# environ.sim.dat<-nupoint.env.simulator(pars=c(60,10,50),
#                                        z.mat=NULL,
#                                        xlim=c(0,200),ylim=c(0,100),
#                                        grid.resolution=1,grad.type='NORM',det.type='HNORM',
#                                        observer.coords=c(100,0),nbr.targets=1000,
#                                        environment.simulator.control=
#                                          list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
#                                        mask.mat=NULL,mask.ang=0,plot=FALSE,
#                                        perp.lines=NULL,n=NULL)

test <- truncate(trunc.prop=trunc.prop, sightings=environ.sim.dat$sighting.mat)
trunc.dist <- test$trunc.radius

# replace sightings inside fitted object with truncated sightings
environ.sim.dat[[1]] <- test$sightings
# parameter estimation
browser()
sim.norm.fit<-nupoint.env.fit(pars=environ.sim.dat$settings$pars,
                              z=environ.sim.dat$sighting.mat$z, 
                              rd=environ.sim.dat$sighting.mat$r,  # is it r or d (data or simulation)
                              dzdy=environ.sim.dat$sighting.mat$dzdy,
                              z.mat=environ.sim.dat$z.mat,
                              dzdy.mat=environ.sim.dat$zGradmat,
                              rd.mat=environ.sim.dat$rd.mat,
                              minz=min(environ.sim.dat$z.mat, na.rm=TRUE),
                              wx=environ.sim.dat$wx, #environ.sim.dat$settings$xlim[2],
                              wy=environ.sim.dat$wy, #environ.sim.dat$settings$ylim[2],
                              wz=environ.sim.dat$wz, #max(environ.sim.dat$z.mat),
                              grad.type=environ.sim.dat$settings$grad.type,
                              det.type=environ.sim.dat$settings$det.type,
                              n=NULL,lower.b=rep(1,length(environ.sim.dat$settings$pars))
                              ,upper.b=rep(100,length(environ.sim.dat$settings$pars)))
#  estimate P for HT
#  truncate the grid at the truncation distance 
new.rdmat <- environ.sim.dat$rd.mat
new.zmat <- environ.sim.dat$z.mat
new.zgrad <- environ.sim.dat$zGradmat
for (i in seq(1:dim(environ.sim.dat$rd.mat)[1])) {
  for (j in seq(1:dim(environ.sim.dat$rd.mat)[2])) {
    if (new.rdmat[i,j]>trunc.dist) {
      new.rdmat[i,j] <- NA
      new.zmat[i,j] <- NA
      new.zgrad[i,j] <- NA
  }
  }
}
gradient.model <- environ.sim.dat$settings$grad.type
detection.model <- environ.sim.dat$settings$det.type 
browser()
#  following two lines need fixing for non-norm,hnorm combination
mat.g <- detectF(new.rdmat[!is.na(new.rdmat)], detection.model, sim.norm.fit$par[3])
mat.pi <- pi.z.f(gradient.model, pars=sim.norm.fit$par[1:2], z=new.zmat[!is.na(new.zmat)], 
                 z.lim=c(min(new.zmat, na.rm=TRUE), max(new.zmat, na.rm=TRUE)))
#  Abundance within truncation zone
Nhat.a <- dim(environ.sim.dat$sightings)[1]/sum(mat.g*mat.pi*abs(new.zgrad[!is.na(new.zgrad)])/(1*environ.sim.dat$settings$xlim[2]))
#   Scale Nhat.a to entire study area by dividing by integral pi(x,y) in region a
divisor <- sum(mat.pi*abs(new.zgrad[!is.na(new.zgrad)])/(1*environ.sim.dat$settings$xlim[2]))
print(divisor)
Nhat.region <- Nhat.a / divisor
return(list(Nhat.covered=Nhat.a, Nhat.region=Nhat.region))
# popn[k] <- Nhat.region
}
