#'Determine if a point is inside or outside of the MBE swath (point transect)
#'@param x and y coordinates
#'@param w truncation distance
#'@param theta swath (point transect) width (radians)
#'@return binary 0 = outside swath; 1= inside
swathInOutF <- function(x,y,w,theta)
{
  ymax=w*cos(theta)
  loc=rep(1,length(x))
  loc[x>(y*tan(theta))]=NA
  loc[y>ymax & x > sqrt(w**2-y**2)]=NA
  return(loc)
}