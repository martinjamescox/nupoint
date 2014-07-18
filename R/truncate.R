#'Truncates all detections beyond radial distance of argument value
#'
#'@param trunc.prop proportion of maximum distance at which to truncate (default 0.9)
#'@param sightings (whether simulated or actual)
#'
#'@return truncation radius
#'@return sightings data frame without detections beyond 0.9*(max detection distance)
#'delicately ensure that radial distances in sightings are called "d" not "r"
#'
truncate <- function(trunc.prop=0.9, sightings) {
  max.detect <- max(sightings$d)
  trunc.radius <- trunc.prop * max.detect
  new.sightings <- sightings[sightings$d<=trunc.radius, ]
  return(list(trunc.radius=trunc.radius, sightings=new.sightings))
}