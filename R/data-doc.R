#' Example data from a parallel density gradient: a multibeam survey of krill
#'
#' Krill data set illustrating a parallel density gradient example. Observations of krill swarms made using a multibeam echosounder.
#' 
#' @format A data frame with 251 rows and 7 variables:
#' \describe{
#'   \item{transect}{Transect number}
#'   \item{x}{Shoal cross-track distance, m}
#'   \item{y}{Shoal depth, m}
#'   \item{r}{Shoal radial distance, m}
#'   \item{theta}{Shoal angle, radians}
#'   \item{bio.g}{Shoal biomass, g}
#'   \item{z}{Seabed depth under shoal, m.}
#' }
#'
#' @source For data collection details see: Cox et al (2009) Multibeam echosounder observations reveal interactions between Antarctic krill and air-breathing predators. *Marine Ecology Progress Series*, 378, 199–209. Acknowledgements: Multi-beam instrument was loaned by J. Condiotty of Simrad USA. Data were collected in association with an NSF-funded (grant #06-OPP-33939) investigation of the Livingston Island nearshore environment. Support provided by US Antarctic Marine Living Resources Program, and engineering by Sea Technology Services.
#'
#' @docType data
#' @keywords datasets
#' @name nupoint_krill
NULL

#' Example beaked whale survey data
#'
#' A synthetic or subsetted dataset illustrating detection functions and environmental gradients for beaked whale surveys. Includes individual detections and environmental matrices used for spatial modeling.
#'
#' @format A list with 14 elements:
#' \describe{
#'   \item{sighting.mat}{A data frame of 191 detections with columns:
#'     \describe{
#'       \item{x}{Sighting x-coordinate (e.g., easting or distance, m)}
#'       \item{y}{Sighting y-coordinate (e.g., northing or distance, m)}
#'       \item{r}{Radial distance from observer (m)}
#'       \item{z}{Depth at sighting (m)}
#'       \item{dzdy}{Vertical gradient of depth at y-location}
#'       \item{obs.period}{Observation period ID (e.g., "P01")}
#'     }
#'   }
#'   \item{x.mat}{Matrix of x coordinates used for interpolation or gridding (132 × 90)}
#'   \item{y.mat}{Matrix of y coordinates used for interpolation or gridding (132 × 90)}
#'   \item{rd.mat}{Matrix of radial distances from observer to grid cells (132 × 90)}
#'   \item{z.mat}{Matrix of depth values (132 × 90), with missing values as NA}
#'   \item{zGradmat}{Matrix of vertical gradients of depth (132 × 90)}
#'   \item{x}{Vector of unique x locations (length 132)}
#'   \item{y}{Vector of unique y locations (length 90)}
#'   \item{obsx}{Observer x-coordinate}
#'   \item{obsy}{Observer y-coordinate}
#'   \item{wx}{Maximum x coordinate of detection window}
#'   \item{wy}{Maximum y coordinate of detection window}
#'   \item{wz}{Maximum depth of detection window}
#'   \item{minz}{Minimum observed depth}
#' }
#'
#' @source Arranz, P., Borchers, D. L., de Soto, N. A., Johnson, M. P., & Cox, M. J. (2014). A new method to study inshore whale cue distribution from land-based observations. Marine Mammal Science, 30(2), 810-818. — subset used for example modeling purposes.
#'
#' @docType data
#' @keywords datasets
#' @name sightings
NULL
