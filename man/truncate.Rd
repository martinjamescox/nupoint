\name{truncate}
\alias{truncate}
\title{Truncates all detections beyond radial distance of argument value}
\usage{
  truncate(trunc.prop = 0.9, sightings)
}
\arguments{
  \item{trunc.prop}{proportion of maximum distance at which
  to truncate (default 0.9)}

  \item{sightings}{(whether simulated or actual)}
}
\value{
  truncation radius

  sightings data frame without detections beyond 0.9*(max
  detection distance) delicately ensure that radial
  distances in sightings are called "d" not "r"
}
\description{
  Truncates all detections beyond radial distance of
  argument value
}

