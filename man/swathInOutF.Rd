\name{swathInOutF}
\alias{swathInOutF}
\title{Determine if a point is inside or outside of the MBE swath (point transect)}
\usage{
  swathInOutF(x, y, w, theta)
}
\arguments{
  \item{x}{and y coordinates}

  \item{w}{truncation distance}

  \item{theta}{swath (point transect) width (radians)}
}
\value{
  binary 0 = outside swath; 1= inside
}
\description{
  Determine if a point is inside or outside of the MBE
  swath (point transect)
}

