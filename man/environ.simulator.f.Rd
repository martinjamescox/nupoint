\name{environ.simulator.f}
\alias{environ.simulator.f}
\title{Simulate data with non-parallel gradient}
\usage{
  environ.simulator.f(pars = c(60, 10, 50), z.mat = NULL,
    xlim = c(0, 200), ylim = c(0, 100),
    grid.resolution = 1, grad.type = "NORM",
    det.type = "HNORM", observer.coords = c(100, 0),
    nbr.targets = 350,
    environment.simulator.control = list(c(X = 50, Y = 10, sd = 60), c(X = 90, Y = 0, sd = 30)),
    mask.mat = NULL, mask.ang = 0, plot = TRUE,
    perp.lines = NULL, n = NULL)
}
\arguments{
  \item{pars}{Parameters of the gradient and detection
  function (respectively)}

  \item{z.mat}{Matrix of depth (or other covariate) values}

  \item{xlim}{Study area limits in x-direction}

  \item{ylim}{Study area limits in y-direction}

  \item{grid.resolution}{Separation distance of covariate
  values}

  \item{grad.type}{Distribution used for animal gradient}

  \item{det.type}{Distribution used for detection function
  (HNORM or HAZARD??)}

  \item{observer.coords}{Location of observer}

  \item{nbr.targets}{Number of animals to simulate}

  \item{environment.simulator.control}{List defining
  location of habitat patches being simulated}

  \item{mask.mat}{Habitat mask (not used)}

  \item{mask.ang}{Habitat maks angles (not used)}

  \item{plot}{Flag to indicating plotting}

  \item{perp.lines}{Structure
  list(nbr.transects=2,min.stop=70,max.stop=100)}

  \item{n}{Number of parameters passed in associated with
  gradient and detection function}
}
\value{
  list consisting of \describe{ \item{sightings}{the
  simulated sightings} \item{rd.mat}{matrix of radial
  distances to grid points} \item{z.mat}{covariate value at
  grid points} \item{zGradMat}{gradient of covariate values
  at grid points} \item{settings}{characteristics
  describing the simulation such as grid description,
  detection and distribution models} }
}
\description{
  Simulate data with non-parallel gradient
}
\seealso{
  detectF
}

