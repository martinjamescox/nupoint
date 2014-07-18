\name{mnorm.f}
\alias{mnorm.f}
\title{Density of mixture of normal distributions}
\usage{
  mnorm.f(z, pars, z.lim, n, seabedAtObs = NULL,
    plot = FALSE, verbose = FALSE)
}
\arguments{
  \item{z}{observation (seabed depth in the case of the
  whales; depth of swarm in case of multibeam)}

  \item{pars}{parameter
  vector:\eqn{\mu_1,\sigma_1,\alpha_1,...,\mu_n,\sd_n},
  there is no alpha_n.}

  \item{z.lim}{z-dimension (whale) or y-dimension
  (multi-beam) truncation distance - used to scale the
  density.}

  \item{n}{number of multivariate normal distributions}

  \item{seabedAtObs}{for multibeam (krill) only. this
  allows density to be rescaled on an
  observation-by-observation basis.  must be the same
  length as z.}

  \item{plot}{TRUE/FALSE plot density of mixture normal
  distributions}

  \item{verbose}{print parameter vectors split into mean,
  variance and weights.}
}
\value{
  scaled density
}
\description{
  Density of mixture of normal distributions
}

