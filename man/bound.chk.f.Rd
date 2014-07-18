\name{bound.chk.f}
\alias{bound.chk.f}
\title{Check starting values are within parameter space bounds when using ,method="L-BFGS-B" in optim()}
\usage{
  bound.chk.f(pars, lower, upper, f.n)
}
\arguments{
  \item{pars}{initial parameter estimates}

  \item{lower}{lower bound}

  \item{upper}{upper bound}

  \item{f.n}{character string of function within which the
  bounds are being checked e.g. 'mbe.fit.f for multibeam.}
}
\value{
  warning flag logical FALSE = estimates are fine. TRUE =
  estimates outside upper or lower bounds.
}
\description{
  Check starting values are within parameter space bounds
  when using ,method="L-BFGS-B" in optim()
}

