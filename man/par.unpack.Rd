\name{par.unpack}
\alias{par.unpack}
\title{Unpack parameters so that optim can run using a single parameter vector.}
\usage{
  par.unpack(grad.type, pars, n)
}
\value{
  list([1]=density gradient parameters; [2]=detection
  function parameters).
}
\description{
  Unpack parameters so that optim can run using a single
  parameter vector.
}

