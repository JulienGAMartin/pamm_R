\name{plot.EAMM}
\alias{plot.EAMM}

\title{ Graphic output of the EAMM function }
\description{
  provide graphic interpretation of the simulation results
}

\usage{
\method{plot}{EAMM}(x, graphtype = "both", vi, vs,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an EAMM object }
  \item{graphtype}{ "VI", "VS","both" or "both.dyn"
    "VI" give graphs with varying variance component of intercept and
    with a fixed variance component for slope specified in vs argument
    "VS" give graphs with varying variance component for slope and with a
    fixed variance component of intercept specified in vi argument
    "both" 3-D plot using lattice package
    "both.dyn" 3D plot using rgl package}
  \item{vi}{ VI for which plots the output. Necessary for "VS" type of graph }
  \item{vs}{ VS for which plots the output. Necessary for "VI" type of graph }
  \item{\dots}{ potentially further arguments to pass to methods }
}

\author{ Julien Martin }

\seealso{\code{\link{EAMM}}, \code{\link{plot.PAMM}}, \code{\link{plot.SSF}} }
\examples{
\dontrun{
  ours=EAMM(numsim=10,group=10,repl=4,VI=seq(0.05,0.95,0.05),VS=c(0.05,0.1))
  plot(ours, "both")
  plot(ours, "VI",vi=0.4)
  plot(ours,"VS",vs=0.1) 
   } 
}

\keyword{ misc }

