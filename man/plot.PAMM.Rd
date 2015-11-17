\name{plot.PAMM}
\alias{plot.PAMM}

\title{ Graphic output of the PAMM function }
\description{
  provide graphic interpretation of the simulation results
}

\usage{
\method{plot}{PAMM}(x, graphtype = "both", nbgroup, nbrepl,\dots)
}

\arguments{
  \item{x}{ a PAMM object }
  \item{graphtype}{ "group", "repl", "both" or "both.dyn"
    "group" give graphs with varying number of ID and with a fixed number
      of replicates specified in nbrepl
    "repl" give graphs with varying number of replicates and with a fixed
      number of ID specified in nbgroup
    "both" 3-D plot using lattice package. Note: useful only with multiple
      group size and multiple number of replicates.
    "both.dyn" 3D plot using rgl package. Note: useful only with multiple group
      size and multiple number of replicates}
  \item{nbgroup}{ number of group for which plots the output. Necessary for "repl" type of graph }
  \item{nbrepl}{ number of replicates for which plots the output. Necessary for "group" type of graph }
  \item{\dots}{ potentially further arguments to pass to methods }
}

\author{ Julien Martin }

\seealso{\code{\link{PAMM}}, \code{\link{plot.EAMM}}, \code{\link{plot.SSF}} }
\examples{
\dontrun{
  ours=PAMM(numsim=10,group=c(10,20,30),repl=c(2,4,6),
     randompart=c(0.4,0.1,0.5,0.1),fixed=c(0,1,0.7))
  plot(ours, "both")
  plot(ours, "group",nbrepl=4)
  plot(ours,"repl",nbgroup=20)
   }
}

\keyword{ misc }

