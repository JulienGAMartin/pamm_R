\name{plot.PAMM}
\alias{plot.PAMM}

\title{ Graphic output of the PAMM function }
\description{
  provide graphic interpretation of the simulation results
}

\usage{
\method{plot}{PAMM}(x, graphtype = "both", nbgroup, nbrepl, fun3D="wireframe",\dots)
}

\arguments{
  \item{x}{ a PAMM object }
  \item{graphtype}{ "group", "repl" or "both"
    "group" give graphs with varying number of ID and with a fixed number
      of replicates specified in nbrepl
    "repl" give graphs with varying number of replicates and with a fixed
      number of ID specified in nbgroup
    "both" 3-D plot. see also fun3D argument. Note: useful only with multiple
      group size and multiple number of replicates.
      }
  \item{nbgroup}{ number of group for which plots the output. Necessary for "repl" type of graph }
  \item{nbrepl}{ number of replicates for which plots the output. Necessary for "group" type of graph }
  \item{fun3D}{plot function used to plot the 3D graph.
  	"wireframe" uses lattice,
  	"persp" uses  graphics
  	and "open3d" uses rgl}
  \item{\dots}{ potentially further arguments to pass to methods }
}

\details{
	Parameters phi, theta, ticktype ("simple" or "detailed"), nticks, nbcol
	from \code{\link{persp}} function could also be specified for 3D plots.
	In addition, color schemes ("grey", "cm.colors" and "rainbow") and 
	coltype ("restricted" or "0-1") parameters could also be specified for
	 3D plots.
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

