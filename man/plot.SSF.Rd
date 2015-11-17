\name{plot.SSF}
\alias{plot.SSF}

\title{ Graphic output of the PAMM function }
\description{
  provide graphic interpretation of the simulation results
}

\usage{
\method{plot}{SSF}(x,\dots)
}

\arguments{
  \item{x}{ an SSF object }
  \item{\dots}{ potentially further arguments to pass to methods }
}

\author{ Julien Martin }

\seealso{\code{\link{SSF}}, \code{\link{plot.PAMM}}, \code{\link{plot.EAMM}} }
\examples{
\dontrun{
   oursSSF=SSF(50,100,10,c(0.4,0.1,0.6,0))
   plot(oursSSF)
   }
}


\keyword{ misc }

