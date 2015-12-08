\name{EAMM}
\alias{EAMM}

\title{Simulation function for exploratory power analysis for random effects}

\description{
 Given a specific sample size, fixed number of group and replicates per group, the function simulate different variance-covariance structure and assess p-values and power of random intercept and random slope
}

\usage{
EAMM(numsim, group, repl, fixed = c(0, 1, 0), VI = seq(0.05, 0.95,0.05),
      VS = seq(0.05, 0.5, 0.05), CoIS = 0, relIS = "cor", n.X, autocorr.X,
          X.dist, intercept = 0, heteroscedasticity = c("null"), mer.sim, mer.model)
}

\arguments{
  \item{numsim}{number of simulation for each step}
  \item{group}{number of group}
  \item{repl}{ number of replicates per group}
  \item{fixed}{vector of lenght 3 with mean, variance and estimate of
    fixed effect to simulate. Default: \code{c(0,1,0)} }
  \item{VI}{variance component of intercept. Could be specified as a
    vector. Default: \code{seq(0.05,0.95,0.05)}}
  \item{VS}{variance component of slope. Could be specified as a vector.
    Default: \code{seq(0.05,0.5,0.05)}}
  \item{CoIS}{value of correlation or covariance between random intercept
    and random slope. Default: 0}
  \item{relIS}{"cor" or "cov" set the type of relation give in CoIS. By
    default the relation is set to correlation}

  \item{n.X}{ number of different values to simulate for the fixed effect (covariate).
   If \code{NA}, all values of X are independent between groups. If the value specified
    is equivalent to the number of replicates per group, \code{repl}, then all groups
     are observed for the same values of the covariate.  Default: \code{NA} }
  \item{autocorr.X}{ correlation between two successive covariate value for a group. Default: \code{0} }
  \item{X.dist}{ specify the distribution of the fixed effect. Only "gaussian" (normal distribution) and
    "unif" (uniform distribution) are accepted actually. Default: \code{"gaussian"} }

  \item{intercept}{a numeric value giving the expected intercept value.
    Default: 0 }
  \item{heteroscedasticity}{a vector specifying heterogeneity in residual
    variance across X.  If \code{c("null")} residual variance is homogeneous
    across X. If \code{c("power",t1,t2)} models heterogeneity with a constant
    plus power variance function. 
    Letting \eqn{v} denote the variance covariate and \eqn{\sigma^2(v)}{s2(v)}
    denote the variance function evaluated at \eqn{v}, the constant plus power
    variance function is defined as \eqn{\sigma^2(v) = (\theta_1 + |v|^{\theta_2})^2}{s2(v) = (t1 + |v|^t2)^2},
    where \eqn{\theta_1,\theta_2}{t1, t2} are the variance function coefficients.
    If \code{c("exp",t)},models heterogeneity with an
    exponential variance function. Letting \eqn{v} denote the variance covariate and \eqn{\sigma^2(v)}{s2(v)}
  denote the variance function evaluated at \eqn{v}, the exponential
  variance function is defined as \eqn{\sigma^2(v) = e^{2 * \theta * v}}{s2(v) = exp(2* t * v)}, where \eqn{\theta}{t} is the variance
  function coefficient. Default:"Null" }
  \item{mer.sim}{Use the simluate.merMod function to simulate the data. Potentially faster for large dataset but more restricted in terms of options}
  \item{mer.model}{Simulate the data based on a existing data and model structure from a lmer object. Should be specified as a list of 3 components: a mer object fitted via lmer, an environmental covariate for which to test the random slope, a random effect (e.g. \code{list(fm1,"Days","Subject"})}
}

\details{
 P-values for random effects are estimated using a log-likelihood ratio
 test between two models with and without the effect. Power represent
 the percentage of simulations providing a significant p-value for a
 given random structure.
 Residual variance, e, is calculted as 1-VI.
}

\value{
  data frame reporting estimated P-values and power with CI for random
  intercept and random slope  
 }

\references{
Martin, Nussey, Wilson and Reale Submitted Measuring between-individual
variation in reaction norms in field and experimental studies: a power
analysis of random regression models. Methods in Ecology and Evolution.  
}
\author{ Julien Martin }

\section{Warning}{the simulation is based on a balanced data set with
  unrelated group}

\seealso{\code{\link{PAMM}}, \code{\link{SSF}}, \code{\link{plot.EAMM}} }

\examples{
\dontrun{
  ours=EAMM(numsim=10,group=10,repl=4,fixed=c(0,1,1),VI=seq(0.05,0.3,0.05),
             VS=seq(0.05,0.2,0.05) )
  plot(ours, "both")
  ours=EAMM(numsim=10, mer.model=list(fm1,"Days","Subject"), VI=seq(0.2,0.5,0.05),
  	 VS=seq(0.05,0.2,0.05) )
   }
}

\keyword{ misc }

