\name{PAMM}
\alias{PAMM}

\title{ Simulation function to assess power of mixed models }
\description{ Given a specific varaince-covariance structure for random
  effect, the function simulate different group size and assess p-values and power of
  random intercept and random slope }

\usage{
PAMM(numsim, group, repl, randompart, fixed, n.X, autocorr.X,
          X.dist, intercept, heteroscedasticity = c("null"), ftype="lmer", mer.sim=FALSE)
}

\arguments{
  \item{numsim}{ number of simulation for each step }
  \item{group}{ number of group. Could be specified as a vector }
  \item{repl}{ number of replicates per group . Could be specified as a vector }
  \item{randompart}{ vector of lenght 4 or 5, with 1: variance component
    of intercept, VI; 2: variance component of slope, VS; 3: residual
    variance, VR; 4: relation between random intercept and random
    slope; 5: "cor" or "cov" determine if the relation 4 between I ans S is a correlation or a covariance. Default: \code{"cor"} }
  \item{fixed}{ vector with mean, variance and estimate of fixed effect to simulate. Default: \code{c(0,1,0)} }

  \item{n.X}{ number of different values to simulate for the fixed effect (covariate).
   If \code{NA}, all values of X are independent between groups. If the value specified
    is equivalent to the number of replicates per group, \code{repl}, then all groups
     are observed for the same values of the covariate.  Default: \code{NA} }
  \item{autocorr.X}{ correlation between two successive covariate value for a group. Default: \code{0} }
  \item{X.dist}{ specify the distribution of the fixed effect. Only "gaussian" (normal distribution) and
    "unif" (uniform distribution) are accepted actually. Default: \code{"gaussian"} }

  \item{intercept}{a numeric value giving the expected intercept value. Default:0 }
  \item{heteroscedasticity}{a vector specifying heterogeneity in residual variance
    across X. If \code{c("null")} residual variance is homogeneous across X. If 
    \code{c("power",t1,t2)} models heterogeneity with a constant plus power variance function. 
    Letting \eqn{v} denote the variance covariate and \eqn{\sigma^2(v)}{s2(v)}
    denote the variance function evaluated at \eqn{v}, the constant plus power
    variance function is defined as \eqn{\sigma^2(v) = (\theta_1 + |v|^{\theta_2})^2}{s2(v) = (t1 + |v|^t2)^2},
    where \eqn{\theta_1,\theta_2}{t1, t2} are the variance function coefficients.
    If \code{c("exp",t)},models heterogeneity with an
    exponential variance function. Letting \eqn{v} denote the variance covariate and \eqn{\sigma^2(v)}{s2(v)}
  denote the variance function evaluated at \eqn{v}, the exponential
  variance function is defined as \eqn{\sigma^2(v) = e^{2 * \theta * v}}{s2(v) = exp(2* t * v)}, where \eqn{\theta}{t} is the variance
  function coefficient.}
  \item{ftype}{character value "lmer", "lme" or "MCMCglmm" specifying the function to use to fit 
    the model. Actually "lmer" only is accepted }
  \item{mer.sim}{simulate the data using simulate.merMod from lme4. Faster for large sample size but not as flexible.}
}

\details{
 P-values for random effects are estimated using a log-likelihood ratio
 test between two models with and without the effect. Power represent
 the percentage of simulations providing a significant p-value for a
 given random structure
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

\seealso{\code{\link{EAMM}}, \code{\link{SSF}}, \code{\link{plot.PAMM}} 
 }

\examples{
\dontrun{
  ours=PAMM(numsim=10,group=c(seq(10,50,10),100),repl=c(2,4,6),
           randompart=c(0.4,0.1,0.5,0.1),fixed=c(0,1,0.7))       
  plot(ours,"both")
   }
}
\keyword{ misc }
