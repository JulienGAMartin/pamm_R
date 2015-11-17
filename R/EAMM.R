`EAMM` <-
function (numsim, group, repl, fixed = c(0, 1, 0), VI = seq(0.05, 
    0.95, 0.05), VS = seq(0.05, 0.5, 0.05), CoIS = 0, relIS = "cor", n.X = NA, autocorr.X = 0,
          X.dist = "gaussian", intercept=0,heteroscedasticity=c("null")) 
{
    Hetero <- heteroscedasticity[[1]]
    het <- as.numeric(heteroscedasticity[-1])

    if (X.dist=="gaussian") {
        FM <- fixed[[1]]
        FV <- fixed[[2]]
        FE <- fixed[[3]]
    }
    if (X.dist=="unif") {
        Xmin <- fixed[[1]]
        Xmax <- fixed[[2]]
        FE <- fixed[[3]]
    }

    vgi <- numeric(length(VI) * length(VS))
    vgs <- numeric(length(VI) * length(VS))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(VI) * length(VS))
    slpowCIlower <- numeric(length(VI) * length(VS))
    slpowCIupper <- numeric(length(VI) * length(VS))
    slpvalestimate <- numeric(length(VI) * length(VS))
    slpvalCIlower <- numeric(length(VI) * length(VS))
    slpvalCIupper <- numeric(length(VI) * length(VS))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(VI) * length(VS))
    intpowCIlower <- numeric(length(VI) * length(VS))
    intpowCIupper <- numeric(length(VI) * length(VS))
    intpvalestimate <- numeric(length(VI) * length(VS))
    intpvalCIlower <- numeric(length(VI) * length(VS))
    intpvalCIupper <- numeric(length(VI) * length(VS))
    kk <- 0
    for (k in VI) {
        for (r in VS) {
            N <- group * repl
            n.x <- ifelse( is.na(n.X)==TRUE, N, n.X)
            VR <- 1 - k
            if (VR >= 0) {
                for (i in 1:numsim) {
                    if (relIS == "cor") { CovIS <- CoIS * sqrt(k) * sqrt(r) }
                    if (relIS == "cov") { CovIS <- CoIS }
                    sigma <- matrix(c(k, CovIS, CovIS, r), ncol = 2)

                    if (X.dist=="gaussian"){
                        if (autocorr.X==0) { ef <- rnorm(n.x, FM, sqrt(FV)) }
                        else {
                            y <- numeric(n.x)
                            phi <- autocorr.X
                            y[1] <- rnorm(1, 0, sd = sqrt(FV))
                            for (t in 2:n.x) { y[t] <- rnorm(1, y[t-1]*phi, sd = sqrt(FV)) }
                            ef <- y+FM
                        }
                    }

                    if (X.dist=="unif"){
                        if (autocorr.X==0) { ef <- runif(n.x, Xmin, Xmax) }
                        else { stop("autocorrelation in fixed effects is not yet implemented for uniform distribution") }
                    }

                    if (n.x!=N) {
                        if (n.x>=repl) {
                            inief <- sample(1:(n.x-repl+1),group,replace=TRUE)
                            EFrk <- rep(inief,repl) + rep (0:(repl-1),each=group)
                            EF <- ef[EFrk]
                        }
                        if (n.x<repl) { 
                            EF <- numeric(N)
                            EF[1:(n.x*group)] <- rep(ef,each=group)
                            EF[(n.x*group+1):N] <- sample(ef,length((n.x*group+1):N),replace=TRUE)
                        }
                    }
                    else { EF <- ef }

                    er <- numeric(length(N))
                    if (Hetero=="null") (er <- rnorm(N, intercept, sqrt(VR)))
                    if (Hetero=="power") (
                    for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*(het[1]+abs(EF[n])^het[2])^2))} )
                    if (Hetero=="exp")  (
                    for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*exp(2*het[1]*EF[n])))} )

                    X <- sort(rep(c(1:repl), group))
                    db <- data.frame(ID = rep(1:group, repl), obs = 1:N, 
                      error = er, X = X, EF = EF)
                    x <- rmvnorm(group, c(0, 0), sigma, method = "chol")
                    db$rand.int <- rep(x[, 1], repl)
                    db$rand.sl <- rep(x[, 2], repl)
                    db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + 
                      db$error
                    
#models
                   if (r > 0) {
                      m.full <- lmer(Y ~ EF + (EF | ID), data = db)
                      m.nocov <- lmer(Y ~ EF + (1 | ID) + (0 + EF | ID), data = db) 
                      m.nosl <- lmer(Y ~ EF + (1 | ID), data = db)
                      m.noint <- lmer(Y ~ EF + (0 + EF | ID), data = db) 
                     
                      #anosl <- anova(m.nocov, m.nosl, refit =FALSE)
                      #powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
                      #pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
                      
                      anoIxE <- anova(m.full, m.nosl, refit =FALSE)
                      powersl[i] <- anoIxE[2, "Pr(>Chisq)"] <= 0.05
                      pvalsl[i] <- anoIxE[2, "Pr(>Chisq)"]
                      
                      anoint <- anova(m.nocov, m.noint, refit =FALSE)
                      powerint[i] <- anoint[2, "Pr(>Chisq)"] <= 0.05
                      pvalint[i] <- anoint[2, "Pr(>Chisq)"]
                   }
                   else {
                     powersl[i] <- 0
                     pvalsl[i] <- 1
                   
                    db$ dummy <- rep(1:2,nrow(db)/2)
                    m.lmer <- lmer(Y ~ EF + (1 | dummy), data = db)
                    m1.lmer <- lmer(Y ~ EF + (1 | ID), data = db) 
                    pvint <- pchisq(-2 * (logLik(m.lmer, REML = TRUE) - 
                      logLik(m1.lmer, REML = TRUE))[[1]], 1, lower.tail = FALSE)
                    powerint[i] <- pvint <= 0.05
                    pvalint[i] <- pvint
                    }
                }

                  kk <- kk + 1
                  vgi[kk] <- k
                  vgs[kk] <- r
                  slCIpow <- ci(powersl)
                  slpowestimate[kk] <- slCIpow["Estimate"]
                  slpowCIlower[kk] <- slCIpow["CI lower"]
                  slpowCIupper[kk] <- slCIpow["CI upper"]
                  slCIpval <- ci(pvalsl)
                  slpvalestimate[kk] <- slCIpval["Estimate"]
                  slpvalCIlower[kk] <- slCIpval["CI lower"]
                  slpvalCIupper[kk] <- slCIpval["CI upper"]
                  intCIpow <- ci(powerint)
                  intpowestimate[kk] <- intCIpow["Estimate"]
                  intpowCIlower[kk] <- intCIpow["CI lower"]
                  intpowCIupper[kk] <- intCIpow["CI upper"]
                  intCIpval <- ci(pvalint)
                  intpvalestimate[kk] <- intCIpval["Estimate"]
                  intpvalCIlower[kk] <- intCIpval["CI lower"]
                  intpvalCIupper[kk] <- intCIpval["CI upper"]
              }
              else {
                kk <- kk + 1
                vgi[kk] <- k
                vgs[kk] <- r
                slpowestimate[kk] <- NA
                slpowCIlower[kk] <- NA
                slpowCIupper[kk] <- NA
                slpvalestimate[kk] <- NA
                slpvalCIlower[kk] <- NA
                slpvalCIupper[kk] <- NA
                intpowestimate[kk] <- NA
                intpowCIlower[kk] <- NA
                intpowCIupper[kk] <- NA
                intpvalestimate[kk] <- NA
                intpvalCIlower[kk] <- NA
                intpvalCIupper[kk] <- NA
            }
        }
    }
    sim.sum <- data.frame(nb.ID = rep(group, length(VI) * length(VS)), 
        nb.repl = rep(repl, length(VI) * length(VS)), VI = vgi, 
        VS = vgs, int.pval = intpvalestimate, CIlow.ipv = intpvalCIlower, 
        CIup.ipv = intpvalCIupper, int.power = intpowestimate, 
        CIlow.ipo = intpowCIlower, CIup.ipo = intpowCIupper, 
        sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
        CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
        CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper)
    class(sim.sum) <- c("EAMM", "data.frame")
    sim.sum
}
