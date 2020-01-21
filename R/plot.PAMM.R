plot.PAMM <- function(x, graphtype = "both", nbgroup, nbrepl, fun3D = "wireframe",
  ...) {

  if (!inherits(x, "PAMM"))
    stop("use only with \"PAMM\" objects")

  extra <- list(...)  #not doing anything for the moment
  
  xs <- c("int.pval", "sl.pval", "int.power", "sl.power")
  ys <- c("ipv", "slpv", "ipo", "slpo")
  labs <- c("int.p-value", "slope.p-value", "int.power", "slope.power")
  mains <- c("Intercept P-value", "Slope P-value", "Intercept Power", "Slope Power")

  if (graphtype == "group") {
    par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 0))
    for (i in 1:4) {

      plot(x$nb.ID[x$nb.repl == nbrepl], x[x$nb.repl == nbrepl, xs[i]], type = "l",
        xlab = "", ylab = "", ylim = c(0, 1), bty = "L")
      if (i <= 2)
        abline(h = 0.05)
      lines(x$nb.ID[x$nb.repl == nbrepl], x[x$nb.repl == nbrepl, paste("CIup.",
        ys[i], sep = "")], lty = 2)
      lines(x$nb.ID[x$nb.repl == nbrepl], x[x$nb.repl == nbrepl, paste("CIlow.",
        ys[i], sep = "")], lty = 2)
    }
    mtext("Number of group", 1, line = 2, outer = TRUE, cex = 1.2)
    mtext(c("P-value", "Power"), 2, at = c(0.75, 0.25), line = 2, outer = TRUE,
      cex = 1.2)
    mtext(c("Intercept", "Slope"), 3, at = c(0.25, 0.75), line = 0.5, outer = TRUE,
      cex = 1.2)
  }
  if (graphtype == "repl") {
    par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 0))
    for (i in 1:4) {
      plot(x$nb.repl[x$nb.ID == nbgroup], x[x$nb.ID == nbgroup, xs[i]], type = "l",
        xlab = "", ylab = "", ylim = c(0, 1), bty = "L")
      if (i <= 2)
        abline(h = 0.05)
      lines(x$nb.repl[x$nb.ID == nbgroup], x[x$nb.ID == nbgroup, paste("CIup.",
        ys[i], sep = "")], lty = 2)
      lines(x$nb.repl[x$nb.ID == nbgroup], x[x$nb.ID == nbgroup, paste("CIlow.",
        ys[i], sep = "")], lty = 2)
    }
    mtext("Number of replicates", 1, line = 2, outer = TRUE, cex = 1.2)
    mtext(c("P-value", "Power"), 2, at = c(0.75, 0.25), line = 2, outer = TRUE,
      cex = 1.2)
    mtext(c("Intercept", "Slope"), 3, at = c(0.25, 0.75), line = 0.5, outer = TRUE,
      cex = 1.2)
  }
  if (graphtype == "both") {
    if (fun3D == "wireframe") {
      par.set <- list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
      zl <- c("P-value", "Power", "P-value", "Power")
      for (i in 1:4) {
        p <- wireframe(as.formula(paste(xs[i], " ~ nb.ID + nb.repl")), x,
          colorkey = FALSE, drape = TRUE, scales = list(arrows = FALSE, distance = c(2,
          2, 2)), xlab = "Group", ylab = "Repl", main = mains[i], zlab = list(zl[i],
          rot = 90), screen = list(z = -50, x = -70, y = 0), par.settings = par.set)
        if (i == 1)
          print(p, split = c(1, 1, 2, 2), more = TRUE) else if (i == 2)
          print(p, split = c(1, 2, 2, 2), more = TRUE) else if (i == 3)
          print(p, split = c(2, 1, 2, 2), more = TRUE) else if (i == 4)
          print(p, split = c(2, 2, 2, 2))
      }
    }
    if (fun3D == "persp") {
      phi <- ifelse(is.null(extra$phi), 25, extra$phi)
      theta <- ifelse(is.null(extra$theta), 30, extra$theta)
      ticktype <- ifelse(is.null(extra$ticktype), "detailed", extra$ticktype)
      nticks <- ifelse(is.null(extra$nticks), 4, extra$nticks)
      nbcol <- ifelse(is.null(extra$nbcol), 100, extra$nbcol)
      color <- ifelse(is.null(extra$color), "grey", extra$color)
      coltype <- ifelse(is.null(extra$coltype), "restricted", extra$coltype)
      par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(0, 4, 2, 0))
      X <- unique(x$nb.ID)
      Y <- unique(x$nb.repl)

      if (color == "grey")
        color <- grey(0:nbcol/nbcol) else if (color == "cm")
        color <- cm.colors(nbcol) else if (color == "rainbow")
        color <- rainbow(nbcol)

      for (i in 1:4) {
        Z <- matrix(x[, xs[i]], nrow = length(X), ncol = length(Y))
        nrz <- nrow(Z)
        ncz <- ncol(Z)
        zfacet <- (Z[-1, -1] + Z[-1, -ncz] + Z[-nrz, -1] + Z[-nrz, -ncz])/4
        if (coltype == "0-1")
          facetcol <- facetcol <- cut(zfacet, seq(0, 1, 1/nbcol)) else if (coltype == "restricted")
          facetcol <- cut(zfacet, nbcol)
        persp(X, Y, Z, col = color[facetcol], zlim = c(0, 1), phi = phi,
          theta = theta, ticktype = ticktype, nticks = nticks, xlab = "group",
          ylab = "repl", zlab = "")
        mtext(c("P-value", "Power"), 2, at = c(0.75, 0.25), line = 2, outer = TRUE,
          cex = 1.2)
        mtext(c("Intercept", "Slope"), 3, at = c(0.25, 0.75), line = 0.5,
          outer = TRUE, cex = 1.2)
      }
    }


    if (fun3D == "open3d") {
      if (!requireNamespace("rgl", quietly = TRUE)) {
        warning("rgl package needed for this function to work. Please install it. ",
          call. = FALSE)
      } else {
        for (i in 1:4) {
          rgl::open3d()
          rgl::bg3d("white")
          rgl::material3d(col = "white")
          rgl::persp3d(unique(x$nb.ID), unique(x$nb.repl), matrix(x[, xs[i]],
          nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))),
          col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "Nb group",
          ylab = "Nb repl", zlab = labs[i])
        }
      }
    }
  }
}
