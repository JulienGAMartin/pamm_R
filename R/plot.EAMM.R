#' Graphic output of the EAMM function 
#' 
#' Provide graphic interpretation of the simulation results
#' 
#'  @param x  an EAMM object
#'  @param graphtype  "VI", "VS", or "both"
#'     "VI" give graphs with varying variance component of intercept and
#'     with a fixed variance component for slope specified in vs argument
#'     "VS" give graphs with varying variance component for slope and with a
#'     fixed variance component of intercept specified in vi argument
#'     "both" 3-D plot see also fun3D argument
#'     
#'  @param vi  VI for which plots the output. Necessary for "VS" type of graph 
#'  @param vs  VS for which plots the output. Necessary for "VI" type of graph 
#'  @param fun3D plot function used to plot the 3D graph.
#'   	"wireframe" uses lattice,
#'   	"persp" uses  graphics,
#'   	and "open3d" uses rgl.
#'  @param \dots  potentially further arguments to pass to methods 
#' 
#' 
#' @examples
#' \dontrun{
#'   ours <- EAMM(numsim=10, group=10, repl=4,
#'                VI=seq(0.1,0.95,0.05), VS=c(0.05,0.1) )
#'   plot(ours, "both")
#'   plot(ours, "VI",vs=0.1)
#'   plot(ours,"VS",vi=0.2)
#'    }
#' 
#' @export


plot.EAMM <- function(x, graphtype = "both", vi, vs, fun3D = "wireframe", ...) {
  if (!inherits(x, "EAMM"))
    stop("use only with \"EAMM\" objects")

  extra <- list(...)  #not doing anything for the moment

  xs <- c("int.pval", "sl.pval", "int.power", "sl.power")
  ys <- c("ipv", "slpv", "ipo", "slpo")
  labs <- c("int.p-value", "slope.p-value", "int.power", "slope.power")
  mains <- c("Intercept P-value", "Slope P-value", "Intercept Power", "Slope Power")

  if (graphtype == "VI") {
    par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 0))
    for (i in 1:4) {
      plot(x$VI[x$VS == vs], x[x$VS == vs, xs[i]], type = "l", xlab = "", ylab = "",
        ylim = c(0, 1), bty = "L")
      if (i <= 2)
        abline(h = 0.05)
      lines(x$VI[x$VS == vs], x[x$VS == vs, paste("CIup.", ys[i], sep = "")],
        lty = 2)
      lines(x$VI[x$VS == vs], x[x$VS == vs, paste("CIlow.", ys[i], sep = "")],
        lty = 2)
    }
    mtext(expression(paste("Intercept Variance ", V[I])), 1, line = 2, outer = TRUE,
      cex = 1.2)
    mtext(c("P-value", "Power"), 2, at = c(0.75, 0.25), line = 2, outer = TRUE,
      cex = 1.2)
    mtext(c("Intercept", "Slope"), 3, at = c(0.25, 0.75), line = 0.5, outer = TRUE,
      cex = 1.2)
  }

  if (graphtype == "VS") {
    par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 0))
    for (i in 1:4) {
      plot(x$VS[x$VI == vi], x[x$VI == vi, xs[i]], type = "l", xlab = "", ylab = "",
        ylim = c(0, 1), bty = "L")
      if (i <= 2)
        abline(h = 0.05)
      lines(x$VS[x$VI == vi], x[x$VI == vi, paste("CIup.", ys[i], sep = "")],
        lty = 2)
      lines(x$VS[x$VI == vi], x[x$VI == vi, paste("CIlow.", ys[i], sep = "")],
        lty = 2)
    }
    mtext(expression(paste("Slope Variance ", V[S])), 1, line = 2, outer = TRUE,
      cex = 1.2)
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
        p <- wireframe(as.formula(paste(xs[i], " ~ VI + VS")), x, colorkey = FALSE,
          drape = TRUE, scales = list(arrows = FALSE, distance = c(2, 2,
          2)), xlab = "Var.Intercept", ylab = "Var.Slope", main = mains[i],
          zlab = list(zl[i], rot = 90), screen = list(z = -50, x = -70, y = 0),
          par.settings = par.set)
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
      X <- unique(x$VI)
      Y <- unique(x$VS)
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
          theta = theta, ticktype = ticktype, nticks = nticks, xlab = expression(V[I]),
          ylab = expression(V[S]), zlab = "")
      }
      mtext(c("P-value", "Power"), 2, at = c(0.75, 0.25), line = 2, outer = TRUE,
        cex = 1.2)
      mtext(c("Intercept", "Slope"), 3, at = c(0.25, 0.75), line = 0.5, outer = TRUE,
        cex = 1.2)
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
          rgl::persp3d(unique(x$VI), unique(x$VS), matrix(x[, xs[i]], nrow = length(unique(x$VI)),
          ncol = length(unique(x$VS))), col = rainbow(10), box = FALSE,
          zlim = c(0, 1), xlab = "VI", ylab = "VS", zlab = labs[i])
        }
      }
    }
  }
}
