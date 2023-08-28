#' Graphic output of the PAMM function
#' 
#' Provide graphic interpretation of the simulation results
#' 
#' 
#' @param x an SSF object 
#' @param \dots potentially further arguments to pass to methods
#' 
#' \examples
#' \dontrun{
#'    oursSSF <- SSF(50,100,10,c(0.4,0.1,0.6,0))
#'    plot(oursSSF)
#'    }
#' 
#' 
#' @export


plot.SSF <- function(x, ...) {
  if (!inherits(x, "SSF"))
    stop("use only with \"SSF\" xs")

  xs <- c("int.pval", "sl.pval", "int.power", "sl.power")
  ys <- c("ipv", "slpv", "ipo", "slpo")
  labs <- c("int.p-value", "int.power", "slope.p-value", "slope.power")
  mains <- c("Intercept P-value", "Intercept Power", "Slope P-value", "Slope Power")
  par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 0))
  for (i in 1:4) {
    plot(x$nb.ID, x[, xs[i]], type = "l", xlab = "", ylab = "", ylim = c(0, 1),
      xaxt = "n", bty = "L")
    axis(1, at = x$nb.ID, labels = paste(x$nb.ID, x$nb.repl, sep = "/"), tick = FALSE)
    if (i <= 2)
      abline(h = 0.05)
    lines(x$nb.ID, x[, paste("CIup.", ys[i], sep = "")], lty = 2)
    lines(x$nb.ID, x[, paste("CIlow.", ys[i], sep = "")], lty = 2)
  }
  mtext("Nb Gr./Repl.", 1, line = 2, outer = TRUE, cex = 1.2)
  mtext(c("P-value", "Power"), 2, at = c(0.75, 0.25), line = 2, outer = TRUE, cex = 1.2)
  mtext(c("Intercept", "Slope"), 3, at = c(0.25, 0.75), line = 0.5, outer = TRUE,
    cex = 1.2)
}
