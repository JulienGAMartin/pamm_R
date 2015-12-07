`plot.PAMM` <-
function (x, graphtype = "both", nbgroup, nbrepl, phi=25, theta=30, ticktype="detailed",nticks= 4, nbcol=100, color="grey",coltype="restricted",...) 
{
    if (!inherits(x, "PAMM")) 
        stop("use only with \"PAMM\" objects")
    if (graphtype == "group") {
        par(mfrow = c(2, 2))
        plot(x$nb.ID[x$nb.repl == nbrepl], x$int.pval[x$nb.repl == 
            nbrepl], type = "l", xlab = "Number of ID", ylab = "P-value", 
            main = "Intercept P-value", ylim = c(0, 1))
        abline(h = 0.05)
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIup.ipv[x$nb.repl == 
            nbrepl], lty = 2)
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIlow.ipv[x$nb.repl == 
            nbrepl], lty = 2)
        plot(x$nb.ID[x$nb.repl == nbrepl], x$sl.pval[x$nb.repl == 
            nbrepl], type = "l", xlab = "Number of ID", ylab = "P-value", 
            main = "Slope P-value", ylim = c(0, 1))
        abline(h = 0.05)
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIup.slpv[x$nb.repl == 
            nbrepl], lty = 2)
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIlow.slpv[x$nb.repl == 
            nbrepl], lty = 2)
        plot(x$nb.ID[x$nb.repl == nbrepl], x$int.power[x$nb.repl == 
            nbrepl], type = "l", ylim = c(0, 1), xlab = "Number of ID", 
            ylab = "Power", main = "Intercept Power Calculations")
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIup.ipo[x$nb.repl == 
            nbrepl], lty = 2)
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIlow.ipo[x$nb.repl == 
            nbrepl], lty = 2)
        plot(x$nb.ID[x$nb.repl == nbrepl], x$sl.power[x$nb.repl == 
            nbrepl], type = "l", ylim = c(0, 1), xlab = "Number of ID", 
            ylab = "Power", main = "Slope Power Calculations")
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIup.slpo[x$nb.repl == 
            nbrepl], lty = 2)
        lines(x$nb.ID[x$nb.repl == nbrepl], x$CIlow.slpo[x$nb.repl == 
            nbrepl], lty = 2)
    }
    if (graphtype == "repl") {
        par(mfrow = c(2, 2))
        plot(x$nb.repl[x$nb.ID == nbgroup], x$int.pval[x$nb.ID == 
            nbgroup], type = "l", ylim = c(0, 1), xlab = "Number of repl", 
            ylab = "P-value", main = "Intercept P-value")
        abline(h = 0.05)
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIup.ipv[x$nb.ID == 
            nbgroup], lty = 2)
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIlow.ipv[x$nb.ID == 
            nbgroup], lty = 2)
        plot(x$nb.repl[x$nb.ID == nbgroup], x$sl.pval[x$nb.ID == 
            nbgroup], type = "l", ylim = c(0, 1), xlab = "Number of repl", 
            ylab = "P-value", main = "Slope P-value")
        abline(h = 0.05)
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIup.slpv[x$nb.ID == 
            nbgroup], lty = 2)
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIlow.slpv[x$nb.ID == 
            nbgroup], lty = 2)
        plot(x$nb.repl[x$nb.ID == nbgroup], x$int.power[x$nb.ID == 
            nbgroup], type = "l", ylim = c(0, 1), xlab = "Number of repl", 
            ylab = "Power", main = "Intercept Power Calculations")
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIup.ipo[x$nb.ID == 
            nbgroup], lty = 2)
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIlow.ipo[x$nb.ID == 
            nbgroup], lty = 2)
        plot(x$nb.repl[x$nb.ID == nbgroup], x$sl.power[x$nb.ID == 
            nbgroup], type = "l", ylim = c(0, 1), xlab = "Number of repl", 
            ylab = "Power", main = "Slope Power Calculations")
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIup.slpo[x$nb.ID == 
            nbgroup], lty = 2)
        lines(x$nb.repl[x$nb.ID == nbgroup], x$CIlow.slpo[x$nb.ID == 
            nbgroup], lty = 2)
    }
    if (graphtype == "both") {

      par(mfrow = c(2, 2),mar=c(1.5,1.5,1.5,1.5))
      X <- unique(x$nb.ID)
      Y <- unique(x$nb.repl)
      z <- c("int.pval", "int.power","sl.pval","sl.power")
      main <- c("Intercept P-value","Intercept Power","Slope P-value","Slope Power")
      if (color == "grey") color <- grey(0:nbcol/nbcol)
      else if (color == "cm") color <- cm.colors(nbcol)
      else if (color == "rainbow") color <- rainbow(nbcol)
      
      for (i in 1:4) {
     	Z <- matrix(x[,z[i]], nrow = length(X), ncol = length(Y))
     	nrz <- nrow(Z)
      	ncz <- ncol(Z)
      	zfacet <- (Z[-1, -1] + Z[-1, -ncz] + Z[-nrz, -1] + Z[-nrz, -ncz]) / 4
     	if (coltype == "0-1") facetcol <- facetcol <- cut(zfacet, seq(0,1,1/nbcol))
     	else if (coltype == "restricted") facetcol <- cut(zfacet,nbcol)
      	persp(X, Y, Z, 
            col = color[facetcol], zlim = c(0, 1),
            phi=phi, theta=theta, ticktype=ticktype, nticks=nticks,
            xlab = "group", ylab = "repl", zlab = "", main=main[i]) 
      }
    }
    if (graphtype == "both.dyn") {
        open3d()
        bg3d("white")
        material3d(col = "white")
        persp3d(unique(x$nb.ID), unique(x$nb.repl), matrix(x$int.pval, 
            nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
            col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "group", 
            ylab = "repl", zlab = "int.p-value")
        open3d()
        bg3d("white")
        material3d(col = "black")
        persp3d(unique(x$nb.ID), unique(x$nb.repl), matrix(x$int.power, 
            nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
            col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "group", 
            ylab = "repl", zlab = "int.power")
        open3d()
        bg3d("white")
        material3d(col = "black")
        persp3d(unique(x$nb.ID), unique(x$nb.repl), matrix(x$sl.pval, 
            nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
            col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "group", 
            ylab = "repl", zlab = "sl.p-value")
        open3d()
        bg3d("white")
        material3d(col = "black")
        persp3d(unique(x$nb.ID), unique(x$nb.repl), matrix(x$sl.power, 
            nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
            col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "group", 
            ylab = "repl", zlab = "sl.power")
    }
}
