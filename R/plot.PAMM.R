`plot.PAMM` <-
function (x, graphtype = "both", nbgroup, nbrepl, ...) 
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
      x11()
      plot(wireframe(int.pval~nb.ID+nb.repl,x,colorkey=FALSE,drape=TRUE,scales = list(arrows = FALSE),
      xlab="group",ylab="repl", main="Intercept P-value",
          screen = list(z =-50, x = -70, y = 0)))    
      x11() 
      plot(wireframe(int.power~nb.ID+nb.repl,x,colorkey=FALSE,drape=TRUE,scales = list(arrows = FALSE),
      xlab="group",ylab="repl", main="Intercept Power Calculations",
          screen = list(z =-50, x = -70, y = 0))) 
      x11()
      plot(wireframe(sl.pval~nb.ID+nb.repl,x,colorkey=FALSE,drape=TRUE,scales = list(arrows = FALSE),
      xlab="group",ylab="repl", main="Slope P-value",
          screen = list(z =-50, x = -70, y = 0)))
      x11()
      plot(wireframe(sl.power~nb.ID+nb.repl,x,colorkey=FALSE,drape=TRUE,scales = list(arrows = FALSE),
      xlab="group",ylab="repl", main="Slope Power Calculations",
          screen = list(z =-50, x = -70, y = 0)))
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
