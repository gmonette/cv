# library(plotrix)

hexagon <- function (x, y, unitcell = 1, col = NA, border = "black", lwd=1)
{
  polygon(c(x, x, x + unitcell/2, x + unitcell, x + unitcell,
            x + unitcell/2), c(y + unitcell * 0.125, y + unitcell *
                                 0.875, y + unitcell * 1.125, y + unitcell * 0.875, y +
                                 unitcell * 0.125, y - unitcell * 0.125), col = col, border = border,
          lwd=lwd)
}


row1 <- seq(0, 0.8, length=5)*5
row2 <- rep(0.2, 5)*5
row3 <- seq(0.8, 0, length=5)*5
heights <- rbind(row1, row2, row3)


png("~/temp/hex.png", width=500, height=500)

par(mar=rep(c(2.5, 0, 3.5, 2.5)))

barplot(heights, xlim=c(-2.5, 7.5), ylim=c(-2.5, 7.5),
        col=c("gray", "blue", "gray"),
        axes=FALSE)

hexagon(-0.65, -1.1, 7.5, lwd=5, border="blue",
        col=scales::alpha("blue", 0.10))

text(3, 6, "cv", cex=4, col="blue")
text(3.1, -0.6, "cross-validation", col="blue", cex=2)

dev.off()
