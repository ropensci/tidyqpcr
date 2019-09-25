## plot helper functions

scale_y_log2nice <- function(name=waiver(),omag=seq(-6,6),scilabels=FALSE,...) {
    # plot y axis on log2-scale with nice defaults
    breaks2 <- 2^omag
    if (scilabels) {
        labels2 <- paste("2^{",omag,"}",sep="")
    } else {
        labels2 <- breaks2
    }
    scale_y_log10(name,breaks=breaks2,labels=parse(text=labels2),...)
}

scale_y_log10nice <- function(name=waiver(),omag=seq(-6,6),scilabels=FALSE,...) {
    # plot y axis on log10-scale with nice defaults
    breaks10 <- 10^omag
    if (scilabels) {
        labels10 <- paste("10^{",omag,"}",sep="")
    } else {
        labels10 <- breaks10
    }
    scale_y_log10(name,breaks=breaks10,labels=parse(text=labels10),...)
}