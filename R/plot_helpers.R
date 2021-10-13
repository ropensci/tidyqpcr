#' Nice tick labels for logarithmic axes in ggplot2.
#' 
#' These functions are wrappers for scale_x_log10, etc., that specify
#' axis tick labels as powers of 10 and powers of 2 with default
#' nice format or easy scientific notation.
#' 
#' @param name axis name.
#' @param omag orders of magnitude or fold-changes for axis labels, 
#' usually an integer sequence.
#' @param scilabels display labels in scientific format, e.g. \eqn{10^2} vs 100.
#' @param ... other arguments to continuous_scale().
#' @return a ggproto object as output by continuous_scale(). 
#' @name log_plot_helpers
NULL

#' @describeIn log_plot_helpers plot x axis on log2-scale with nice defaults
#' @export
#' 
#' @examples 
#' library(ggplot2)
#' 
#' # create example plot with ggplot2 dataset 
#' p1 <- ggplot(mpg, aes(displ, hwy)) +
#' geom_point()
#' 
#' p1 + scale_x_log2nice()
#' 
#' p1 + scale_y_log10nice()
#' 
scale_x_log2nice <- function(name=ggplot2::waiver(),
                             omag=seq(-10,10),
                             scilabels=FALSE,
                             ...) {
    breaks2 <- 2^omag
    if (scilabels) {
        labels2 <- paste("2^{",omag,"}",sep="")
    } else {
        labels2 <- breaks2
    }
    ggplot2::scale_x_log10(name,breaks=breaks2,labels=parse(text=labels2),...)
}

#' @describeIn log_plot_helpers plot x axis on log10-scale with nice defaults
#' @export
#' 
scale_x_log10nice <- function(name=ggplot2::waiver(),omag=seq(-10,10),scilabels=FALSE,...) {
    # @ewallace: ideally would also create minor breaks from 2:9.
    breaks10 <- 10^omag
    if (scilabels) {
        labels10 <- paste("10^{",omag,"}",sep="")
    } else {
        labels10 <- breaks10
    }
    ggplot2::scale_x_log10(name,breaks=breaks10,labels=parse(text=labels10),...)
}

#' @describeIn log_plot_helpers plot y axis on log2-scale with nice defaults
#' @export
#' 
scale_y_log2nice <- function(name=ggplot2::waiver(),omag=seq(-10,10),scilabels=FALSE,...) {
    breaks2 <- 2^omag
    if (scilabels) {
        labels2 <- paste("2^{",omag,"}",sep="")
    } else {
        labels2 <- breaks2
    }
    ggplot2::scale_y_log10(name,breaks=breaks2,labels=parse(text=labels2),...)
}

#' @describeIn log_plot_helpers plot y axis on log10-scale with nice defaults
#' @export
#' 
scale_y_log10nice <- function(name=ggplot2::waiver(),omag=seq(-10,10),scilabels=FALSE,...) {
    # @ewallace: ideally would also create minor breaks from 2:9.
    breaks10 <- 10^omag
    if (scilabels) {
        labels10 <- paste("10^{",omag,"}",sep="")
    } else {
        labels10 <- breaks10
    }
    ggplot2::scale_y_log10(name,breaks=breaks10,labels=parse(text=labels10),...)
}

#' @describeIn log_plot_helpers plot x AND y axes on log10-scale with nice defaults
#' @export
#' 
scale_loglog10 <- function(...) {
    list(scale_x_log10nice(...),scale_y_log10nice(...))
}