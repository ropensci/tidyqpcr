
debaseline <- function(plateamp,maxcycle=10) {
    baseline <- 
        plateamp %>%
        group_by(Well) %>%
        filter(Program == 2, Cycle <= maxcycle) %>%
        summarize(Base = median(Fluor))
    plateamp %>% 
        left_join(baseline) %>%
        mutate(Signal=Fluor-Base)
        
}

getdRdT <- function(TT,RR,method=c("spline","diff"),...) {
    if (method == "diff") {
       return( -c(diff(RR)/diff(TT),NA) )
    } else if (method == "spline") {
        fit <- smooth.spline(x = TT, y=RR,...)
        return(-1 * predict(object = fit, x = TT, deriv = 1)$y )
    }
}

demelt <- function(platemelt) {
    platemelt %>%
        arrange(Well,Temperature) %>%
        group_by(Well) %>%
        mutate(dRdT=getdRdT(Temperature,Signal)) %>%
        ungroup()
}
