
create_blank_plate <- function(WellR=LETTERS[1:16],WellC=1:24) {
    ## create blank plate data frame
    plate <- expand.grid(WellR=WellR,WellC=WellC)
    plate$Well   <- with(plate,paste0(WellR,WellC))
    # plate$Sample <- NA
    # plate$Probe  <- NA
    return(plate)
}

create_colkey_6in24 <- function(...) {
    ## creates a 24-column key with 3+RT Techreps and 1-RT
    ## suitable for 6 pieces (samples or target/probe sets)
    ## example: create_colkey_6in24(Sample=LETTERS[1:6])
    
    colkey <- tibble(WellC=1:24,
                          Type=c(rep("+RT",18),rep("-RT",6)) %>% 
                              factor(levels=c("+RT","-RT")),
                          TechRep=rep(c(1,2,3,1),each=6) %>% 
                              factor(levels=1:3) )
    if( !missing(...) ) {
        pieces6 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces6) == 6)
        pieces24 <- bind_rows(pieces6,pieces6,pieces6,pieces6)
        colkey <- bind_cols(colkey, pieces24)
    }
    return(colkey)
}

create_rowkey_4in16 <- function(...) {
    ## creates a 16-row key suitable for 4 pieces (samples or target/probe sets)
    ## example: create_rowkey_4in16(Sample=c("me","you","them","him","her","dog","cat","monkey"))
    rowkey <- tibble(WellR=LETTERS[1:16],
                     Type=c(rep("+RT",12),rep("-RT",4)) %>% 
                         factor(levels=c("+RT","-RT")),
                     TechRep=rep(c(1,2,3,1),each=4) %>% 
                         factor(levels=1:3) )
    if( !missing(...) ) {
        pieces4 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces4) == 4)
        pieces16 <- bind_rows(pieces4,pieces4)
        rowkey <- bind_cols(rowkey, pieces16)
    }
    return(rowkey)
}

create_rowkey_8in16_plain <- function(...) {
    ## creates a 16-row key suitable for 8 pieces (samples or target/probe sets)
    ## example: create_rowkey_8in16(Sample=c("me","you","them","him","her","dog","cat","monkey"))
    rowkey <- tibble(WellR=LETTERS[1:16])
    if( !missing(...) ) {
        pieces8 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces8) == 8)
        pieces16 <- bind_rows(pieces8,pieces8)
        rowkey <- bind_cols(rowkey, pieces16)
    }
    return(rowkey)
}


label_plate_rowcol <- function(plate,rowkey=NULL,colkey=NULL) {
    ## label a plate by row and column keys
    if (!is.null(colkey)) {
        plate <- merge(plate,colkey,by="WellC")
    }
    if (!is.null(rowkey)) {
        plate <- merge(plate,rowkey,by="WellR")
    }
    return(plate[order(plate$WellR,plate$WellC),])
}

display_plate <- function(plate) {
    ggplot(data=plate,aes(x=factor(WellC),
                          y=factor(WellR,levels=rev(LETTERS)))) +
        geom_tile(aes(fill=Probe),alpha=0.3) +
        geom_text(aes(label=paste(Probe,Sample,Type,sep="\n")),size=2.5,lineheight=1) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        coord_equal() +
        theme_void() + 
        theme(axis.text=element_text(angle=0),
              panel.grid.major=element_blank(),
              legend.position="none",
              plot.margin=unit(rep(0.01,4),"npc"),
              panel.border = element_blank())
}


getNormCt <- function(df,mycolumn="Ct",normProbes="ALG9",probename="Probe") {
    ### function to take data frame and attach a column to normalize things by.
    
    # make subset of df where gene is one of normGenes
    subdf <- df[df[[probename]] %in% normProbes,]
    
    # assign median of mycolumn to df$normct
    # note this is the same value for every row, a waste of space technically
    df$norm.by <- median(subdf[[mycolumn]],na.rm=TRUE)
    return(df)
}


normalizeqPCR <- function(df,mycolumn="Ct",normProbes="ALG9",probename="Probe") {
    # make normed count, grouped by Sample (biological rep)
    dfout <- 
        group_by(df,Sample) %>% # group by Sample
        do(getNormCt(.,mycolumn,normProbes,probename)) %>%      # get norm value for each Sample
        ungroup()                 # combine/ungroup again
    
    # Assign normalized values by dividing by normby
    dfout$Value.norm <- dfout[[mycolumn]] - dfout$norm.by
    dfout$Value.normexp <- 2^-dfout$Value.norm
                             
    return(dfout)
}
