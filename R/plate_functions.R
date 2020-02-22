
#' Create a blank plate template as a tibble
#' 
#' @param WellR Vector of Row labels, usually LETTERS
#' @param WellC Vector of Column labels, usually numbers
#' @return tibble (data frame) with columns WellR, WellC, Well. This contains
#'   all pairwise combinations of WellR and WellC, as well as individual Well
#'   names. Both WellR and WellC are coerced to factors (even if WellC
#'   is supplied as numbers), to ensure order is consistent.
#'   
#'   However, Well is a character vector as that is the default behaviour of 
#'   "unite", and display order doesn't matter.
#'   
#'   Default value describes a full 384-well plate.
#'   
#' @examples
#' create_blank_plate(WellR=LETTERS[1:2],WellC=1:3)
#' create_blank_plate(WellR=LETTERS[1:8],WellC=1:12)
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom forcats as_factor
#' 
create_blank_plate <- function(WellR=LETTERS[1:16],WellC=1:24) {
    plate <- tidyr::crossing(WellR=as_factor(WellR),
                             WellC=as_factor(WellC)) %>%
        as_tibble() %>%
        tidyr::unite(Well,WellR,WellC,sep="",remove=FALSE)
    return(plate)
}

#' Create a 6-value, 24-column key for plates
#' 
#' Create a 24-column key with 6 values repeated over 24 plate columns. 
#' Each of the 6 values is repeated over 3x +RT Techreps and 1x -RT.
#' 
#' @param ... Vectors of length 6 describing well contents, e.g. sample or probe. 
#' @return tibble (data frame) with 24 rows, and columns WellC, Type, TechRep, and supplied values. 
#' @examples
#' create_colkey_6in24(Sample=LETTERS[1:6])
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#' 
create_colkey_6in24 <- function(...) {
    colkey <- tibble(WellC=factor(1:24),
                     Type=c(rep("+RT",18),rep("-RT",6)) %>% 
                         factor(levels=c("+RT","-RT")),
                     TechRep=rep(c(1,2,3,1),each=6) %>% 
                         factor(levels=1:3) )
    if( !missing(...) ) {
        pieces6 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces6) == 6)
        pieces24 <- dplyr::bind_rows(pieces6,pieces6,pieces6,pieces6)
        colkey <- dplyr::bind_cols(colkey, pieces24)
    }
    return(colkey)
}

#' Create a 4-dilution column key for primer calibration
#'
#' Creates a 24-column key for primer calibration, with 2x BioReps and 2x
#' TechReps, and 5-fold dilution until 5^4 of +RT; then -RT, NT controls. That
#' is a total of 6 versions of each sample replicate.
#'
#' @param Dilution Numeric vector of length 6 describing sample dilutions
#' @param DilutionNice Character vector of length 6 with nice labels for sample
#'   dilutions
#' @param Type Character vector of length 6 describing type of sample (+RT, -RT,
#'   NT)
#' @param BioRep Character vector of length 6 describing biological replicates
#' @param TechRep Character vector of length 6 describing technical replicates
#' @return tibble (data frame) with 24 rows, and columns WellC, Dilution,
#'   DilutionNice, Type, BioRep, TechRep.
#' @examples
#' create_colkey_4dilutions_mRTNT_in24()
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble
#' 
create_colkey_4dilutions_mRTNT_in24 <- function(
                     Dilution=c(1,1/5,1/25,1/125,1,1),
                     DilutionNice=c("1x","5x","25x","125x","-RT","NT"),
                     Type=c(rep("+RT",4),"-RT","NT"),
                     BioRep=rep(c("A","B"),each=12,length.out=24),
                     TechRep = rep(1:2,each=6,length.out=24)) {
    colkey <- tibble(WellC=factor(1:24),
                     Dilution=rep(Dilution,4),
                     DilutionNice=rep(DilutionNice,4),
                     Type=rep(Type,4) %>%
                         factor(levels=c("+RT","-RT","NT")),
                     BioRep=factor(BioRep),
                     TechRep=factor(TechRep))
    return(colkey)
}

#' Create a 6-dilution column key for primer calibration
#'
#' Creates a 24-column key for primer calibration, with 1x BioReps and 3x
#' TechReps, and 5-fold dilution until 5^6 of +RT; then -RT, NT controls. That
#' is a total of 8 versions of each replicate.
#'
#' @param Dilution Numeric vector of length 8 describing sample dilutions
#' @param DilutionNice Character vector of length 8 with nice labels for sample
#'   dilutions
#' @param Type Character vector of length 8 describing type of sample (+RT, -RT,
#'   NT)
#' @param TechRep Character vector of length 8 describing technical replicates
#' @return tibble (data frame) with 24 rows, and variables WellC, Dilution,
#'   DilutionNice, Type, BioRep, TechRep.
#' @examples
#' create_colkey_6dilutions_mRTNT_in24()
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble
#' 
create_colkey_6dilutions_mRTNT_in24 <- function(
                     Dilution=c(5^{0:-5},1,1),
                     DilutionNice=c("1x","5x","25x","125x",
                                    "625x","3125x","-RT","NT"),
                     Type=c(rep("+RT",6),"-RT","NT"),
                     TechRep = rep(1:3,each=8,length.out=24)) {
    colkey <- tibble(WellC=factor(1:24),
                     Dilution=rep(Dilution,3),
                     DilutionNice=rep(DilutionNice,3),
                     Type=rep(Type,3) %>%
                         factor(levels=c("+RT","-RT","NT")),
                     TechRep=factor(TechRep))
    return(colkey)
}

#' Create a 4-value, 16-row key for plates
#'
#' Create a 16-row key with 4 values repeated over 16 plate rows. Each of the 4
#' values is repeated over 3x +RT Techreps and 1x -RT.
#'
#' @param ... Vectors of length 4 describing well contents, e.g. sample or
#'   probe.
#' @return tibble (data frame) with 16 rows, and variables WellR, Type, TechRep,
#'   and supplied values.
#' @examples
#' create_rowkey_4in16(Sample=c("sheep","goat","cow","chicken"))
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#' 
create_rowkey_4in16 <- function(...) {
    rowkey <- tibble(WellR=LETTERS[1:16],
                     Type=c(rep("+RT",12),rep("-RT",4)) %>% 
                         factor(levels=c("+RT","-RT")),
                     TechRep=rep(c(1,2,3,1),each=4) %>% 
                         factor(levels=1:3) )
    if( !missing(...) ) {
        pieces4 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces4) == 4)
        pieces16 <- dplyr::bind_rows(pieces4,pieces4,pieces4,pieces4)
        rowkey <- dplyr::bind_cols(rowkey, pieces16)
    }
    return(rowkey)
}

#' Create a plain 8-value, 16-row key for plates
#'
#' Create a 16-row key with 8 values repeated over 16 plate rows. No other
#' information is included by default, hence "plain".
#'
#' @param ... Vectors of length 8 describing well contents, e.g. sample or
#'   probe.
#' @return tibble (data frame) with 16 rows, and variables WellC, and supplied
#'   values.
#' @examples
#' create_rowkey_8in16_plain(Sample=c("me","you","them","him",
#'                                    "her","dog","cat","monkey"))
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' 
create_rowkey_8in16_plain <- function(...) {
    rowkey <- tibble(WellR=LETTERS[1:16])
    if( !missing(...) ) {
        pieces8 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces8) == 8)
        pieces16 <- dplyr::bind_rows(pieces8,pieces8)
        rowkey <- dplyr::bind_cols(rowkey, pieces16)
    }
    return(rowkey)
}

#' Label a plate with sample and probe information
#'
#' See vignettes for further examples
#'
#' @param plate tibble (data frame) with variables WellR, WellC, Well. This
#'   would usually be produced by create_blank_plate(). It is possible to
#'   include other information in additional variables
#' @param rowkey tibble (data frame) describing plate rows, with variables WellR
#'   and others.
#' @param colkey tibble (data frame) describing plate columns, with variables
#'   WellC and others.
#' @return tibble (data frame) with variables WellR, WellC, Well. This contains
#'   all combinations of WellR and WellC found in the input plate, and all
#'   information supplied in rowkey and colkey distributed across every 
#'   well of the plate. Return plate is ordered by row WellR then column WellC. 
#'   Note this may cause a problem if WellC is a character (1,10,11,...), 
#'   instead of a factor or integer (1,2,3,...)
#' @examples
#' label_plate_rowcol(plate = create_blank_plate()) # returns blank plate
#' @family plate creation functions
#' 
#' @export
#' 
label_plate_rowcol <- function(plate,rowkey=NULL,colkey=NULL) {
    if (!is.null(colkey)) {
        assertthat::assert_that(has_name(colkey,"WellC"))
        # Note: should this if clause be a freestanding function?
        # coerce_column_to_factor(df, col, warn=TRUE) ?
        if( !is.factor(colkey$WellC) ) {
            warning("coercing WellC to a factor")
            colkey <- dplyr::mutate(colkey,WellC=as_factor(WellC))
        }
        plate <- dplyr::left_join(plate,colkey,by="WellC")
    }
    if (!is.null(rowkey)) {
        assertthat::assert_that(has_name(rowkey,"WellR"))
        if( !is.factor(rowkey$WellR) ) {
            warning("coercing WellR to a factor")
            rowkey <- dplyr::mutate(rowkey,WellR=as_factor(WellR))
        }
        plate <- dplyr::left_join(plate,rowkey,by="WellR")
    }
    return( dplyr::arrange( plate, WellR, WellC ) )
}


#' Display plate plan with Sample, Probe, Type per Well
#'
#' @param plate tibble with variables WellC, WellR, Sample, Probe, Type. 
#'   Output from label_plate_rowcol. 
#'
#' @return ggplot object; major output is to plot it
#'
#' @examples # !Needs a labeled example from label_plate_rowcol...
#' @family plate creation functions
#' 
#' @export
#' @importFrom forcats as_factor
#' 
display_plate <- function(plate) {
    rowlevels <- plate %>%
        pull(WellR) %>%
        as_factor %>%
        levels
                        
    ggplot2::ggplot(data=plate,
                    aes(x=as_factor(WellC),
                        y=as_factor(WellR))) +
        ggplot2::geom_tile(aes(fill=Probe),alpha=0.3) +
        ggplot2::geom_text(aes(label=paste(Probe,Sample,Type,sep="\n")),
                           size=2.5,lineheight=1) +
        ggplot2::scale_x_discrete(expand=c(0,0)) +
        ggplot2::scale_y_discrete(expand=c(0,0),
                                  limits=rev(rowlevels)) +
        ggplot2::coord_equal() +
        ggplot2::theme_void() + 
        ggplot2::theme(axis.text=ggplot2::element_text(angle=0),
                       panel.grid.major=ggplot2::element_blank(),
                       legend.position="none",
                       plot.margin=grid::unit(rep(0.01,4),"npc"),
                       panel.border=ggplot2::element_blank())
}


#' @describeIn normalizeqPCR get the median value of a set of normalization
#'   (reference) probes, for each sample.
#' 
#' @export
#' @importFrom magrittr %>%
#' 
getNormCt <- function(ct_df,value="Ct",normProbes="ALG9",probename="Probe") {
    # make subset of ct_df where gene is one of normProbes
    norm.by <- dplyr::filter(ct_df, 
                             !!dplyr::sym(probename) %in% normProbes) %>%
        .[[value]] %>%
        stats::median(na.rm=TRUE)
    
    # assign median of value to ct_df$norm.by
    # note this is the same value for every row, a waste of space technically
    ct_df %>%
        dplyr::mutate(norm.by = norm.by) %>%
        return()
}

#' Normalize cycle count (log2-fold) data within Sample
#'
#' @param ct_df a data frame containing columns "Sample", value (default Ct) and
#'   probe (default Probe). Crucially, Sample name should be the same for
#'   different technical replicates measuring identical reactions in different
#'   wells of the plate, but differ for different biological replicates.
#' @param value the column name of the value that will be normalized
#' @param normProbes names of PCR probes (or primer sets) to normalize by, i.e.
#'   reference genes
#' @param probename the column name for probe sets
#'   
#' @return data frame like ct_df with three additional columns:
#' 
#' \tabular{ll}{
#'   norm.by       \tab the median value of the reference probes  \cr
#'   Value.norm    \tab the normalized value, \eqn{\Delta Ct} \cr
#'   Value.normexp \tab the normalized ratio, \eqn{2^(-\Delta Ct)}
#'   }
#' 
#' @export
#' @importFrom magrittr %>%
#' 
normalizeqPCR <- function(ct_df,value="Ct",normProbes="ALG9",probename="Probe") {
    ct_df %>%
        dplyr::group_by(Sample) %>%
        dplyr::do(getNormCt(.,value,normProbes,probename)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(.Value = !!dplyr::sym(value), # a tidyeval trick
               Value.norm = .Value - norm.by, 
               Value.normexp =2^-Value.norm ) %>%
        dplyr::select(-.Value) %>%
        return()
}
