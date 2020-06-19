
#' Create a blank plate template as a tibble
#' 
#' For more help, examples and explanations, see the plate setup vignette:
#' \code{vignette("platesetup_vignette", package = "tidyqpcr")}
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
#' create_blank_plate_96well()
#' 
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr %>%
#' @importFrom forcats as_factor
#' 
create_blank_plate <- function(WellR=LETTERS[1:16],WellC=1:24) {
    plate <- tidyr::crossing(WellR=as_factor(WellR),
                             WellC=as_factor(WellC)) %>%
        as_tibble() %>%
        tidyr::unite(Well,WellR,WellC,sep="",remove=FALSE)
    return(plate)
}

#' @describeIn create_blank_plate create blank 96-well plate
#' @export
#' 
create_blank_plate_96well <- function() {
    return( create_blank_plate(WellR=LETTERS[1:8],WellC=1:12) )
}

#' @describeIn create_blank_plate Row names for 1536-well plates on Lightcycler 1536 Aa,Ab,Ac,Ad,Ba,...,Hd.
#' @export
#' 
make_row_names_LC_1536 <- function() {
    return( paste0(rep(LETTERS[1:8],each=4),letters[1:4]) )
}

#' @describeIn create_blank_plate Row names for 1536-well plates on Labcyte Echo A,B,...,Z,AA,AB,...,AF.
#' @export
#' 
make_row_names_Echo_1536 <- function() {
    c(LETTERS[1:26],paste0("A",LETTERS[1:6]))
}

#' @describeIn create_blank_plate create blank 1536-well plate
#' @export
#' 
create_blank_plate_1536well <- function(
    WellR=make_row_names_LC_1536(),
    WellC=1:48) {
    return( create_blank_plate(WellR,WellC) )
}

#' Create a 6-value, 24-column key for plates
#' 
#' Create a 24-column key with 6 values repeated over 24 plate columns. 
#' Each of the 6 values is repeated over 3x +RT Techreps and 1x -RT.
#' 
#' @param ... Vectors of length 6 describing well contents, e.g. sample or probe. 
#' @return tibble (data frame) with 24 rows, and columns WellC, Type, TechRep, and supplied values. 
#' @examples
#' create_colkey_6in24(SampleID=LETTERS[1:6])
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr %>%
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
#' create_rowkey_4in16(SampleID=c("sheep","goat","cow","chicken"))
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr %>%
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
#' create_rowkey_8in16_plain(SampleID=c("me","you","them","him",
#'                                    "her","dog","cat","monkey"))
#' @family plate creation functions
#' 
#' @export
#' @importFrom tibble tibble
#' @importFrom tidyr %>%
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
#'   include other information in additional variables.
#' @param rowkey tibble (data frame) describing plate rows, with variables WellR
#'   and others.
#' @param colkey tibble (data frame) describing plate columns, with variables
#'   WellC and others.
#' @param coercefactors if TRUE, coerce WellR in rowkey and WellC in colkey to factors
#' @return tibble (data frame) with variables WellR, WellC, Well, and others. 
#'   
#'   This tibble contains all combinations of WellR and WellC found in the input
#'   plate, and all information supplied in rowkey and colkey distributed across
#'   every well of the plate. Return plate is ordered by row WellR then column
#'   WellC.
#'   
#'   Note this ordering may cause a problem if WellC is supplied as a character
#'   (1,10,11,...), instead of a factor or integer (1,2,3,...). For this reason,
#'   the function by default converts WellR in `rowkey`, and WellC in `colkey`,
#'   to factors, taking factor levels from `plate`, and warns the user.
#'   
#'   Other tidyqpcr functions require plate plans to contain variables SampleID,
#'   TargetID, and Type, so `label_plate_rowcol` will warn if any of these are
#'   missing. This is a warning, not an error, because these variables can be
#'   added by users later.
#'   
#' @examples
#' label_plate_rowcol(plate = create_blank_plate()) # returns blank plate
#' @family plate creation functions
#' 
#' @export
#' 
label_plate_rowcol <- function(plate,rowkey=NULL,colkey=NULL,coercefactors=TRUE) {
    if (!is.null(colkey)) {
        assertthat::assert_that(has_name(colkey,"WellC"))
        # Note: should this if clause be a freestanding function?
        # coerce_column_to_factor(df, col, warn=TRUE) ?
        if( !is.factor(colkey$WellC) & coercefactors ) {
            warning("coercing WellC to a factor with levels from plate$WellC")
            colkey <- dplyr::mutate(colkey,
                                    WellC=factor(WellC,
                                                 levels=levels(plate$WellC))
                                    )
        }
        plate <- dplyr::left_join(plate,colkey,by="WellC")
    }
    if (!is.null(rowkey)) {
        assertthat::assert_that(has_name(rowkey,"WellR"))
        if( !is.factor(rowkey$WellR) & coercefactors ) {
            warning("coercing WellR to a factor with levels from plate$WellR")
            rowkey <- dplyr::mutate(rowkey,
                                    WellR=factor(WellR,
                                                 levels=levels(plate$WellR))
                                    )
        }
        plate <- dplyr::left_join(plate,rowkey,by="WellR")
    }
    # check that plate contains SampleID, TargetID, Type, warn if not
    if( ! "SampleID" %in% names(plate) ) {
        warning("plate does not contain variable SampleID")
    }
    if( ! "TargetID" %in% names(plate) ) {
        warning("plate does not have variable TargetID")
    }
    if( ! "Type" %in% names(plate) ) {
        warning("plate does not have variable Type")
    }
    return( dplyr::arrange( plate, WellR, WellC ) )
}


#' Display plate plan with SampleID, TargetID, Type per Well
#'
#' @param plate tibble with variables WellC, WellR, SampleID, TargetID, Type. 
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
        ggplot2::geom_tile(aes(fill=TargetID),alpha=0.3) +
        ggplot2::geom_text(aes(label=paste(TargetID,SampleID,Type,sep="\n")),
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
#' @param normby.function Function to use to calculate the value to 
#' normalise by on log2/Cq scale. 
#' Default value is median, alternatively could use mean.
#' 
#' @export
#' @importFrom tidyr %>%
#' 
getNormCq <- function(cq_df,value="Cq",normTargetIDs="ALG9",probename="TargetID",
                      normby.function=median) {
    # make subset of cq_df where gene is one of normTargetIDs
    norm.by <- dplyr::filter(cq_df, 
                             !!dplyr::sym(probename) %in% normTargetIDs) %>%
        .[[value]] %>%
        normby.function(na.rm=TRUE)
    
    # assign median of value to cq_df$norm.by
    # note this is the same value for every row, a waste of space technically
    cq_df %>%
        dplyr::mutate(norm.by = norm.by) %>%
        return()
}

#' Normalize cycle count (log2-fold) data within SampleID
#'
#' @param cq_df a data frame containing columns "SampleID", value (default Cq) and
#'   probename (default TargetID). Crucially, SampleID should be the same for
#'   different technical replicates measuring identical reactions in different
#'   wells of the plate, but differ for different biological and experimental 
#'   replicates.
#' @param value the column name of the value that will be normalized
#' @param normTargetIDs names of PCR probes (or primer sets) to normalize by, i.e.
#'   reference genes
#' @param probename the column name for probe sets
#'   
#' @return data frame like cq_df with three additional columns:
#' 
#' \tabular{ll}{
#'   norm.by       \tab the median value of the reference probes  \cr
#'   Value.norm    \tab the normalized value, \eqn{\Delta Cq} \cr
#'   Value.normexp \tab the normalized ratio, \eqn{2^(-\Delta Cq)}
#'   }
#' 
#' @export
#' @importFrom tidyr %>%
#' 
normalizeqPCR <- function(cq_df,value="Cq",normTargetIDs="ALG9",probename="TargetID") {
    cq_df %>%
        dplyr::group_by(SampleID) %>%
        dplyr::do(getNormCq(.,value,normTargetIDs,probename)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(.Value = !!dplyr::sym(value), # a tidyeval trick
               Value.norm = .Value - norm.by, 
               Value.normexp =2^-Value.norm ) %>%
        dplyr::select(-.Value) %>%
        return()
}

#' @describeIn normalizeqPCR Normalise cycle count (log2-fold) data within SampleID.
#' Synonym for normalizeqPCR.
#' 
#' @export
#' @importFrom tidyr %>%
#' 
normaliseqPCR <- function(cq_df,value="Cq",normTargetIDs="ALG9",probename="TargetID") {
    normalizeqPCR(cq_df=cq_df,value=value,normTargetIDs=normTargetIDs,probename=probename)
} 
