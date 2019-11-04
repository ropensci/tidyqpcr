## ----setup,warning=FALSE,message=FALSE,echo=FALSE------------------------

## knitr options for report generation
knitr::opts_chunk$set(warning=FALSE,message=FALSE,echo=TRUE,cache=FALSE,
                      results="show")
library(tidyverse)
library(cowplot)
library(tidyqpcr)

 # set default theme for graphics
theme_set(theme_cowplot(font_size=11) %+replace% 
              theme(panel.border=element_rect(colour = "grey50",
                                            linetype = "solid",size=0.5),
                    strip.background = element_blank()))

## ----label_plates,dependson="plate_functions"----------------------------

# list Targets /Probes
Probes <- c("HOR7", "HSP12", "HSP26", "HSP78", 
            "HSP104", "RTC3", "SSA4", "PGK1", 
            "ALG9", " HHT2", "HTB2", "RPS3", 
            "RPS13", "RPS15", "RPS30A", "RPL39")

rowkey <- tibble(WellR=LETTERS[1:16],
                 Probe=factor(Probes,levels=Probes) 
                 ) 

# Set up experimental samples
HSlevels <- c("-","+")
HSvalues <- factor(rep(HSlevels,each=3),levels=HSlevels)
Druglevels <- c("C","P","T")
Drugvalues <- factor(rep(Druglevels,times=2),levels=Druglevels)
Conditionlevels <- paste0(Druglevels,rep(HSlevels,each=3))
Conditionvalues <- factor(Conditionlevels,levels=Conditionlevels)

colkey <- create_colkey_6in24(HS=HSvalues,
                              Drug=Drugvalues,
                              Condition=Conditionvalues)

plateplan <-     
    label_plate_rowcol(create_blank_plate(WellR = LETTERS[1:16],WellC=1:24),
                       rowkey,colkey)


## ----display_plates,fig.height=9,fig.width=10,dependson="label_plates"----
display_plate( plateplan %>%
    mutate(Sample=Condition) )

## ----load_plates,dependson="label_plates",results="show"-----------------
# set working directory, if needed
# setwd("~/Dropbox (Drummond Lab)/EdinburghLab/RT-qPCR/2018-06/")
# read my plates

plate1 <- read_tsv("../inst/extdata/Edward_qPCR_TxnInhibitors_HS_2018-06-15_plate1_Ct.txt",skip=1) %>%
    mutate(Well=Pos,Ct=Cp) %>%
    left_join(plateplan) %>%
    mutate(BiolRep="1",Plate="1")


plate2 <- read_tsv("../inst/extdata/Edward_qPCR_TxnInhibitors_HS_2018-06-15_plate2_Ct.txt",skip=1) %>%
    mutate(Well=Pos,Ct=Cp) %>%
    left_join(plateplan) %>%
    mutate(BiolRep="2",Plate="2")

plates <- bind_rows(plate1,plate2) %>%
    mutate( Sample=paste0(Condition,BiolRep) )

summary(plates)


## ----plot_unnormalized,dependson="load_plates",fig.height=6,fig.width=8----
ggplot(data=plates) +
    geom_point(aes(x=Probe,y=Ct,shape=Condition,colour=Condition),
               position=position_jitter(width = 0.2,height=0)) +
    labs(y="Cycle count to threshold",
         title="All reps, unnormalized") +
    scale_shape_manual(values=c(15:18,5:6)) + 
    facet_grid(BiolRep~Type) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5),
          panel.border=element_rect(fill = NA,linetype=1,
                                    colour = "grey50",size=0.5))


## ----normalize_counts,dependson="load_plates"----------------------------
platesnorm <- plates %>% 
    filter(Type=="+RT") %>%
    normalizeqPCR(normProbes = "PGK1")

platesmed <- platesnorm %>%
    group_by(Sample,Condition,BiolRep,HS,Drug,Sample,Probe) %>%
    summarize(Ct=median(Value.norm,na.rm=TRUE),
              Abund=median(Value.normexp,na.rm=TRUE))

filter(platesmed,Probe=="HSP26")


## ----plot_normalized,dependson="normalize_counts",fig.height=6,fig.width=6----
ggplot(data=platesnorm) +
    geom_point(aes(x=Probe,y=Value.norm,shape=Condition,colour=Condition),
               position=position_jitter(width = 0.2,height=0)) +
    labs(y="Ct relative to PGK1") +
    scale_shape_manual(values=c(15:18,5:6)) + 
    facet_grid(BiolRep~.) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5),
          panel.border=element_rect(fill = NA,linetype=1,
                                    colour = "grey50",size=0.5))


## ----plot_normalizedsummarized1,dependson="normalize_counts",fig.height=3,fig.width=4----
ggplot(data=platesmed) +
    geom_point(aes(x=Probe,y=Abund,shape=BiolRep,colour=Condition),
               position=position_jitter(width = 0.2,height=0)) +
    scale_shape_manual(values=c(15:18,5:6)) + 
    scale_y_log10nice("mRNA relative detection") + 
    theme(axis.text.x=element_text(angle=90,vjust=0.5),
          panel.border=element_rect(fill = NA,linetype=1,
                                    colour = "grey50",size=0.5))


## ----plot_normalizedsummarizedbyHS,dependson="normalize_counts",fig.height=3,fig.width=9----
ggplot(data=platesmed) +
    geom_point(aes(x=Probe,y=Abund,shape=BiolRep,colour=HS),
               position=position_jitter(width = 0.2,height=0)) +
    facet_wrap(~Drug,ncol=3) +
    scale_colour_manual(values=c("-"="grey50","+"="red")) +
    scale_y_log10nice("mRNA relative detection") + 
    theme(axis.text.x=element_text(angle=90,vjust=0.5),
          panel.border=element_rect(fill = NA,linetype=1,
                                    colour = "grey50",size=0.5))


## ----plot_normalizedsummarizedbyDrug,dependson="normalize_counts",fig.height=3,fig.width=9----
ggplot(data=platesmed) +
    geom_point(aes(x=Probe,y=Abund,shape=BiolRep,colour=Drug),
               position=position_jitter(width = 0.2,height=0)) +
    facet_wrap(~HS,ncol=3) +
    scale_y_log10nice("mRNA relative detection") + 
    theme(axis.text.x=element_text(angle=90,vjust=0.5),
          panel.border=element_rect(fill = NA,linetype=1,
                                    colour = "grey50",size=0.5))


## ----load_amp,dependson="label_plates",results="show"--------------------

plate1curve <- read_tsv("../inst/extdata/Edward_qPCR_TxnInhibitors_HS_2018-06-15_plate1.txt",
                       skip=2,
                        col_names=c("Well","SID","Program","Segment",
                                    "Cycle","Time","Temperature","Fluor") 
                      ) %>%
    debaseline() %>%
    left_join(plateplan) %>%
    mutate(BiolRep=1,Plate=1)

plate2curve <- read_tsv("../inst/extdata/Edward_qPCR_TxnInhibitors_HS_2018-06-15_plate2.txt",
                       skip=2,
                        col_names=c("Well","SID","Program","Segment",
                                    "Cycle","Time","Temperature","Fluor") 
                      ) %>%
    debaseline() %>%
    left_join(plateplan)  %>%
    mutate(BiolRep=2,Plate=2)

# platesamp  <- plate1curve %>% filter(Program == 2)

# plate1curve
# plate2curve

platesamp <- bind_rows(plate1curve,plate2curve) %>% 
    filter(Program == 2)

platesmelt <- bind_rows(plate1curve,plate2curve) %>% 
    filter(Program != 2) %>% 
    demelt() %>% 
    filter(Temperature >= 61)  


## ----plotmelt_Rep1,dependson="load_amp",fig.width=12,fig.height=6--------
ggplot(data=platesmelt %>% 
           filter(TechRep==1,BiolRep==1),
       aes(x=Temperature,y=dRdT,linetype=Type)) + 
    facet_grid(Condition~Probe) + 
    geom_line() +
    scale_linetype_manual(values=c("-RT"="dashed","+RT"="solid")) +
    scale_x_continuous(breaks=seq(60,100,10),minor_breaks=seq(60,100,5)) + 
    labs(title="Melt curves, BiolRep 1, Techrep 1") + panel_border()

## ----plotmelt_Rep2,dependson="load_amp",fig.width=12,fig.height=6--------
ggplot(data=platesmelt %>% 
           filter(TechRep==1,BiolRep==2),
       aes(x=Temperature,y=dRdT,linetype=Type)) + 
    facet_grid(Condition~Probe) + 
    geom_line() +
    scale_linetype_manual(values=c("-RT"="dashed","+RT"="solid")) +
    scale_x_continuous(breaks=seq(60,100,10),minor_breaks=seq(60,100,5)) + 
    labs(title="Melt curves, BiolRep 2, TechRep 1") + panel_border()

## ----plotamp_Rep1,dependson="load_amp",fig.width=12,fig.height=6---------
ggplot(data=platesamp %>% 
           filter(TechRep==1,BiolRep==1),
       aes(x=Cycle,y=Fluor,linetype=Type)) + 
    facet_grid(Condition~Probe) + 
    geom_line() +
    scale_linetype_manual(values=c("-RT"="dashed","+RT"="solid")) +
    expand_limits(y=0) + 
    labs(title="Amp. curves, BiolRep 1, Techrep 1") + panel_border()

## ----plotamp_Rep2,dependson="load_amp",fig.width=12,fig.height=6---------
ggplot(data=platesamp %>% 
           filter(TechRep==1,BiolRep==2),
       aes(x=Cycle,y=Fluor,linetype=Type)) + 
    facet_grid(Condition~Probe) + 
    geom_line() +
    scale_linetype_manual(values=c("-RT"="dashed","+RT"="solid")) +
    expand_limits(y=0) + 
    labs(title="Amp. curves, BiolRep 2, Techrep 1") + panel_border()

