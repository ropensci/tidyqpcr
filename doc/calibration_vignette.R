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
# Names of target genes
Names <- c("ECM38","FET5","GPT2","ILV5","NRD1","RDL1","TFS1") 
# ORF ids of target genes
Targets  <- c("YLR299W","YFL041W","YKR067W","YLR355C","YNL251C","YOR285W","YLR178C")
# Repeats of gene names to account for resting multiple probesets
Namesrep <- c(rep(Names[1:2],each=3),rep(Names[3:7],each=2))
# id numbers of multiple probesets (reflecting IDs as ordered)
Probes <- paste(Namesrep,
                c(1,2,3,1,3,4,1,4,1,4,1,2,4,5,1,5),sep="_")


rowkey <- data.frame(WellR=LETTERS[1:16],
                     Names=Namesrep,
                     # Targets=rep(Targets,each=2),
                     Probe=factor(Probes, levels=Probes)
)


plate1plan <- 
    label_plate_rowcol(
        create_blank_plate(),
        rowkey,
        create_colkey_4dilutions_mRTNT_in24()) %>%
    mutate(Sample=paste(BioRep,DilutionNice,sep="_"))



## ----display_plates,fig.height=8,fig.width=12,dependson="label_plates"----
display_plate(plate1plan) 


## ----load_plates,dependson="label_plates",results="show"-----------------

# read my plates

plates <- read_tsv("../inst/extdata/Edward_qPCR_Nrd1_calibration_2019-02-02_Ct.txt",
                       skip=1) %>%
    mutate(Well=Pos,Ct=Cp) %>%
    right_join(plate1plan)


## ----show_plates,dependson="load_plates",results="show"------------------
plates

summary(plates)

## ----plot_unnormalized,dependson="load_plates",fig.height=6,fig.width=9----

ggplot(data=plates) +
    geom_point(aes(x=Probe,y=Ct,colour=DilutionNice,shape=Type),
               position=position_jitter(width = 0.2,height=0)) +
    labs(y="Cycle count to threshold",
         title="All reps, unnormalized") +
    facet_wrap(~BioRep) +
    panel_border() +
    theme(axis.text.x=element_text(angle=90,vjust=0.5))

## ----plot_dilutions,dependson="load_plates",fig.height=11,fig.width=6----
ggplot(data=filter(plates,Type=="+RT"),aes(x=Dilution,y=Ct)) +
    geom_point() +
    stat_smooth(formula=y~x,method="lm",se=FALSE,
                aes(colour="fit",linetype="fit")) + 
    stat_smooth(formula=y~1+offset(-x*log(10)/log(2)),method="lm",se=FALSE,
                aes(colour="theory",linetype="theory")) + 
    scale_x_log10(breaks=10^-(0:3)) +
    scale_y_continuous(breaks=seq(0,30,2)) + 
    labs(y="Cycle count to threshold",
         title="All reps, unnormalized",
         colour="Dilution") +
    facet_grid(Probe~BioRep,scales="free_y") + 
    theme(axis.text.x=element_text(angle=90,vjust=0.5))

## ----plot_dilutions_nice,dependson="load_plates",fig.height=6,fig.width=4----
Probesnice <- c("ECM38_3","FET5_1","GPT2_4","ILV5_4","NRD1_1","RDL1_4","TFS1_1")

ggplot(data=filter(plates,Type=="+RT",Probe %in% Probesnice),
       aes(x=Dilution,y=Ct)) +
    geom_point() +
    stat_smooth(formula=y~x,method="lm",se=FALSE,
                aes(colour="fit",linetype="fit")) + 
    stat_smooth(formula=y~1+offset(-x*log(10)/log(2)),method="lm",se=FALSE,
                aes(colour="theory",linetype="theory")) + 
    scale_x_log10(breaks=10^-(0:3)) +
    scale_y_continuous(breaks=seq(0,30,2)) + 
    labs(y="Cycle count to threshold",
         title="All reps, unnormalized",
         colour="Dilution") +
    facet_grid(Probe~BioRep,scales="free_y") + 
    theme(axis.text.x=element_text(angle=90,vjust=0.5))

## ----load_amp,dependson="label_plates",results="show"--------------------

plate1curve <- read_tsv("../inst/extdata/Edward_qPCR_Nrd1_calibration_2019-02-02.txt",
                       skip=2,
                        col_names=c("Well","SID","Program","Segment",
                                    "Cycle","Time","Temperature","Fluor") 
                      ) %>%
    debaseline() %>%
    left_join(plate1plan) 

# amplification curve is program 2
platesamp  <- plate1curve %>% 
  filter(Program == 2)

# melt curve is program 3 or 4, depending on cycle
platesmelt <- plate1curve %>% 
  filter(Program > 2) %>% 
  demelt() %>% 
  filter(Temperature >= 61)


## ----testplottraj,dependson="load_amp",fig.width=4,fig.height=3----------
ggplot(data=platesamp %>% filter(Well=="A1"),
       aes(x=Cycle,y=Signal)) + 
    geom_line() + 
    scale_y_log10()

## ----print_techreps,results="show",echo=TRUE,cache=FALSE,eval="FALSE"----
plate1plan %>% 
           filter(TechRep=="1",Probe==Probes[1],DilutionNice=="1x")

plate1plan %>% 
           filter(TechRep=="2",Probe==Probes[1],DilutionNice=="1x")

## ----plotamp_all,dependson="load_amp",fig.height=11,fig.width=7----------
ggplot(data=platesamp %>% 
           filter(TechRep=="1"),
       aes(x=Cycle,y=Signal,colour=factor(Dilution),linetype=Type)) + 
    facet_grid(Probe~BioRep,scales="free_y") + 
    scale_linetype_manual(values=c("+RT"="solid","-RT"="dashed","NT"="dotted")) + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,10),minor_breaks=seq(60,100,5)) + 
    labs(title="All Amp Curves, TechRep A") 

ggplot(data=platesamp %>% 
           filter(TechRep=="2"),
       aes(x=Cycle,y=Signal,colour=factor(Dilution),linetype=Type)) + 
    facet_grid(Probe~BioRep,scales="free_y") + 
    scale_linetype_manual(values=c("+RT"="solid","-RT"="dashed","NT"="dotted")) + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,10),minor_breaks=seq(60,100,5)) + 
    labs(title="All Amp Curves, TechRep B") 


## ----plotmelt_A1,dependson="load_amp",fig.width=4,fig.height=1.5,eval=FALSE----
#  ggplot(data=platesmelt %>%
#             filter(Well=="A1"),
#         aes(x=Temperature,y=dRdT)) +
#      facet_wrap(~Probe) +
#      geom_line()
#  

## ----plotmelt_all,dependson="load_amp",fig.height=11,fig.width=7---------
ggplot(data=platesmelt %>% 
           filter(TechRep=="1"),
       aes(x=Temperature,y=dRdT,colour=factor(Dilution),linetype=Type)) + 
    facet_grid(Probe~BioRep,scales="free_y") + 
    scale_linetype_manual(values=c("+RT"="solid","-RT"="dashed","NT"="dotted")) + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,10),minor_breaks=seq(60,100,5)) + 
    labs(title="All Melt Curves, TechRep A") 

ggplot(data=platesmelt %>% 
           filter(TechRep=="2"),
       aes(x=Temperature,y=dRdT,colour=factor(Dilution),linetype=Type)) + 
    facet_grid(Probe~BioRep,scales="free_y") + 
    scale_linetype_manual(values=c("+RT"="solid","-RT"="dashed","NT"="dotted")) + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,10),minor_breaks=seq(60,100,5)) + 
    labs(title="All Melt Curves, TechRep B") 


## ----plotmelt_SS_zoomed,dependson="load_amp",fig.height=11,fig.width=7----
ggplot(data=platesmelt %>% 
           filter(TechRep=="1",Type=="+RT"),
       aes(x=Temperature,y=dRdT,colour=factor(Dilution))) + 
    facet_grid(Probe~BioRep,scales="free_y") + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,5),minor_breaks=seq(60,100,1),
                       limits=c(73,87)) + 
    labs(title="Melt curves, zoomed, Superscript, TechRep A") +
    theme(panel.grid.major.x=element_line(colour="grey50",size=0.4),
          panel.grid.minor.x=element_line(colour="grey70",size=0.1))

ggplot(data=platesmelt %>% 
           filter(TechRep=="2",Type=="+RT"),
       aes(x=Temperature,y=dRdT,colour=factor(Dilution))) + 
    facet_grid(Probe~BioRep,scales="free_y") + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,5),minor_breaks=seq(60,100,1),
                       limits=c(73,87)) + 
    labs(title="Melt curves, zoomed, Superscript, TechRep B") +
    theme(panel.grid.major.x=element_line(colour="grey50",size=0.4),
          panel.grid.minor.x=element_line(colour="grey70",size=0.1))


## ----plotmelt_SS_zoomed_nice,dependson="load_amp",fig.height=6,fig.width=4----
Probesnice <- c("ECM38_3","FET5_1","GPT2_4","ILV5_4","NRD1_1","RDL1_4","TFS1_1")
ggplot(data=platesmelt %>% 
           filter(TechRep=="1",Type=="+RT",DilutionNice=="1x",
                  Probe %in% Probesnice),
       aes(x=Temperature,y=dRdT,colour=BioRep)) + 
    facet_grid(Probe~.,scales="free_y") + 
    geom_line() +
    scale_x_continuous(breaks=seq(60,100,5),minor_breaks=seq(60,100,1),
                       limits=c(73,87)) + 
    labs(title="Nice probes, TechRep A") +
    theme(panel.grid.major.x=element_line(colour="grey50",size=0.4),
          panel.grid.minor.x=element_line(colour="grey70",size=0.1))

