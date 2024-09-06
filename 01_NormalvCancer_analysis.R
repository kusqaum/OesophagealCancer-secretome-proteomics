library(janitor)
library(plyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(cowplot)
library(patchwork)
library(grid)
library(gridExtra)
library(broom)
library(clusterProfiler)

Master_pep <- clean_names(read.delim("raw/Master_pep_KA.txt"))
colnames(Master_pep)

# Subset data into Normal data-set and cancer data-set
normal <- dplyr::filter(Master_pep, type == "Oesophagus")
length(unique(normal$leading_proteins))
#767

cancer <- dplyr::filter(Master_pep, type!= "Oesophagus")
length(unique(cancer$leading_proteins))
#1700


# Focusing on normal cells
# create a Unique identifier (extra column) of Gene Name + UniProt ID:
normal$Unique_ID <- paste(normal$gene_names, normal$leading_proteins, sep = " - ")
cancer$Unique_ID <- paste(cancer$gene_names, cancer$leading_proteins, sep = " - ")
# create unique peptide count per protein plus number of time-points that we have RIA, t data:
pcountN <- normal %>% group_by(Unique_ID) %>% summarise(pep_count=n_distinct(peptide_sequence))
tcountN <- normal %>% group_by(Unique_ID) %>% summarise(time_count=n_distinct(time))

# merge peptide count and time count data with main data:
normal_1 <- merge(normal,pcountN,by = "Unique_ID")
normal_1 <- merge(normal_1,tcountN,by = "Unique_ID")
length(unique(normal_1$Unique_ID))
#768

# filter to keep the IDs (Proteins) where we have RIA,t data at all 5 time-points
# i.e.the "best" IDs. This is quite stringent filtering of the data.
filtN <- dplyr::filter(normal_1, time_count>4)
length(unique(filtN$Unique_ID))
# 258

# write a nice plotting function for RIA,t data:
plot_RIA_n=function(d){
  ggplot(d,aes(x=time,y=ria)) +
    geom_point(pch=4,size=0.75) + 
    ylab("RIA") + xlab("Time (h)") + ylim(0,1) + xlim(0,24) +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "tomato", size = 0.45) +
    facet_wrap( ~ Unique_ID) +
    theme(legend.position="none",
          panel.border = element_rect(colour = "black",fill = NA),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black"))  
}
# run this function on our subsetted data set of high qaulity hits (IDs).
p1 <- dlply(filtN,"Unique_ID",plot_RIA_n)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q1 <- marrangeGrob(p1,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/RIA_Plots_For_Normal_Oesophageal_myofibroblasts.pdf",q1,width=6,height=9)
# clear memory
gc()

# to write out .png but only last page:
png("output/test.png")
q1
dev.off()

### do the same for the cancer cells with slightly different appearance of plots to help
# make them distinct:

# create unique peptide count per protein plus number of time-points that we have RIA, t data:
pcountC <- cancer %>% group_by(Unique_ID) %>% summarise(pep_count=n_distinct(peptide_sequence))
tcountC <- cancer %>% group_by(Unique_ID) %>% summarise(time_count=n_distinct(time))
# merge peptide count and time count data with main data:
cancer_1 <- merge(cancer,pcountC,by = "Unique_ID")
cancer_1 <- merge(cancer_1,tcountC,by = "Unique_ID")

# filter to keep the IDs (Proteins) where we have RIA,t data at all 5 time-points
# i.e.the "best" IDs. This is quite stringent filtering of the data.
filtC <- dplyr::filter(cancer_1, time_count>4)
length(unique(filtC$Unique_ID))
# 895

plot_RIA_c=function(d){
  ggplot(d,aes(x=time,y=ria)) +
    geom_point(pch=4,size=0.75) + 
    ylab("RIA") + xlab("Time (h)") + ylim(0,1) + xlim(0,24) +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "dodgerblue", size = 0.45) +
    facet_wrap( ~ Unique_ID) +
    theme(legend.position="none",
          panel.border = element_rect(colour = "black",fill = NA),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black"))  
}
# run this function on our subsetted data set of high qaulity hits (IDs).
p2 <- dlply(filtC,"Unique_ID",plot_RIA_c)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q2 <- marrangeGrob(p2,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/RIA_Plots_For_Oesophageal_cancer_cells.pdf",q2,width=6,height=9)
# clear memory
gc()

# this reports errors indicative of failed fitting (nls). So, the plots
# have no curve/fit-line. We looked at thes multiplots and found these
# so we need to filter them out as they are poor qaulity hits, in terms of
# ability to plot RIA with time
# these the hits to filter:
filtC <- dplyr::filter(filtC, Unique_ID!= "CSTF2 - P33240")
filtC <- dplyr::filter(filtC, Unique_ID!= "DYNLRB1 - Q9NP97")
filtC <- dplyr::filter(filtC, Unique_ID!= "MTDH - Q86UE4")
filtC <- dplyr::filter(filtC, Unique_ID!= "PSMD11 - O00231")
filtC <- dplyr::filter(filtC, Unique_ID!= "PSMD14 - O00487")
filtC <- dplyr::filter(filtC, Unique_ID!= "RPS14 - P62263")
filtC <- dplyr::filter(filtC, Unique_ID!= "SEC23A - Q15436")
filtC <- dplyr::filter(filtC, Unique_ID!= "SETD3 - Q86TU7")
filtC <- dplyr::filter(filtC, Unique_ID!= "ZC3H6 - P61129")
filtC <- dplyr::filter(filtC, Unique_ID!= "CHMP3 - Q9Y3E7")
filtC <- dplyr::filter(filtC, Unique_ID!= "DYNC1LI2 - O43237")


p2 <- dlply(filtC,"Unique_ID",plot_RIA_c)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q2 <- marrangeGrob(p2,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/RIA_Plots_For_Oesophageal_cancer_cells.pdf",q2,width=6,height=9)
# clear memory
gc()



# we now have our multi-page plots for normal Supplementary data. These have all data points
# on them. Let's see if they look better if we first mean the RIA,t data?
filtN_mean_RIA <- filtN %>%
  group_by(Unique_ID, time) %>%
  summarise(mean_RIA=mean(ria,na.rm=TRUE))
filtN_mean_RIA <- as.data.frame(filtN_mean_RIA)

plot_RIA_n_mean=function(d){
  ggplot(d,aes(x=time,y=mean_RIA)) +
    geom_point(pch=1,size=0.75) + 
    ylab("RIA") + xlab("Time (h)") + ylim(0,1) + xlim(0,24) +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "blue", size = 0.45, linetype="longdash") +
    facet_wrap( ~ Unique_ID) +
    theme(legend.position="none",
          panel.border = element_rect(colour = "black",fill = NA),
          panel.background = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black"))  
}
# run this function on our subsetted data set of high qaulity hits (IDs).
p3 <- dlply(filtN_mean_RIA,"Unique_ID",plot_RIA_n_mean)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q3 <- marrangeGrob(p3,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/Meaned_RIA_Plots_Normals_myofibroblasts.pdf",q3,width=6,height=9)
# clear memory
gc()

#Now we can do this for cancer cells as well
# we now have our multi-page plots for normal Supplementary data. These have all data points
# on them. Let's see if they look better if we first mean the RIA,t data?
filtC_mean_RIA <- filtC %>%
  group_by(Unique_ID, time) %>%
  summarise(mean_RIA=mean(ria,na.rm=TRUE))
filtC_mean_RIA <- as.data.frame(filtC_mean_RIA)

plot_RIA_c_mean=function(d){
  ggplot(d,aes(x=time,y=mean_RIA)) +
    geom_point(pch=1,size=0.75) + 
    ylab("RIA") + xlab("Time (h)") + ylim(0,1) + xlim(0,24) +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "tomato", size = 0.45, linetype="longdash") +
    facet_wrap( ~ Unique_ID) +
    theme(legend.position="none",
          panel.border = element_rect(colour = "black",fill = NA),
          panel.background = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black"))  
}
# run this function on our subsetted data set of high quality hits (IDs).
p3c <- dlply(filtC_mean_RIA,"Unique_ID",plot_RIA_c_mean)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q3c <- marrangeGrob(p3c,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/Meaned_RIA_Plots_For_Oesophageal_cancer_cells.pdf",q3c,width=6,height=9)
# clear memory
gc()


# we will ignore plotting all abundance data points. Instead, just focus on the total
# abundance.
# normals first:
#filtN_mean_Ab=filtN %>%
#  group_by(Unique_ID, time) %>%
#  summarise(mean_ab=mean(intensity,na.rm=TRUE))
filtN_total_Ab <- filtN %>%
  group_by(Unique_ID, time) %>%
  summarise(total_ab=sum(intensity,na.rm=TRUE))
filtN_total_Ab <- as.data.frame(filtN_total_Ab)

plot_Ab_n=function(d){
  ggplot(d, aes(time, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "darkgreen", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
}
# run this function on our subsetted data set of high quality hits (IDs).
p4 <- dlply(filtN_total_Ab,"Unique_ID",plot_Ab_n)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q4 <- marrangeGrob(p4,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/Abundance_Plots_Normals_Meaned.pdf",q4,width=6,height=9)
# clear memory
gc()

# redo this for cancer cells.
# match mean RIA and mean Ab multiplots for each ID across normal and cancer to pick out
# proteins to discuss in report

filtC_total_Ab <- filtC %>%
  group_by(Unique_ID, time) %>%
  summarise(total_ab=sum(intensity,na.rm=TRUE))
filtC_total_Ab <- as.data.frame(filtC_total_Ab)

plot_Ab_c=function(d){
  ggplot(d, aes(time, total_ab)) +
    geom_point(pch = 19, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "violetred2", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
}

# run this function on our subsetted data set of high qaulity hits (IDs).
p4c <- dlply(filtC_total_Ab,"Unique_ID",plot_Ab_c)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
q4c <- marrangeGrob(p4c,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/Abundance_Plots_Cancers_Meaned.pdf",q4c,width=6,height=9)
# clear memory
gc()


#back to normals#

# what's the dynamic range of abundance at 24 hours (where, if a protein is secreted,
# this will be at its max)
range(filtN_total_Ab$total_ab)
#  1563900 10593890000
# re=draw abundance plots specifying range of abundance on y-axis. Create another plot
# abundance function:
plot_Ab_n_fullrange=function(d){
  ggplot(d, aes(time, total_ab)) +
    geom_point(pch = 8, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "darkgreen",
                linetype="dotdash",size = 0.35) +
    xlim(0, 24) + ylim(0, 10593890000) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
}
# run this function on our subsetted data set of high qaulity hits (IDs).
pa <- dlply(filtN_total_Ab,"Unique_ID",plot_Ab_n_fullrange)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
qa <- marrangeGrob(pa,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/Abundance_Plots_Normals_Fullrange_Meaned.pdf",qa,width=6,height=9)
# clear memory
gc()


# OK. We want to get dataframes (tables) of:
# 1. non-linear curve-fitting = the rate of synthesis and secretion of newly-made proteins (k)
# obtained from RIA of protein (y) with time (x) data.
# 2. slope of the linear regression fits of how abundance in secretome (y) changes with time (x)

lm_norms <- filtN_total_Ab %>% nest(data = -Unique_ID) %>%
  mutate(fit=map(data, ~ lm(total_ab ~ time,data=.x)),
         tidied=map(fit, broom::tidy), augmented=map(fit,broom::augment))
# need to unnest the nested fit data to get slope of linear regression line:
lm_norms_data <- unnest(lm_norms,cols = tidied)
lm_norms_data_augmented <- unnest(lm_norms,cols = augmented)
# just get the slope info (which is annoyingly called "time")
lm_norms_data_2 <- dplyr::filter(lm_norms_data,term=="time")
# extract only the columns we are interested in:
lm_norms_data_2 <- lm_norms_data_2[,c(1,5)]
# rename estimate column to slope:
colnames(lm_norms_data_2)[2] <- "slope"

# to do the fitted curves on one plot, unnest the augment data:
lm_norms_data_augment <- unnest(lm_norms,cols = augmented)
colnames(lm_norms_data_augment)[7] <- "Fitted"
fits_norms <- dplyr::select(lm_norms_data_augment,Unique_ID,time,total_ab,Fitted)
# need only accessions:
fits_norms$Unique_ID <- sub(".*- ","",fits_norms$Unique_ID)
# merge this with our Outcyte, SigP and SecP data web tools output:
Norms_web_tools <- read.delim("raw/norms_accessions.txt",header=TRUE)
# ignore SecP data:
Norms_web_tools <- Norms_web_tools[-3]
Norms_web_tools$Type <- ifelse((Norms_web_tools$SigP_YN=="TRUE" &
                               Norms_web_tools$Outcyte_UPS_YN=="FALSE" |
                               Norms_web_tools$SigP_YN=="TRUE" &
                               Norms_web_tools$Outcyte_UPS_YN=="TRUE"),"Classical",
                            ifelse(Norms_web_tools$SigP_YN=="FALSE" &
                                     Norms_web_tools$Outcyte_UPS_YN=="TRUE","UPS",
                                   "Leaked / Dead"))
# merge with the fitted line data:
Norms_merged <- merge(Norms_web_tools,fits_norms,by="Unique_ID")
# draw the plot object now, making each line "class" transparent a little bit so that if
# they overlap then you can see them on top of one another:
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "PRDX1")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "P10599")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "P23528")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "P04075")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "P14618")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "P21333")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "P60709")
Norms_merged <- dplyr::filter(Norms_merged, Unique_ID!= "Q05682")

Normsplot <- ggplot(Norms_merged,aes(x=time,y=Fitted,group=Unique_ID,colour=Type)) +
  geom_line(stat="smooth",method = "lm", fullrange = TRUE, se=FALSE, size=0.25,alpha=0.75) +
  ylab("Abundance") + ylim(0, 2.0E+9) +
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black",fill = NA),
        axis.ticks = element_line(colour = "black"))
# write a .pdf out of the plot object
pdf("output/Norms_classifications.pdf")
Normsplot
dev.off()

# get rid of the Dead / Leaked to see if it makes it look better.
Norms_merged_9 <- dplyr::filter(Norms_merged, Unique_ID!= "PRDX1")
Norms_merged_9 <- dplyr::filter(Norms_merged_9, Unique_ID!= "P10599")
Norms_merged_9 <- dplyr::filter(Norms_merged_9, Unique_ID!= "P23528")
Norms_merged_9 <- dplyr::filter(Norms_merged_9,Type!="Leaked / Dead")
Normsplot_9 <- ggplot(Norms_merged_9,aes(x=time,y=Fitted,group=Unique_ID,colour=Type)) +
  geom_line(stat="smooth",method = "lm", fullrange = TRUE, se=FALSE, size=0.25,alpha=0.75) +
  ylab("Abundance") +
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black",fill = NA),
        axis.ticks = element_line(colour = "black"))
# write a .pdf out of the plot object
pdf("output/UPS_and_Classical_only_Norms_classifications.pdf")
Normsplot_9
dev.off()


##
# get rid of the UPS to see if it makes it look better.
Norms_merged_5 <- dplyr::filter(Norms_merged, Unique_ID!= "P04075")
Norms_merged_5 <- dplyr::filter(Norms_merged_5, Unique_ID!= "P14618")
Norms_merged_5 <- dplyr::filter(Norms_merged_5, Unique_ID!= "P21333")
Norms_merged_5 <- dplyr::filter(Norms_merged_5, Unique_ID!= "P60709")
Norms_merged_5 <- dplyr::filter(Norms_merged_5, Unique_ID!= "Q05682")
Norms_merged_5 <- dplyr::filter(Norms_merged_5,Type!="UPS")

Normsplot_5 <- ggplot(Norms_merged_5,aes(x=time,y=Fitted,group=Unique_ID,colour=Type)) +
  geom_line(stat="smooth",method = "lm", fullrange = TRUE, se=FALSE, size=0.25,alpha=0.75) +
  ylab("Abundance") +
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black",fill = NA),
        axis.ticks = element_line(colour = "black"))
# write a .pdf out of the plot object
pdf("output/UPS_and_Classical_only_Norms_classifications_p.pdf")
Normsplot_5
dev.off()



# to do the fitted curves on one plot, unnest the augment data:
lm_norms_data_augmented <- unnest(lm_norms,cols = augmented)
colnames(lm_norms_data_augmented)[7] <- "Fitted"
fits_norms <- dplyr::select(lm_norms_data_augmented,Unique_ID,time,total_ab,Fitted)
#need only accessions
fits_norms$Unique_ID <- sub(".*- ", "", fits_norms$Unique_ID)
#merge this with our outvyte, sigp, and secp data web tools output:
Norms_web_tools <- read.delim("raw/norms_accessions.txt", header = TRUE)
#ignore SecP daat
Norms_web_tools <- Norms_web_tools[-3]
Norms_web_tools$Type <- ifelse((Norms_web_tools$SigP_YN=="TRUE" &
                               Norms_web_tools$Outcyte_UPS_YN=="FALSE" |
                               Norms_web_tools$SigP_YN=="TRUE" &
                               Norms_web_tools$Outcyte_UPS_YN=="TRUE"),"Classical",
                            ifelse(Norms_web_tools$SigP_YN=="FALSE" &
                                     Norms_web_tools$Outcyte_UPS_YN=="TRUE","UPS",
                                   "Leaked / Dead"))
#MERGE with the fitted line data:


Norms_merged <- merge(Norms_web_tools,fits_norms,by="Unique_ID")
Norms_merged %>% count(SigP_YN)

Norms_merged <- Norms_merged[-2]
Norms_merged <- Norms_merged[-4]
N_subset <- Norms_merged[1:2]
Norms_merged <- Norms_merged[-2]

#bind rows :
Norms_merged_2 <- bind_rows(Norms_merged, N_zeroes)
#re-add "Type" column:
Norms_merged_3 <- (unique(Norms_merged_2$Type))
view(Norms_merged_3)
#create t=0 abundance=0 df:
N_zeroes <- unique(dplyr::select(Norms_merged, Unique_ID))

N_zeroes$Time=0
N_zeroes$Fitted=0
#merge again
Norms_merged_2t <- bind_rows(Norms_merged_3, N_zeroes) 
# draw the plot object now, making each line "class" transparent a little bit so that if
# they overlap thenyou can see them on top of one another:
Normsplot_2 <- ggplot(Norms_merged_2t,aes(x=time,y=Fitted,group=Unique_ID,colour=Type)) +
  geom_line(stat="smooth",method = "lm", fullrange = TRUE, se=FALSE, size=0.5,alpha=0.75) +
  ylab("Abundance") +
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black",fill = NA),
        axis.ticks = element_line(colour = "black"))
# write a .pdf out of the plot object
#pdf("Norms_classification.pdf")
Normsplot_2
dev.off()


# let's do the same for the RIA,t non-linear curve-fits now:
nls_norms <- filtN_mean_RIA %>% nest(data= -Unique_ID) %>%
  mutate(fit = map(data, ~nls(mean_RIA ~ (0 + (1 - 0) * (1 - exp(-k * time))),
                              start=list(k=0.01),data=.x)),
         tidied=map(fit,broom::tidy),
         augmented=map(fit,broom::augment))
# need to unnest the nested nls fit data to get kinetic rate of synthesis and 
# secretion (k):
nls_norms_data <- unnest(nls_norms,cols = tidied)
# extract only the columns we are interested in:
nls_norms_data_2 <- nls_norms_data[,c(1,5)]
# rename estimate column to k = rate of synthesis and secretion:
colnames(nls_norms_data_2)[2] <- "k"

# filter out proteins whose abundance does not increase with time in the secretome.
# These are any with a negative slope value from the linear regression fit of abundance (y)
# with time (x) data (= in lm_norms_data_2 output)
# create a smaller version of lm_norms_data_2 with the negative slope hits removed:
lm_norms_data_3 <- dplyr::filter(lm_norms_data_2,slope>0)

# merge this lm output dataframe with our abundance data ignoring those with negative
# slope data:
filtN_total_Ab_2 <- merge(filtN_total_Ab,lm_norms_data_3,by="Unique_ID",all=FALSE)
length(unique(filtN_total_Ab_2$Unique_ID))
# down to 135 final filtered hits.

# we can filter the RIA data to only give us RIA of these proteins too as our FINAL
# filtered datasets to work on:
nls_norms_data_3 <- merge(nls_norms_data_2,lm_norms_data_3,by="Unique_ID",all=FALSE)
length(unique(nls_norms_data_3$Unique_ID))
# correct length 135

# merge with RIA data now:
filtN_mean_RIA_2 <- merge(filtN_mean_RIA,nls_norms_data_3,by="Unique_ID",all=FALSE)
length(unique(filtN_total_Ab_2$Unique_ID))
# correct length = 135.

# FINAL Normal merged data:
FINAL_norms <- merge(filtN_total_Ab_2,filtN_mean_RIA_2,by="Unique_ID",all=FALSE)
# let's just get the slope and fit data now from this:
S_and_k_norms <- dplyr::select(FINAL_norms,c(Unique_ID,slope.x,k))
S_and_k_norms <- unique(S_and_k_norms)
# rename column names:
colnames(S_and_k_norms)=c("Unique_ID","slope_n","k_n")



#snapshot of part of our code- This is the double plot function use to create the 
# We can now run our double plot function:



#Open up a function call
double_plot <- function(d){
  # RIA-t plot first:
  #Open up ggplot object, pointing to RIA data
  ria <- ggplot(d, aes(time.y, mean_RIA)) +
    #add points to that plot
    geom_point(pch = 23, size = 0.7, fill = "#F0808085") +
    #non-linear curve fitting with all parameters of the line
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#F0808085", size = 0.35) +
    #add limits and labels to axes
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #do multiplot by splitting proteins up on the basis of Unique_ID
    facet_wrap(~ Unique_ID) +
    #theme of the plot
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  #Open up ggplot object, pointing to abundance data
  ab <- ggplot(d, aes(time.x, total_ab)) +
    #add points to that plot
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    #linear fitting with all parameters of the line 
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#6CA6CD85", size = 0.35) +
    #add limits and labels to axes
    xlim(0, 24)
  ylab("Abundance") + xlab("Time (h)") +
    #do multiplot by splitting proteins up on the basis of Unique_ID
    facet_wrap(~ Unique_ID) +
    #theme of the plot
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  #tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(ria + ab)
  #close function call
}
# do the multiplots and save:
plot_list_N <- dlply(FINAL_norms,"Unique_ID",double_plot)
N_grob <- marrangeGrob(plot_list_N,nrow=5,ncol=1)
ggsave("output/Side_by_Side_Multiplots_Normals.pdf",N_grob,width=5,height=10)

#now
#back to cancers#
# what's the dynamic range of abundance at 24 hours (where, if a protein is secreted,
# this will be at its max)
range(filtC_total_Ab$total_ab)
## 1129200 24076293000
# re=draw abundance plots specifying range of abundance on y-axis. Create another plot
# abundance function:
plot_Ab_c_fullrange=function(d){
  ggplot(d, aes(time, total_ab)) +
    geom_point(pch = 8, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "darkgreen",
                linetype="dotdash",size = 0.35) +
    xlim(0, 24) + ylim(0, 24076293000) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
}
# run this function on our subsetted data set of high quality hits (IDs).
pc <- dlply(filtC_total_Ab,"Unique_ID",plot_Ab_c_fullrange)
# create a multipanel page to fill with plots for each ID over multiple pages, all
# scaled the same
qc <- marrangeGrob(pc,nrow=4,ncol=3)
# save these as .pdf files.
ggsave("output/Abundance_Plots_Cancers_Fullrange_Meaned.pdf",qc,width=6,height=9)
# clear memory
gc()


# OK. We want to get dataframes (tables) of:
# 1. non-linear curve-fitting = the rate of synthesis and secretion of newly-made proteins (k)
# obtained from RIA of protein (y) with time (x) data.
# 2. slope of the linear regression fits of how abundance in secretome (y) changes with time (x)

lm_cancers <- filtC_total_Ab %>% nest(data = -Unique_ID) %>%
  mutate(fit=map(data, ~ lm(total_ab ~ time,data=.x)),
         tidied=map(fit, broom::tidy))
#glanced=map(fit,broom::glance),
#augmented=map(fit,broom::augment))
# need to unnest the nested fit data to get slope of linear regression line:
lm_cancers_data <- unnest(lm_cancers,cols = tidied)
# just get the slope info (which is annoyingly called "time")
lm_cancers_data_2 <- dplyr::filter(lm_cancers_data,term=="time")
# extract only the columns we are interested in:
lm_cancers_data_2 <- lm_cancers_data_2[,c(1,5)]
# rename estimate column to slope:
colnames(lm_cancers_data_2)[2] <- "slope"

# let's do the same for the RIA,t non-linear curve-fits now:
nls_cancers <- filtC_mean_RIA %>% nest(data= -Unique_ID) %>%
  mutate(fit = map(data, ~nls(mean_RIA ~ (0 + (1 - 0) * (1 - exp(-k * time))),
                              start=list(k=0.01),data=.x)),
         tidied=map(fit,broom::tidy), augmented=map(fit,broom::augment))
# need to unnest the nested nls fit data to get kinetic rate of synthesis and 
# secretion (k):
nls_cancers_data <- unnest(nls_cancers,cols = tidied)
# extract only the columns we are interested in:
nls_cancers_data_2 <- nls_cancers_data[,c(1,5)]
# rename estimate column to k = rate of synthesis and secretion:
colnames(nls_cancers_data_2)[2] <- "k"

# filter out proteins whose abundance does not increase with time in the secretome.
# These are any with a negative slope value from the linear regression fit of abundance (y)
# with time (x) data (= in lm_norms_data_2 output)
# create a smaller version of lm_cancers_data_2 with the negative slope hits removed:
lm_cancers_data_3 <- dplyr::filter(lm_cancers_data_2,slope>0)

# merge this lm output dataframe with our abundance data ignoring those with negative
# slope data:
filtC_total_Ab_2 <- merge(filtC_total_Ab,lm_cancers_data_3,by="Unique_ID",all=FALSE)
length(unique(filtC_total_Ab_2$Unique_ID))
# down to 560 final filtered hits.
# we can filter the RIA data to only give us RIA of these proteins too as our FINAL
# filtered datasets to work on:
nls_cancers_data_3 <- merge(nls_cancers_data_2,lm_cancers_data_3,by="Unique_ID",all=FALSE)
length(unique(nls_cancers_data_3$Unique_ID))
# it is correct length 562
# merge with RIA data now:
filtC_mean_RIA_2 <- merge(filtC_mean_RIA,nls_cancers_data_3,by="Unique_ID",all=FALSE)
length(unique(filtC_total_Ab_2$Unique_ID))
# correct length = 560.

# FINAL Cancer merged data:
FINAL_cancers <- merge(filtC_total_Ab_2,filtC_mean_RIA_2,by="Unique_ID",all=FALSE)
length(unique(FINAL_cancers$Unique_ID))

# let's just get the slope and fit data now from this:
S_and_k_Cancers <- dplyr::select(FINAL_cancers,c(Unique_ID,slope.x,k))
S_and_k_Cancers <- unique(S_and_k_Cancers)
# rename column names:
colnames(S_and_k_Cancers)=c("Unique_ID","slope_c","k_c")

# We can now run our double plot function on the cancer data now:
double_plot_c <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "blue") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "black", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "darkgreen") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "coral", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "lightgrey", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + ab)
}
# do the multiplots and save:
plot_list_c <- dlply(FINAL_cancers,"Unique_ID",double_plot_c)
C_grob <- marrangeGrob(plot_list_c,nrow=5,ncol=1)
ggsave("Side_by_Side_Multiplots_Cancers.pdf",C_grob,width=5,height=10)

# how about merging slope and k data for each cell type to then look
# at common proteins between the 2 cell types:
S_and_k <- merge(S_and_k_Cancers,S_and_k_norms,by="Unique_ID")
# keep non common protein:
S_and_k_all <- merge(S_and_k_Cancers,S_and_k_norms,by="Unique_ID",all=TRUE)

# let's have a quick look at k vs k across cancer and normals of the
# common proteins:
ggplot(S_and_k,aes(k_n,k_c)) + geom_point() +
  theme(aspect.ratio = 1)
max(S_and_k$k_n)
ggplot(S_and_k,aes(k_n,k_c)) + geom_point() +
  xlim(0,0.03) +ylim(0,0.03) +
  theme(aspect.ratio = 1)

ggplot(S_and_k,aes(k_n,k_c)) + geom_point() +
  xlab(expression(italic(" k") * " (Normal cell secretome)")) +
  ylab(expression(italic(" k") * " (Cancer cell secretome)")) +
  stat_smooth(method = "lm") +
  xlim(0,0.01) +ylim(0,0.01) +
  theme(aspect.ratio = 1)


#svs
ggplot(S_and_k,aes(slope_n,slope_c)) + geom_point() +
  theme(aspect.ratio = 1)
max(S_and_k$slope_n)
ggplot(S_and_k,aes(slope_n,slope_c)) + geom_point() +
  xlim(0,2.0e+7) +ylim(0,2.0e+8) +
  theme(aspect.ratio = 1)

ggplot(S_and_k,aes(slope_n,slope_c)) + geom_point() +
  stat_smooth(method = "lm") +
  xlim(0,2.0e+7) +ylim(0,2.0e+8) +
  theme(aspect.ratio = 1)



d_n<- density(filtN_mean_RIA_2$k)
plot(d_n, col = "blue", main = expression("KDE of" * italic(" k") * " (Normal cell secretome)"), 
     ylab = "Density (d)", xlab = expression(italic("k"))) +
  (abline(v = 0.0235, col = 'grey', lwd=3, lty=2) +
     text(x = 0.05, y = 25,
          "cut-off = 0.025", cex = 0.95))

#geom_vline(xintercept=c(0.025), linetype="dotted") #, 
#type = "l", lty = "dashed"))

d_c<- density(filtC_mean_RIA_2$k)
plot(d_c, col = "red", main = expression("KDE of" * italic(" k") * " (Cancer cell secretome)"),
     xlim=c(0,0.2),
     ylab = "Density (d)", xlab = expression(italic("k"))) +
  (abline(v = 0.015, col = 'grey', lwd=3, lty=2) +
     text(x = 0.04, y = 25,
          "cut-off = 0.015", cex = 0.95))


#lines(d_n, col = "blue", lwd=2)        


# let's plot RIA and abundance plots of these "interesting proteins" to
# discuss in main project report:
# need to filter these from the main dataframe and create individual plot
# objects. Normals first:
prot1=dplyr::filter(FINAL_norms,Unique_ID=="TGFBI - Q15582")
prot2=dplyr::filter(FINAL_cancers,Unique_ID=="TGFBI - Q15582")
#prot3=dplyr::filter(FINAL_norms,Unique_ID=="DPYSL3 - Q14195")
#prot4=dplyr::filter(FINAL_norms,Unique_ID=="IGFBP5 - P24593")

# what are the abundance limits for these?
range(prot1$total_ab)
range(prot2$total_ab)
#range(prot3$total_ab)
#range(prot4$total_ab)
subgroup1=bind_rows(prot1,prot2) #,prot3,prot4)
range(subgroup1$total_ab)
#1924032000

dev.off()
plot1_norm <- ggplot()
double_plot_subset_norm <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 0.1,
             label = expression(italic ("k") * " = 0.1044 h"^-1), size = 3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 9,
             y = 1.75E+9,
             label = expression("slope = 20009668"), size = 2.6, col = "blue")
  #finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

#canc
plot2_canc=ggplot()
double_plot_subset_cancer <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 0.1,
             label = expression(italic ("k") * " = 0.1723 h"^-1), size = 3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 9,
             y = 1.75E+9,
             label = expression("slope = 78066032"), size = 2.6, col = "blue")
  #finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

# the famous multiplot function:
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

protplot1 <- double_plot_subset_norm(prot1)
protplot2 <- double_plot_subset_cancer(prot2)
dev.off()
#protplot3=double_plot_subset(prot3)
#protplot4=double_plot_subset(prot4)
multiplot(protplot1, protplot2, cols=1)
("TGFBI_proteins.pdf") #, protplot1, protplot2)
pdf("output/TGFBI_proteins.pdf")

pdf("output/TGFBI_proteins.pdf")
dev.off()
gc()



#NOW FOR A DIFFERENT PROTEIN: BAG3
prot_bagn <- dplyr::filter(FINAL_norms,Unique_ID=="BAG3 - O95817")
prot_bagc <- dplyr::filter(FINAL_cancers,Unique_ID=="BAG3 - O95817")

bag_group <- bind_rows(prot_bagn,prot_bagc) #,prot3,prot4)
range(bag_group$total_ab)
#2xE+8

plot_bag <- ggplot()
double_plot_bag <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#F0808085") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#F0808085", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#6CA6CD85", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+08) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}
# the famous multiplot function:
multiplot

bag_plotn <- double_plot_bag(prot_bagn)
bag_plotc <- double_plot_bag(prot_bagc)
multiplot(bag_plotn, bag_plotc) #, protplot3, protplot4, cols=1)
pdf(file = "output/BAG3_proteins.pdf")
dev.off()

#NOW FOR A DIFFERENT PROTEIN:CTNNB1
prot_ctnnb1 <- dplyr::filter(FINAL_cancers,Unique_ID=="CTNNB1 - P35222")
range(prot_ctnnb1$total_ab)
plot_ctnnb1=ggplot()
double_plot_ctnnb1 <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 5,
             y = 0.125,
             label = expression(italic ("k") * " = 0.00096 h"^-1), size = 3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + #ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 5,
             y = 2.9E+7,
             label = expression("slope = 352944"), size = 2.6, col = "blue")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}
multiplot

ctnnb1_plot <- double_plot_ctnnb1(prot_ctnnb1)
multiplot(ctnnb1_plot)

#NOW FOR THE TWO PROTEINS COMMON TO BOTH

prot_pdia3n <- dplyr::filter(FINAL_norms,Unique_ID=="PDIA3 - P30101")
prot_pdia3c <- dplyr::filter(FINAL_cancers,Unique_ID=="PDIA3 - P30101")
prot_plecN <- dplyr::filter(FINAL_norms,Unique_ID=="PLEC - Q15149")
prot_plecC <- dplyr::filter(FINAL_cancers,Unique_ID=="PLEC - Q15149")

PDIA3_PLEC_group <- bind_rows(prot_pdia3n, prot_pdia3c, prot_plecN, prot_plecC)
range(PDIA3_PLEC_group$total_ab)
#2.5xE+9
range(prot_pdia3n$total_ab)
range(prot_pdia3c$total_ab)
range(prot_plecC$total_ab)

plot_PDIA3N=ggplot()
double_plot_PDIA3N <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 0.87,
             label = expression(italic ("k") * " = 0.0014 h"^-1), size = 2.3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.2e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 14,
             y = 2.0E+9,
             label = expression("slope = 3079056"), size = 2.3, col = "blue")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

plot_PDIA3C=ggplot()
double_plot_PDIA3C <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 0.87,
             label = expression(italic ("k") * " = 0.0036 h"^-1), size = 2.3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.2e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 14,
             y = 2.5E+8,
             label = expression("slope = 49563558"), size = 2.3, col = "blue")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

#PLEC
plot_PLECn=ggplot()
double_plot_PLECn <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 0.87,
             label = expression(italic ("k") * " = 0.0099 h"^-1), size = 2.3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 14,
             y = 1.8E+9,
             label = expression("slope = 5881075"), size = 2.3, col = "blue")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

plot_PLECc=ggplot()
double_plot_PLECc <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 0.87,
             label = expression(italic ("k") * " = 0.0048 h"^-1), size = 2.3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 14,
             y = 2.5E+8,
             label = expression("slope = 51986164"), size = 2.3, col = "blue")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

# the famous multiplot function:
multiplot

PDIA3_plotn <- double_plot_PDIA3N(prot_pdia3n)
PDIA3_plotc <- double_plot_PDIA3C(prot_pdia3c)
PLEC_PLOTN <- double_plot_PLECn(prot_plecN)
PLEC_PLOTC <- double_plot_PLECc(prot_plecC)
multiplot(PDIA3_plotn, PDIA3_plotc, PLEC_PLOTN, PLEC_PLOTC, cols = 2) #, protplot3, protplot4, cols=1)
pdf("output/FABP5_proteins.pdf")
#since BAG3 and FABP5 plots are quite similar, maybe showcase them together:
multiplot(fabp_plotn, fabp_plotc, bag_plotn, bag_plotc, cols = 1)
pdf("output/FABP5_and_BAG3_proteins.pdf")
dev.off()



#NOW FOR A DIFFERENT PROTEIN (compare TIMP1 with TXN)

prot_timp <- dplyr::filter(FINAL_norms,Unique_ID=="TIMP1 - P01033")
prot_txn <- dplyr::filter(FINAL_norms,Unique_ID=="TXN - P10599")

timp_txn_group <- bind_rows(prot_timp,prot_txn)
range(timp_txn_group$total_ab)
#2xE+9

plot_timp_txn=ggplot()
double_plot_timp <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#FF0000") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#FF0000", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 15,
             y = 0.125,
             label = expression(italic ("k") * " = 0.1852 h"^-1), size = 3, col = "red")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#FF0000") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#FF0000", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 16,
             y = 2.5e+8,
             label = "slope = 67786381", size = 2.5, col = "red")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

double_plot_txn <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#0000FF") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#0000FF", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 15,
             y = 0.125,
             label = expression(italic ("k") * " = 0.0006 h"^-1), size = 3, col = "blue")
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#0000FF") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#0000FF", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1) +
    annotate("text",
             x = 17,
             y = 2.5e+8,
             label = "slope = 88392", size = 2.5, col = "blue")
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}

#

multiplot
# the famous multiplot function:
multiplot 
timp_txn_plot1=double_plot_timp(prot_timp)
timp_txn_plot2=double_plot_txn(prot_txn)
multiplot(timp_txn_plot1, timp_txn_plot2) #, protplot3, protplot4, cols=1)
pdf("output/timp_txn_proteins_with_values.pdf")

#NOW FOR A DIFFERENT PROTEIN: BAG3
prot_egfr <- dplyr::filter(FINAL_cancers,Unique_ID=="EGFR ??? P00533")
prot_mapk1 <- dplyr::filter(FINAL_cancers,Unique_ID=="MAPK1 ??? P28482")

egfr_group <- bind_rows(prot_egfr,prot_mapk1) #,prot3,prot4)
range(bag_group$total_ab)
#2xE+8

plot_egfr=ggplot()
double_plot_egfr <- function(d){
  # RIA-t plot first:
  RIA <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "indianred1") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "indianred1", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    #facet_wrap(~ Unique_ID) +
    ggtitle(d$Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  Ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#1e90ff") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#1e90ff", size = 0.35) +
    xlim(0, 24) + ylim(0, 2.0e+08) +
    ylab("Abundance") + xlab("Time (h)") +
    ggtitle(d$Unique_ID) +
    #facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(RIA + Ab)
}
# the famous multiplot function:
multiplot

egfr_plotc=double_plot_egfr(prot_egfr)
mapk1_plotc=double_plot_egfr(prot_mapk1)
multiplot(egfr_plotc, mapk1_plotc) #, protplot3, protplot4, cols=1)
pdf(file = "output/egfr_mapk1_proteins.pdf")
dev.off()













































#####################testing




# FINAL Normal merged data:
FINAL_norms_new <- merge(filtN_total_Ab_2,filtN_mean_RIA_2,by="Unique_ID",all=FALSE)
Final_norms_new_2 <-  dplyr::filter(FINAL_norms, k > "0.025")
# let's just get the slope and fit data now from this:
S_and_k_norms_2 <- dplyr::select(Final_norms_new_2,c(Unique_ID,slope.x,k))
S_and_k_norms_2 <- unique(S_and_k_norms_2)
# rename column names:
colnames(S_and_k_norms_2) <- c("Unique_ID","slope_n","k_n")



# We can now run our double plot function:
double_plot_new <- function(d){
  # RIA-t plot first:
  ria <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#F0808085") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#F0808085", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#6CA6CD85", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(ria + ab)
}
# do the multiplots and save:
plot_list_N_new=dlply(Final_norms_new_2,"Unique_ID",double_plot_new)
N_grob_new=marrangeGrob(plot_list_N_new,nrow=5,ncol=1)
ggsave("Above_cutoff_Side_by_Side_Multiplots_Normals.pdf",N_grob_new,width=5,height=10)


###  now for those below the cut off
Final_norms_new_3 <-  dplyr::filter(FINAL_norms, k <= "0.025")
# let's just get the slope and fit data now from this:
S_and_k_norms_3 <- dplyr::select(Final_norms_new_3,c(Unique_ID,slope.x,k))
S_and_k_norms_3 <- unique(S_and_k_norms_3)
# rename column names:
colnames(S_and_k_norms_3) <- c("Unique_ID","slope_n","k_n")



# We can now run our double plot function:
double_plot_new_3 <- function(d){
  # RIA-t plot first:
  ria <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#F0808085") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#F0808085", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#6CA6CD85", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(ria + ab)
}
# do the multiplots and save:
plot_list_N_new_3=dlply(Final_norms_new_3,"Unique_ID",double_plot_new_3)
N_grob_new_3=marrangeGrob(plot_list_N_new_3,nrow=5,ncol=1)
ggsave("Below_cutoff_Side_by_Side_Multiplots_Normals.pdf",N_grob_new_3,width=5,height=10)





###cancer
plot(density(filtC_mean_RIA_2$k))


# FINAL Normal merged data:
FINAL_cancers_new=merge(filtC_total_Ab_2,filtC_mean_RIA_2,by="Unique_ID",all=FALSE)
Final_cancers_new_2= dplyr::filter(FINAL_cancers, k > "0.05")
# let's just get the slope and fit data now from this:
S_and_k_cancers_2=dplyr::select(Final_cancers_new_2,c(Unique_ID,slope.x,k))
S_and_k_cancers_2=unique(S_and_k_cancers_2)
# rename column names:
colnames(S_and_k_cancers_2)=c("Unique_ID","slope_c","k_c")



# We can now run our double plot function:
double_plot_new_2c <- function(d){
  # RIA-t plot first:
  ria <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#F0808085") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#F0808085", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#6CA6CD85", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(ria + ab)
}
# do the multiplots and save:
plot_list_C_new=dlply(Final_cancers_new_2,"Unique_ID",double_plot_new_2c)
C_grob_new <- marrangeGrob(plot_list_C_new,nrow=5,ncol=1)
ggsave("output/Above_cutoff_Side_by_Side_Multiplots_Cancers.pdf",C_grob_new,width=5,height=10)


###  now for those below the cut off
Final_cancers_new_3 <- dplyr::filter(FINAL_cancers, k <= "0.05")
# let's just get the slope and fit data now from this:
S_and_k_cancers_3 <- dplyr::select(Final_cancers_new_3,c(Unique_ID,slope.x,k))
S_and_k_cancers_3 <- unique(S_and_k_cancers_3)
# rename column names:
colnames(S_and_k_cancers_3) <- c("Unique_ID","slope_c","k_c")



# We can now run our double plot function:
double_plot_new_3c <- function(d){
  # RIA-t plot first:
  ria <- ggplot(d, aes(time.y, mean_RIA)) +
    geom_point(pch = 23, size = 0.7, fill = "#F0808085") +
    stat_smooth(method = "nls", formula = y ~ (0 + (1 - 0) * (1 - exp(-k * x))),
                method.args = list(start = list(k = 0.01)), se = F,
                fullrange = T, colour = "#F0808085", size = 0.35) +
    xlim(0, 24) + ylim(0, 1) +
    ylab("RIA") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#F0808060", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # next the abundance-t plot:
  ab <- ggplot(d, aes(time.x, total_ab)) +
    geom_point(pch = 23, size = 0.7, fill = "#6CA6CD85") +
    stat_smooth(method = "lm", se = F, fullrange = T, colour = "#6CA6CD85", size = 0.35) +
    xlim(0, 24) +# ylim(0, 2.0e+09) +
    ylab("Abundance") + xlab("Time (h)") +
    facet_wrap(~ Unique_ID) +
    theme(axis.text = element_text(size = 5.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.35, linetype = 'solid'),
          panel.border = element_rect(colour = "black", size = 0.35,fill=NA),
          strip.background = element_rect(fill = "#B0E2FF80", color = "black", size = 0.35),
          strip.text = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 1)
  # finally, tell cowplot of the pairs of plots we want to draw:
  cowplot::plot_grid(ria + ab)
}
# do the multiplots and save:
plot_list_C_new_3c <- dlply(Final_cancers_new_3,"Unique_ID",double_plot_new_3c)
C_grob_new_3 <- marrangeGrob(plot_list_C_new_3c,nrow=5,ncol=1)
ggsave("output/Below_cutoff_Side_by_Side_Multiplots_Cancers.pdf",C_grob_new_3,width=5,height=10)



# 
# organism= "Horg.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
# ?repositories
# library(organism, character.only = TRUE)
# library(GO.db)
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("org.Hs.eg.db", force = TRUE)

# data(filtN, package = "DOSE")

# citation(package = "tidyverse")

#relax some of the filters: only use petides that were at 6 and 24 hours.
#what to do next

#gen in cancenorm cells sec/exch of label from light to evy is quicker

# citation(package = "tidyverse")


##
