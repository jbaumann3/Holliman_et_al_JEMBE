# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# HPLC and Chlorophyll a

# Last Updated: 6/4/23

#################### GETTING STARTED ####################

# Set working directory
setwd("~/Desktop/JEMBE Paper Data Analysis/Data")

# Load libraries 
library(lubridate)
library (tidyverse)
library(xts) 
library(nls.multstart)
library(broom)
library(rTPC)
library(ggsci)
library(modelr)
library(performance)
library(patchwork)
library(car)
library(minpack.lm)
library(nlstools)
library(MuMIn)
library(ggrepel)

#################### ORGANIZE DATA ####################

# Read in data
chla <- read.csv("chla_hplc.csv")

#delete MAV, outlier
chla <- chla[-c(31),]

# Take out LAHP and DAHP 
hp <- chla[c(1:6),]

# Delete LAHP and DAHP from chla
chla <- chla[-c(1:6),]

chla$treatment_temp[is.na(chla$treatment_temp)] <- "Control"

#################### REPLICATES ####################

# Average replicates FIRST, then do stats
# Replicates are NCV, NCW, NCX, and RIF

rep <- chla[c(55:60,70,71),]
rep_mean <- rep %>%
  group_by(coral_ID) %>%
  summarize(mean=mean(peak_area), sd=sd(peak_area), n=n(), se=(sd/sqrt(n)))

chla_no_reps <- chla[-c(55:60,70,71),]

rep <- rep[c(1,3,5,7),]

rep_mean_only <- rep_mean[,c(2)]
rep <- cbind(rep_mean_only, rep)

rep<- rep[,-c(5)]
rep$peak_area <- rep$mean
rep<- rep[,-c(1)]

chla <- rbind(rep, chla_no_reps)

#################### MAKE A FEW MORE COLUMNS ####################

# HPLC samples contained 18% of the total symbionts for each coral
# Therefore, chla/cell = peak_area / .18(sym_number_sample)

# For units of picomoles, need to multiply injection volume by 0.888

chla$chla_per_cell <- ((chla$peak_area*0.888) / (.18 * (chla$sym_number_sample)))

# Need to add population column 

chla$population <- ifelse(grepl("NC", chla$coral_ID), "NC", 
                          ifelse(grepl("MA", chla$coral_ID), "MA", 
                                 ifelse(grepl("RI", chla$coral_ID), "RI", "blank")))

# Chla density 
chla$chla_density <- ((chla$peak_area*0.888) / (.18 * (chla$surface_area)))

#chla$treatment_temp <- factor(chla$treatment_temp, levels = c("NA", "18", "22", "28","32"))

#################### SUMMARY STATS ####################

chla_cell_stats <- chla %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(chla_per_cell), sd=sd(chla_per_cell), n=n(), se=(sd/sqrt(n)))

chla_area_stats <- chla %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(chla_density), sd=sd(chla_density), n=n(), se=(sd/sqrt(n)))

#################### MAKE A QUICK FIGURE ####################

pd=position_dodge(width=0.4)

chla$population <- factor(chla$population, levels = c("MA", "RI", "NC"))
chla_cell_stats$population <- factor(chla_cell_stats$population, levels = c("MA", "RI", "NC"))
chla_area_stats$population <- factor(chla_area_stats$population, levels = c("MA", "RI", "NC"))
chla_cell_stats$treatment_temp <- as.factor(chla_cell_stats$treatment_temp)
chla_cell_stats$treatment_temp <- factor(chla_cell_stats$treatment_temp, levels = c("Control","18","22","28","32"))

ggplot(chla_area_stats, aes(x=treatment_temp, y=mean, color=population))+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  geom_point(size=3, position=pd)+
  geom_point(data=chla, aes(x=treatment_temp, y=chla_density, color=population), size=3, alpha=0.2, position=pd)+
  theme_classic()+
  ylab("Chlorophyll a (pmol/cell)")+
  xlab("Treatment (°C)")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  scale_x_discrete(labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  labs(color="Population")

##### Minus Control 
chla_area_stats_nocontrol <- chla_area_stats[-c(5,10,15),]
chla_nocontrol <- subset(chla,treatment_temp!="Control")

ggplot(chla_area_stats_nocontrol, aes(x=treatment_temp, y=mean, color=population))+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  geom_point(size=3, position=pd)+
  #geom_point(data=chla_nocontrol, aes(x=treatment_temp, y=chla_density, color=population), size=3, alpha=0.2, position=pd)+
  theme_classic(base_size = 12)+
  ylab("Chlorophyll a (pmol/cm^2)")+
  xlab("Treatment (°C)")+
  labs(color="Population")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))



  
#################### SOME MORE STATS ####################

anova1 <- aov(chla_per_cell ~ treatment_temp*population, data = chla)
summary(anova1)
tukey1<-TukeyHSD(anova1)
tukey1

anova2 <- aov(chla_density ~ treatment_temp*population, data = chla)
summary(anova2)
tukey2<-TukeyHSD(anova2)
tukey2


#stats without control!
chla2<- chla %>%
  filter(treatment_temp!='Control')

anova3 <- aov(chla_per_cell ~ treatment_temp*population, data = chla2)
summary(anova3)
tukey3<-TukeyHSD(anova3)
tukey3

