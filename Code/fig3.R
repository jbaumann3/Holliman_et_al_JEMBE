# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# Generating Figure 3: Fv/Fm, Symbiont Density, and Chlorophyll a 

# Last Updated: 4/9/25

#################### GETTING STARTED ####################

# set wd
setwd("~/Library/CloudStorage/OneDrive-BowdoinCollege/Desktop/Holliman_et_al_JEMBE_edited042025/Data")

# load libraries
library(tidyverse)

pd=position_dodge(width=0.6)

##### Fv/Fm

# Read in data 
pam <- read.csv("pam.csv")
pam <- pam[c(1:568), c(1:14)]

# Make some new columns
pam$meanfvfm <- (pam$Fv.Fm_1 + pam$Fv.Fm_2 + pam$Fv.Fm_3)/3

pam$population <- ifelse(grepl("NC", pam$coral_ID), "NC", 
                         ifelse(grepl("MA", pam$coral_ID), "MA", "RI"))

pam_end <- pam[c(421:488,490:505),]

end_stats <- pam_end %>%
  group_by(population, treatment_temp) %>%
  dplyr::summarize(mean=mean(meanfvfm), sd=sd(meanfvfm), n=n(), se=(sd/sqrt(n)))

end_stats$population <- factor(end_stats$population, levels = c("MA", "RI", "NC"))
end_stats$treatment_temp <- as.factor(end_stats$treatment_temp)

pam_end$population <- factor(pam_end$population, levels = c("MA", "RI", "NC"))
pam_end$treatment_temp <- as.factor(pam_end$treatment_temp)

A <- ggplot(end_stats, aes(x=treatment_temp, y=mean, color=population)) + 
  geom_point(size=3, position=pd)+
  geom_point(data=pam_end, aes(x=treatment_temp, y=meanfvfm, color=population), size=3, alpha=0.2, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(title = "A", y = expression(F[v]/F[m]), x = "Treatment (°C)", color = "Population")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  geom_text(aes(x = 0.77, y = 0, label = "Aab"), size = 3, color="black")+ 
  geom_text(aes(x = 1.0, y = 0, label = "Aab"), size = 3, color="black")+
  geom_text(aes(x = 1.23, y = 0, label = "Bab"), size = 3, color="black")+
  ##
  geom_text(aes(x = 1.77, y = 0, label = "Aa"), size = 3, color="black")+ 
  geom_text(aes(x = 2.0, y = 0, label = "Aa"), size = 3, color="black")+
  geom_text(aes(x = 2.23, y = 0, label = "Ba"), size = 3, color="black")+
  ##
  geom_text(aes(x = 2.77, y = 0, label = "Ab"), size = 3, color="black")+ 
  geom_text(aes(x = 3.0, y = 0, label = "Ab"), size = 3, color="black")+
  geom_text(aes(x = 3.23, y = 0, label = "Bb"), size = 3, color="black")+
  ##
  geom_text(aes(x = 3.77, y = 0, label = "Ac"), size = 3, color="black")+ 
  geom_text(aes(x = 4.0, y = 0, label = "Ac"), size = 3, color="black")+
  geom_text(aes(x = 4.23, y = 0, label = "Bc"), size = 3, color="black")

A

##### Symbiont Density 

# Load in data
sym <- read.csv("sym_count.csv")

options(scipen=0)

#New columns to calculate 
sym$total_sym_number = sym$mean_sym*(sym$total_solution_volume_ul/sym$hemo_vol_ul)

sym
# Create new columns
sym$sym_sci <- sym$total_sym_number / 1000000

# Make new dfs
sym_exp <- sym[c(25:108), c(1,2,7,14,15)]
# delete MAH (outlier) from sym
sym_exp <- sym_exp[-c(56),]
sym_T0 <- sym[c(1:24),c(1,2,7,14,15)]

sym_exp$population <- ifelse(grepl("NC", sym_exp$coral_ID), "NC", 
                             ifelse(grepl("MA", sym_exp$coral_ID), "MA","RI"))

sym_T0$population <- ifelse(grepl("NC", sym_T0$coral_ID), "NC", 
                            ifelse(grepl("MA", sym_T0$coral_ID), "MA","RI"))

Tend_sym_density <- sym_exp %>%
  group_by(treatment_temp, population) %>%
  summarize(mean=mean(sym_sci), sd=sd(sym_sci), n=n(), se=(sd/sqrt(n)))

Tend_sym_density$timepoint <- "End"

T0_sym_density <- sym_T0 %>%
  group_by(treatment_temp, population) %>%
  summarize(mean=mean(sym_sci), sd=sd(sym_sci), n=n(), se=(sd/sqrt(n)))

T0_sym_density$timepoint <- "Start"

# Combine tables 
bleach <- rbind(Tend_sym_density, T0_sym_density)

bleach$treatment_temp[is.na(bleach$treatment_temp)] <- "T0 Control"

bleach$timepoint <- factor(bleach$timepoint, levels = c("Start", "End"))
bleach$population <- factor(bleach$population, levels = c("MA", "RI", "NC"))
bleach$treatment_temp <- factor(bleach$treatment_temp, levels = c("T0 Control", "18", "22", "28", "32"))

sym$treatment_temp[is.na(sym$treatment_temp)] <- "T0 Control"

sym$population <- ifelse(grepl("NC", sym$coral_ID), "NC", 
                         ifelse(grepl("MA", sym$coral_ID), "MA","RI"))

sym$population <- factor(sym$population, levels = c("MA", "RI", "NC"))
sym$treatment_temp <- factor(sym$treatment_temp, levels = c("T0 Control", "18", "22", "28", "32"))

pd=position_dodge(width=0.4)

# Graph this
ggplot(bleach, aes(x=population, y=mean, color=treatment_temp)) + 
  geom_point(size=3, position=pd)+
  geom_point(data=sym, aes(x=population, y=sym_sci, color=treatment_temp), size=3, alpha=0.2, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  theme_classic()+
  ylab("Symbiont Density (10^6 cells / cm^2)")+
  xlab("Population")+
  labs(color="Treatment (°C)")+
  scale_color_manual(labels=c("T0 Control" = "Control", "18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"), values=c("grey57", "purple4","cyan4", "springgreen4", "goldenrod"))

anova_sym <- aov(sym_sci ~ treatment_temp*population, data = sym)
summary(anova_sym)
tukey_sym<-TukeyHSD(anova_sym)
tukey_sym

####### SYM but no control

bleach_end <- bleach[c(1:12),]
sym_end <- sym[-c(1:24),]

B <- ggplot(bleach_end, aes(x=treatment_temp, y=mean, color=population)) + 
  geom_point(size=3, position=pd)+
  geom_point(data=sym_end, aes(x=treatment_temp, y=sym_sci, color=population), size=3, alpha=0.2, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(title = "B", y = expression("Sym Density" ~ (10^6 ~ "cells" ~ "/" ~ cm^2)), x = "Treatment (°C)", color = "Population")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  geom_text(aes(x = 0.77, y = 0, label = "ab"), size = 3, color="black")+ 
  geom_text(aes(x = 1.0, y = 0, label = "ab"), size = 3, color="black")+
  geom_text(aes(x = 1.23, y = 0, label = "ab"), size = 3, color="black")+
  ##
  geom_text(aes(x = 1.77, y = 0, label = "a"), size = 3, color="black")+ 
  geom_text(aes(x = 2.0, y = 0, label = "a"), size = 3, color="black")+
  geom_text(aes(x = 2.23, y = 0, label = "a"), size = 3, color="black")+
  ##
  geom_text(aes(x = 2.77, y = 0, label = "ab"), size = 3, color="black")+ 
  geom_text(aes(x = 3.0, y = 0, label = "ab"), size = 3, color="black")+
  geom_text(aes(x = 3.23, y = 0, label = "ab"), size = 3, color="black")+
  ##
  geom_text(aes(x = 3.77, y = 0, label = "b"), size = 3, color="black")+ 
  geom_text(aes(x = 4.0, y = 0, label = "b"), size = 3, color="black")+
  geom_text(aes(x = 4.23, y = 0, label = "b"), size = 3, color="black")

plot(B)

anova_sym2 <- aov(sym_sci ~ treatment_temp*population, data = sym_end)
summary(anova_sym2)
tukey_sym2<-TukeyHSD(anova_sym2)
tukey_sym2

##### Chlorophyll a 

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


##### Minus Control 
chla_area_stats_nocontrol <- chla_area_stats[-c(5,10,15),]
chla_nocontrol <- subset(chla,treatment_temp!="Control")

C <- ggplot(chla_area_stats_nocontrol, aes(x=treatment_temp, y=mean, color=population))+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  geom_point(size=3, position=pd)+
  geom_point(data=chla_nocontrol, aes(x=treatment_temp, y=chla_density, color=population), size=3, alpha=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(title = "C", y = expression("Chlorophyll a" ~ (pmol / cm^2)), x = "Treatment (°C)", color = "Population")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  geom_text(aes(x = 0.77, y = 0, label = "a"), size = 3, color="black")+ 
  geom_text(aes(x = 1.0, y = 0, label = "a"), size = 3, color="black")+
  geom_text(aes(x = 1.23, y = 0, label = "a"), size = 3, color="black")+
  ##
  geom_text(aes(x = 1.77, y = 0, label = "a"), size = 3, color="black")+ 
  geom_text(aes(x = 2.0, y = 0, label = "a"), size = 3, color="black")+
  geom_text(aes(x = 2.23, y = 0, label = "a"), size = 3, color="black")+
  ##
  geom_text(aes(x = 2.77, y = 0, label = "a"), size = 3, color="black")+ 
  geom_text(aes(x = 3.0, y = 0, label = "a"), size = 3, color="black")+
  geom_text(aes(x = 3.23, y = 0, label = "a"), size = 3, color="black")+
  ##
  geom_text(aes(x = 3.77, y = 0, label = "a"), size = 3, color="black")+ 
  geom_text(aes(x = 4.0, y = 0, label = "a"), size = 3, color="black")+
  geom_text(aes(x = 4.23, y = 0, label = "a"), size = 3, color="black")


C

######### MAKING THE FINAL FIGURE
library(patchwork)

final_plot <- A + B + C + 
  plot_layout(nrow = 3) 

final_plot
