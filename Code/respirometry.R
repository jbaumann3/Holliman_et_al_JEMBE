# DEVA HOLLIMAN
# BOWDOIN COLLEGE
# HONORS PROJECT 2022-23
# DATA ANALYSIS

# Last Updated: 11/17/23

#################### GETTING STARTED ####################

# Set working directory
#setwd("~/Desktop/Honors Project/Winter 2023 Analysis/Data")

# Load libraries
library(lubridate)
library(ggpubr)
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
library(Matrix)

######################################################
#################### RESPIROMETRY ####################
######################################################

#################### ORGANIZE DATA ####################

# Read in data
resp <- read.csv("respirometry.csv")
resp <- resp[,c(1:24)] # Delete empty row 

sym <- read.csv("sym_count.csv")
sym$sym_sci <- sym$sym.cm2 / 1000000 # Convert to scientific notation (# x 10^6 cells/cm2)
sym_exp <- sym[c(25:108), c(1,2,7,15,16)] # Symbiont density of just experimental corals
sym_T0 <- sym[c(1:24),c(1,2,7,15,16)] # Symbiont density of just initial control corals

# Make some new columns
resp$population <- ifelse(grepl("NC", resp$coral_ID), "NC", 
                          ifelse(grepl("MA", resp$coral_ID), "MA", 
                                 ifelse(grepl("RI", resp$coral_ID), "RI", "blank")))

resp$coral_volume <- resp$empty_chamber_volume - resp$full_chamber_volume # Unit is cm3
resp$coral_volume[is.na(resp$coral_volume)] <- 0

resp$mean_bw <- (resp$BW_1 + resp$BW_2 + resp$BW_3)/3 # Unit is grams

resp$coral_ID[is.na(resp$coral_ID)] <- "blank"

#################### Separate dfs by timepoint ####################

## Going to do all the calculations, and then compare duplicates. 

# INITIAL 
resp_initial <- resp[(1:101),]

# Make start/stop times seconds
resp_initial$dark_time_start <- ms(resp_initial$dark_time_start)
resp_initial$dark_time_start <- as.duration(resp_initial$dark_time_start)
resp_initial$dark_time_start <- as.numeric(resp_initial$dark_time_start)

resp_initial$dark_time_stop <- ms(resp_initial$dark_time_stop)
resp_initial$dark_time_stop <- as.duration(resp_initial$dark_time_stop)
resp_initial$dark_time_stop <- as.numeric(resp_initial$dark_time_stop)

resp_initial$light_time_start <- ms(resp_initial$light_time_start)
resp_initial$light_time_start <- as.duration(resp_initial$light_time_start)
resp_initial$light_time_start <- as.numeric(resp_initial$light_time_start)

resp_initial$light_time_stop <- ms(resp_initial$light_time_stop)
resp_initial$light_time_stop <- as.duration(resp_initial$light_time_stop)
resp_initial$light_time_stop <- as.numeric(resp_initial$light_time_stop)

# Calculate resp_initial respiration rate 
resp_initial$dark_time_elapsed_hr <- (resp_initial$dark_time_stop - resp_initial$dark_time_start) / 3600
resp_initial$dark_DO_produced_umol <- (resp_initial$dark_DO_final - resp_initial$dark_DO_initial) * (1000/31.999)

# If coral volume = 0, then __. Else, __.
resp_initial$resp_rate <- ifelse(grepl(0, resp_initial$coral_volume), (resp_initial$dark_DO_produced_umol * (1/resp_initial$empty_chamber_volume) * (1/resp_initial$dark_time_elapsed_hr)) ,
                                 (resp_initial$dark_DO_produced_umol * (1/resp_initial$mean_bw) * (1/resp_initial$full_chamber_volume) * (1/resp_initial$dark_time_elapsed_hr)))

# Calculate resp_initial photosynthesis rate 
resp_initial$light_time_elapsed_hr <- (resp_initial$light_time_stop - resp_initial$light_time_start) / 3600
resp_initial$light_DO_produced_umol <- (resp_initial$light_DO_final - resp_initial$light_DO_initial) * (1000/31.999)

resp_initial$photo_rate <- ifelse(grepl(0, resp_initial$coral_volume), (resp_initial$light_DO_produced_umol * (1/resp_initial$empty_chamber_volume) * (1/resp_initial$light_time_elapsed_hr)) ,
                                  (resp_initial$light_DO_produced_umol * (1/resp_initial$mean_bw) * (1/resp_initial$full_chamber_volume) * (1/resp_initial$light_time_elapsed_hr)))

# FINAL TREATMENT
resp_final_treatments <- resp[c(102:216, 221:222, 227),]

# Make start/stop times seconds
resp_final_treatments$dark_time_start <- ms(resp_final_treatments$dark_time_start)
resp_final_treatments$dark_time_start <- as.duration(resp_final_treatments$dark_time_start)
resp_final_treatments$dark_time_start <- as.numeric(resp_final_treatments$dark_time_start)

resp_final_treatments$dark_time_stop <- hms(resp_final_treatments$dark_time_stop)
resp_final_treatments$dark_time_stop <- as.duration(resp_final_treatments$dark_time_stop)
resp_final_treatments$dark_time_stop <- as.numeric(resp_final_treatments$dark_time_stop)

resp_final_treatments$light_time_start <- ms(resp_final_treatments$light_time_start)
resp_final_treatments$light_time_start <- as.duration(resp_final_treatments$light_time_start)
resp_final_treatments$light_time_start <- as.numeric(resp_final_treatments$light_time_start)

resp_final_treatments$light_time_stop <- hms(resp_final_treatments$light_time_stop)
resp_final_treatments$light_time_stop <- as.duration(resp_final_treatments$light_time_stop)
resp_final_treatments$light_time_stop <- as.numeric(resp_final_treatments$light_time_stop)

# Calculate resp_final_treatments respiration rate 
resp_final_treatments$dark_time_elapsed_hr <- (resp_final_treatments$dark_time_stop - resp_final_treatments$dark_time_start) / 3600
resp_final_treatments$dark_DO_produced_umol <- (resp_final_treatments$dark_DO_final - resp_final_treatments$dark_DO_initial) * (1000/31.999)

resp_final_treatments$resp_rate <- ifelse(grepl(0, resp_final_treatments$coral_volume), (resp_final_treatments$dark_DO_produced_umol * (1/resp_final_treatments$empty_chamber_volume) * (1/resp_final_treatments$dark_time_elapsed_hr)) ,
                                          (resp_final_treatments$dark_DO_produced_umol * (1/resp_final_treatments$mean_bw) * (1/resp_final_treatments$full_chamber_volume) * (1/resp_final_treatments$dark_time_elapsed_hr)))

# If there is missing data:
resp_final_treatments$light_time_start[is.na(resp_final_treatments$light_time_start)] <- 0
resp_final_treatments$light_time_stop[is.na(resp_final_treatments$light_time_stop)] <- 4000

# Calculate resp_final_treatments photosynthesis rate                                       
resp_final_treatments$light_time_elapsed_hr <- (resp_final_treatments$light_time_stop - resp_final_treatments$light_time_start) / 3600
resp_final_treatments$light_DO_produced_umol <- (resp_final_treatments$light_DO_final - resp_final_treatments$light_DO_initial) * (1000/31.999)

resp_final_treatments$photo_rate <- ifelse(grepl(0, resp_final_treatments$coral_volume), (resp_final_treatments$light_DO_produced_umol * (1/resp_final_treatments$empty_chamber_volume) * (1/resp_final_treatments$light_time_elapsed_hr)) ,
                                           (resp_final_treatments$light_DO_produced_umol * (1/resp_final_treatments$mean_bw) * (1/resp_final_treatments$full_chamber_volume) * (1/resp_final_treatments$light_time_elapsed_hr)))

# FINAL AMBIENT
resp_final_ambient <- resp[c(186:197, 204:293),]

# Make start/stop times seconds
resp_final_ambient$dark_time_start <- ms(resp_final_ambient$dark_time_start)
resp_final_ambient$dark_time_start <- as.duration(resp_final_ambient$dark_time_start)
resp_final_ambient$dark_time_start <- as.numeric(resp_final_ambient$dark_time_start)

resp_final_ambient$dark_time_stop <- hms(resp_final_ambient$dark_time_stop)
resp_final_ambient$dark_time_stop <- as.duration(resp_final_ambient$dark_time_stop)
resp_final_ambient$dark_time_stop <- as.numeric(resp_final_ambient$dark_time_stop)

resp_final_ambient$light_time_start <- ms(resp_final_ambient$light_time_start)
resp_final_ambient$light_time_start <- as.duration(resp_final_ambient$light_time_start)
resp_final_ambient$light_time_start <- as.numeric(resp_final_ambient$light_time_start)

resp_final_ambient$light_time_stop <- hms(resp_final_ambient$light_time_stop)
resp_final_ambient$light_time_stop <- as.duration(resp_final_ambient$light_time_stop)
resp_final_ambient$light_time_stop <- as.numeric(resp_final_ambient$light_time_stop)

# Calculate resp_final_ambient respiration rate 
resp_final_ambient$dark_time_elapsed_hr <- (resp_final_ambient$dark_time_stop - resp_final_ambient$dark_time_start) / 3600
resp_final_ambient$dark_DO_produced_umol <- (resp_final_ambient$dark_DO_final - resp_final_ambient$dark_DO_initial) * (1000/31.999)

resp_final_ambient$resp_rate <- ifelse(grepl(0, resp_final_ambient$coral_volume), (resp_final_ambient$dark_DO_produced_umol * (1/resp_final_ambient$empty_chamber_volume) * (1/resp_final_ambient$dark_time_elapsed_hr)) ,
                                       (resp_final_ambient$dark_DO_produced_umol * (1/resp_final_ambient$mean_bw) * (1/resp_final_ambient$full_chamber_volume) * (1/resp_final_ambient$dark_time_elapsed_hr)))

# Fix missing data:
resp_final_ambient$light_time_stop[is.na(resp_final_ambient$light_time_stop)] <- 4000

# Calculate resp_final_ambient photosynthesis rate       
resp_final_ambient$light_time_elapsed_hr <- (resp_final_ambient$light_time_stop - resp_final_ambient$light_time_start) / 3600
resp_final_ambient$light_DO_produced_umol <- (resp_final_ambient$light_DO_final - resp_final_ambient$light_DO_initial) * (1000/31.999)

resp_final_ambient$photo_rate <- ifelse(grepl(0, resp_final_ambient$coral_volume), (resp_final_ambient$light_DO_produced_umol * (1/resp_final_ambient$empty_chamber_volume) * (1/resp_final_ambient$light_time_elapsed_hr)) ,
                                        (resp_final_ambient$light_DO_produced_umol * (1/resp_final_ambient$mean_bw) * (1/resp_final_ambient$full_chamber_volume) * (1/resp_final_ambient$light_time_elapsed_hr)))

#################### BLANK CORRECTION ####################

# Blank-correct data, then compare duplicates. 

# INITIAL

#step 1: make df for corals and df for blanks
ri_blanks<-subset(resp_initial, coral_ID=='blank')
ri_corals<-subset(resp_initial, coral_ID!='blank')

#step 2: rename columns
ri_blanks<- ri_blanks %>%
  rename('blank_resp_rate'='resp_rate', 'blank_photo_rate'='photo_rate')

#step 3: Mean blank values per trial
mean_blanks0<- ri_blanks %>%
  group_by(trial) %>%
  summarize(mean_blank_R= mean (blank_resp_rate), mean_blank_P=mean(blank_photo_rate))

#step 4: merge blank P and R to corals
ri_2<-merge(ri_corals, mean_blanks0, by='trial')
head(ri_2)

#step 5: create corrected P/R columns
ri_2$p_correct <- ri_2$photo_rate - ri_2$mean_blank_P
ri_2$r_correct <- ri_2$resp_rate - ri_2$mean_blank_R

# FINAL TREATMENTS

#step 1: make df for corals and df for blanks
rft_blanks<-subset(resp_final_treatments, coral_ID=='blank')
rft_corals<-subset(resp_final_treatments, coral_ID!='blank')

#step 2: rename columns
rft_blanks<- rft_blanks %>%
  rename('blank_resp_rate'='resp_rate', 'blank_photo_rate'='photo_rate')

#step 3: Mean blank values per trial
mean_blanks<- rft_blanks %>%
  group_by(trial) %>%
  summarize(mean_blank_R= mean (blank_resp_rate), mean_blank_P=mean(blank_photo_rate))

#step 4: merge blank P and R to corals
rft_2<-merge(rft_corals, mean_blanks, by='trial')
head(rft_2)

#step 5: create corrected P/R columns
rft_2$p_correct <- rft_2$photo_rate - rft_2$mean_blank_P
rft_2$r_correct <- rft_2$resp_rate - rft_2$mean_blank_R

# FINAL AMBIENT

#step 1: make df for corals and df for blanks
rfa_blanks<-subset(resp_final_ambient, coral_ID=='blank')
rfa_corals<-subset(resp_final_ambient, coral_ID!='blank')

#step 2: rename columns
rfa_blanks<- rfa_blanks %>%
  rename('blank_resp_rate'='resp_rate', 'blank_photo_rate'='photo_rate')

#step 3: Mean blank values per trial
mean_blanks2<- rfa_blanks %>%
  group_by(trial) %>%
  summarize(mean_blank_R= mean (blank_resp_rate), mean_blank_P=mean(blank_photo_rate))

#step 4: merge blank P and R to corals
rfa_2<-merge(rfa_corals, mean_blanks2, by='trial')
head(rfa_2)

#step 5: create corrected P/R columns
rfa_2$p_correct <- rfa_2$photo_rate - rfa_2$mean_blank_P
rfa_2$r_correct <- rfa_2$resp_rate - rfa_2$mean_blank_R

#################### DUPLICATES ####################

# INITIAL
initial_dup <- ri_2 %>% 
  group_by(coral_ID) %>% 
  filter(n() != 1) # Return coral_IDs which have more than one row of data

ri_2 <- ri_2 %>% 
  group_by(coral_ID) %>% 
  filter(coral_ID != "NCK")

# FINAL TREATMENTS
treat_dup <- rft_2 %>% 
  group_by(coral_ID) %>% 
  filter(n() != 1) # Return coral_IDs which have more than one row of data (NCK)

treat_dup2 <- treat_dup %>%
  group_by(coral_ID) %>%
  summarize(r_correct_mean= mean(r_correct), p_correct_mean=mean(p_correct))

rft_2 <- rft_2 %>% 
  group_by(coral_ID) %>% 
  filter(coral_ID != "NCK")

# FINAL AMBIENT

ambient_dup <- rfa_2 %>% 
  group_by(coral_ID) %>% 
  filter(n() != 1) # No duplicates here! 

rfa_2 <- rfa_2 %>% 
  group_by(coral_ID) %>% 
  filter(coral_ID != "NCK")

# RECAP:
# ri_2 <- initial
# rft_2 <- final treatment
# rfa_2 <- final ambient

# Want to remove MAH because it was a huge outlier. 
ri_2 <- ri_2 %>% 
  group_by(coral_ID) %>% 
  filter(coral_ID != "MAH")

rft_2 <- rft_2 %>% 
  group_by(coral_ID) %>% 
  filter(coral_ID != "MAH")

rfa_2 <- rfa_2 %>% 
  group_by(coral_ID) %>% 
  filter(coral_ID != "MAH")

#################### FIGURES ####################

# INITIAL 

# Convert wide to long

ri_2$Rdark <- ri_2$r_correct
ri_2$Pnet <- ri_2$p_correct
ri_2$Pgross <- (ri_2$p_correct - ri_2$r_correct)

ri_3 <- ri_2 %>% pivot_longer(cols=c('Pgross', 'Rdark','Pnet'),
                                names_to='light_treatment',
                                values_to='rate')

ri_sum<- ri_3 %>%
  group_by(population, light_treatment) %>%
  summarize(mean=mean(rate), sd=sd(rate), n=n(), se=sd/sqrt(n))

ri_3$population <- factor(ri_3$population, levels = c("MA", "RI", "NC"))
ri_sum$population <- factor(ri_sum$population, levels = c("MA", "RI", "NC"))

pd=position_dodge(width=0.4)

ri_gross <- ggplot()+
  geom_point(aes(x=population, y=mean, color=light_treatment), ri_sum, size=3, position=pd)+
  geom_point(aes(x=population, y=rate, color=light_treatment), ri_3, alpha=0.1, size=3, position=pd)+
  geom_errorbar(aes(x=population, ymin=mean-se, ymax=mean+se, color=light_treatment), ri_sum, width=0.2, position=pd)+
  theme_classic()+
  ylab(expression(ΔDO~(umol~O[2]~g~coral^{"-1"}~mL~H[2]~O^{"-1"}~hr^{"-1"}))) +
  xlab("Metabolic Rate")+
  geom_hline(yintercept = 0, alpha=0.2, linetype="dashed")+
  scale_color_manual(values=c("#1F449C","#A8B6CC","#F05039"), name="Metabolic Rate", labels=c("Rdark"=expression(R[dark]), "Pgross"=expression(P[gross]), "Pnet"=expression(P[net])))

ri_gross

# NET P

# FINAL TREATMENTS

rft_2$treatment_temp <- as.character(rft_2$treatment_temp)
rft_2$population <- factor(rft_2$population, levels = c("MA", "RI", "NC"))

rft_2$Photosynthesis <- rft_2$p_correct
rft_2$Respiration <- rft_2$r_correct

rft_3 <- rft_2 %>% pivot_longer(cols=c('Photosynthesis', 'Respiration'),
                              names_to='light_treatment',
                              values_to='rate')

rft_3$population <- factor(rft_3$population, levels = c("MA", "RI", "NC"))
rft_3$light_treatment <- factor(rft_3$light_treatment, levels = c("Photosynthesis", "Respiration"))

rft_sum<- rft_3 %>%
  group_by(population, light_treatment, treatment_temp) %>%
  summarize(mean=mean(rate), sd=sd(rate), n=n(), se=sd/sqrt(n))

rft <- ggplot()+
  geom_point(aes(x=treatment_temp, y=mean, color=population, shape=light_treatment), rft_sum, size=3)+
  geom_point(aes(x=treatment_temp, y=rate, color=population, shape=light_treatment), rft_3, alpha=0.2, size=3)+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), rft_sum, width=0.2)+
  theme_classic()+
  ylab(expression(ΔDO~(umol~O[2]~g~coral^{"-1"}~mL~H[2]~O^{"-1"}~hr^{"-1"}))) +
  xlab("Temperature Treatment (ºC)")+
  geom_hline(yintercept = 0, alpha=0.2, linetype="dashed")+
  facet_grid(cols=vars(population))+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  guides(color = FALSE)+
  ggtitle("A")+
  scale_shape(name="Metabolic Rate", labels=c("Respiration"=expression(R[dark]), "Photosynthesis"=expression(P[net])))

rft

# FINAL AMBIENT

rfa_2$treatment_temp <- as.character(rfa_2$treatment_temp)
rfa_2$population <- factor(rfa_2$population, levels = c("MA", "RI", "NC"))

rfa_2$Photosynthesis <- rfa_2$p_correct
rfa_2$Respiration <- rfa_2$r_correct

rfa_3 <- rfa_2 %>% pivot_longer(cols=c('Photosynthesis', 'Respiration'),
                                names_to='light_treatment',
                                values_to='rate')

rfa_3$population <- factor(rfa_3$population, levels = c("MA", "RI", "NC"))
rfa_3$light_treatment <- factor(rfa_3$light_treatment, levels = c("Photosynthesis", "Respiration"))

rfa_sum<- rfa_3 %>%
  group_by(population, light_treatment, treatment_temp) %>%
  summarize(mean=mean(rate), sd=sd(rate), n=n(), se=sd/sqrt(n))

rfa <- ggplot()+
  geom_point(aes(x=treatment_temp, y=mean, color=population, shape=light_treatment), rfa_sum, size=3)+
  geom_point(aes(x=treatment_temp, y=rate, color=population, shape=light_treatment), rfa_3, alpha=0.2, size=3)+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), rfa_sum, width=0.2)+
  theme_classic()+
  ylab(expression(ΔDO~(umol~O[2]~g~coral^{"-1"}~mL~H[2]~O^{"-1"}~hr^{"-1"}))) +
  xlab("Temperature Treatment (ºC)")+
  geom_hline(yintercept = 0, alpha=0.2, linetype="dashed")+
  facet_grid(cols=vars(population))+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  guides(color = FALSE)+
  ggtitle("B")+
  scale_shape(name="Metabolic Rate", labels=c("Respiration"=expression(R[dark]), "Photosynthesis"=expression(P[net])))

rfa

combined_plot <- ggarrange(rft, rfa, nrow = 2, ncol = 1)
combined_plot


##### GROSS P

# FINAL TREATMENTS

rft_2$GrossP <- (rft_2$p_correct - rft_2$r_correct)

rft_4 <- rft_2 %>% pivot_longer(cols=c('GrossP', 'Respiration'),
                                names_to='light_treatment',
                                values_to='rate')

rft_sum4<- rft_4 %>%
  group_by(population, light_treatment, treatment_temp) %>%
  summarize(mean=mean(rate), sd=sd(rate), n=n(), se=sd/sqrt(n))

rft_gross <- ggplot()+
  geom_point(aes(x=treatment_temp, y=mean, color=population, shape=light_treatment), rft_sum4, size=3)+
  geom_point(aes(x=treatment_temp, y=rate, color=population, shape=light_treatment), rft_4, alpha=0.2, size=3)+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), rft_sum4, width=0.2)+
  theme_classic()+
  theme(legend.text = element_text(size = 13))+
  ylab(expression(ΔDO~(umol~O[2]~g~coral^{"-1"}~mL~H[2]~O^{"-1"}~hr^{"-1"}))) +
  xlab("Temperature Treatment (ºC)")+
  geom_hline(yintercept = 0, alpha=0.2, linetype="dashed")+
  facet_grid(cols=vars(population))+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  guides(color = FALSE)+
  ggtitle("A - Treatment")+
  scale_shape(name="Metabolic Rate", labels=c("Respiration"=expression(R[dark]), "GrossP"=expression(P[gross])))+
  ylim(-1.9,2.4)

rft_gross

# FINAL AMBIENT

rfa_2$GrossP <- (rfa_2$p_correct - rfa_2$r_correct)

rfa_4 <- rfa_2 %>% pivot_longer(cols=c('GrossP', 'Respiration'),
                                names_to='light_treatment',
                                values_to='rate')

rfa_sum4<- rfa_4 %>%
  group_by(population, light_treatment, treatment_temp) %>%
  summarize(mean=mean(rate), sd=sd(rate), n=n(), se=sd/sqrt(n))

rfa_gross <- ggplot()+
  geom_point(aes(x=treatment_temp, y=mean, color=population, shape=light_treatment), rfa_sum4, size=3)+
  geom_point(aes(x=treatment_temp, y=rate, color=population, shape=light_treatment), rfa_4, alpha=0.2, size=3)+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), rfa_sum4, width=0.2)+
  theme_classic()+
  theme(legend.text = element_text(size = 13))+
  ylab(expression(ΔDO~(umol~O[2]~g~coral^{"-1"}~mL~H[2]~O^{"-1"}~hr^{"-1"}))) +
  xlab("Temperature Treatment (ºC)")+
  geom_hline(yintercept = 0, alpha=0.2, linetype="dashed")+
  facet_grid(cols=vars(population))+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  guides(color = FALSE)+
  ggtitle("B - Ambient")+
  scale_shape(name="Metabolic Rate", labels=c("Respiration"=expression(R[dark]), "GrossP"=expression(P[gross])))+
  ylim(-1.9,2.4)

rfa_gross

# Combine Plots
gross_plot <- ggarrange(rft_gross, rfa_gross, nrow = 2, ncol = 1)
gross_plot

#################### STATS ####################

anova_sym <- aov(sym_sci ~ treatment_temp*population, data = sym)
summary(anova_sym)
tukey_sym<-TukeyHSD(anova_sym)
tukey_sym


