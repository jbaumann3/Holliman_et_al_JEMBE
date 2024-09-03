# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# Symbiont Density

# Last Updated: 6/4/23

#################### GETTING STARTED ####################

# Set working directory
setwd("~/Desktop/JEMBE Paper Data Analysis/Data")

library(dplyr)

# Load in data
sym <- read.csv("sym_count.csv")

# Create new columns
sym$total_sym <- sym$total_sym_number*10
sym$sym_sci <- sym$total_sym / 1000000

# Make new dfs
sym_exp <- sym[c(25:108), c(1,2,7,17)]
# delete MAH (outlier) from sym
sym_exp <- sym_exp[-c(56),]
sym_T0 <- sym[c(1:24),c(1,2,7,17)]

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

ggplot(bleach_end, aes(x=treatment_temp, y=mean, color=population)) + 
  geom_point(size=3, position=pd)+
  #geom_point(data=sym_end, aes(x=treatment_temp, y=sym_sci, color=population), size=3, alpha=0.2, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  theme_classic(base_size = 12)+
  ylab("Symbiont Density (10^6 cells / cm^2)")+
  xlab("Treatment (°C)")+
  labs(color="Population")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))

anova_sym2 <- aov(sym_sci ~ treatment_temp*population, data = sym_end)
summary(anova_sym2)
tukey_sym2<-TukeyHSD(anova_sym2)
tukey_sym2
