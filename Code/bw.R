# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# Buoyant Weight

# Last Updated: 6/4/23

#################### GETTING STARTED ####################

# Set working directory
setwd("~/Desktop/JEMBE Paper Data Analysis/Data")

# Read in data 
bw <- read.csv("bw.csv")
bw <- bw[c(1:236, 238:253),c(1:11)]

# Make some new columns
bw$population <- ifelse(grepl("NC", bw$coral_ID), "NC", 
                        ifelse(grepl("MA", bw$coral_ID), "MA", 
                               ifelse(grepl("RI", bw$coral_ID), "RI", "blank")))

bw$mean_bw <- (bw$BW_1 + bw$BW_2 + bw$BW_3)/3

# Run some stats
bw_stats <- bw %>%
  group_by(treatment_temp, population, timepoint) %>%
  summarize(mean=mean(mean_bw), sd=sd(mean_bw), n=n(), se=(sd/sqrt(n)))

# Transform data long to wide, so I can calculate percent change
bw_wide <- bw %>%
  pivot_wider(c(treatment_temp, coral_ID, population), names_from="timepoint", values_from="mean_bw")

# Delete row 57... RIU being finicky 
bw_wide <- bw_wide[-c(57),]

# Calculate percent change
bw_wide$percent_change <- (bw_wide$end - bw_wide$start) / bw_wide$start

bw_stats2 <- bw_wide %>%
  group_by(treatment_temp, population) %>%
  summarise(mean=mean(percent_change), sd=sd(percent_change), n=n(), se=(sd/sqrt(n)))

# Make treatment_temp a chr
bw_stats2$treatment_temp <- as.character(bw_stats2$treatment_temp)

# Graph bw percent change

pd=position_dodge(width=0.3)

ggplot(bw_stats2, aes(x=treatment_temp, y=mean, color=population))+
  geom_errorbar(aes(x=treatment_temp, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  geom_point(position=pd, size=3)+
  geom_line(aes(group=population), position=pd)+
  theme_classic()+
  ylab("% Change Buoyant Weight")+
  xlab("Treatment Temperature (°C)")+
  scale_color_manual(values=c("blue", "firebrick3","forestgreen"))

ggplot(bw_stats2, aes(x=treatment_temp, y=mean, color=population))+
  geom_bar()