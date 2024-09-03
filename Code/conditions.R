# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# Lab Conditions

# Last Updated: 6/4/23

#################### GETTING STARTED ####################

# Set working directory
setwd("~/Desktop/JEMBE Paper Data Analysis/Data")

### In-situ HOBO temps for Deva's Astrangia thermal acclimation study - summer 2022####

#load packages
library(tidyverse)
library(ggsci)
library(lubridate)

#load data
coraltemp<-read.csv('HOBO.csv')

head(coraltemp)
str(coraltemp)

coraltemp$treatment=as.factor(as.character(coraltemp$treatment))
coraltemp$position=as.factor(as.character(coraltemp$position))
coraltemp$datetime=as.factor(as_datetime(coraltemp$datetime))


#convert F to C
coraltemp$temp_c = (coraltemp$temp_f-32)*(5/9)

#parse date and time
coraltemp$datetime2= mdy_hm(coraltemp$datetime)

coraltemp$Date = format(coraltemp$datetime2, "%d/%m/%Y")
coraltemp$Time = format(coraltemp$datetime2, "%H:%M")

head(coraltemp)

#preliminary plots
ggplot(coraltemp, aes(x=datetime2, y=temp_c, color=treatment, linetype=position))+
  geom_line()+
  scale_color_aaas()+
  theme_classic()

#trim data to start date (6/24)
#trim to end date
# 32 ended on 8/4
# 28 ended 8/6
# 22 ended 8/7
# 18 ended on 8/13

coraltemp2<-coraltemp[coraltemp$datetime2 > "2022-06-25" & coraltemp$datetime2 < "2022-08-04",]

coraltemp2$treatment=as.factor(as.character(coraltemp2$treatment))
coraltemp2$position=as.factor(as.character(coraltemp2$position))

ggplot(coraltemp2, aes(x=datetime2, y=temp_c, color=treatment, linetype=position))+
  geom_line(size=1)+
  scale_color_aaas()+
  theme_classic()+
  facet_grid(position~treatment)


#####summary stats#####

#overall means#
meantemp<- coraltemp2 %>%
  group_by(treatment, position) %>%
  na.omit() %>%
  summarize(mean=mean(temp_c), sd=sd(temp_c), n=n(), se=sd/sqrt(n))
meantemp

pd=position_dodge(width=0.1)

ggplot(meantemp, aes(x=treatment, y=mean, color=treatment, shape=position,group=position))+
  geom_point(size=3.5, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic()+
  scale_color_aaas()


#weekly means#
wktemp<- coraltemp2 %>%
  group_by(treatment, position, week = week(datetime2)) %>%
  na.omit() %>%
  summarize(mean=mean(temp_c), sd=sd(temp_c), n=n(), se=sd/sqrt(n))
wktemp

ggplot(wktemp, aes(x=week, y=mean, color=treatment, shape=position))+
  geom_point(size=3.5, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic()+
  scale_color_aaas()+
  geom_hline(yintercept = 18, linetype=2)+
  geom_hline(yintercept = 22, linetype=2)+
  geom_hline(yintercept = 28, linetype=2)+
  geom_hline(yintercept = 32, linetype=2)

#### Temperature from HOBO Loggers (Sumps & Cups)

# Read in data
hobo <- read.csv("HOBO.csv")

# Make "treatment" and "position" factors
hobo$treatment <- as.factor(hobo$treatment)
hobo$position <- as.factor(hobo$position)

# Make new column to convert temp_f to temp_c
# C = (F-32)/1.8
hobo$temp_c <- (hobo$temp_f - 32) / 1.8

# Delete light_lum_ft2 column 

hobo$light_lum_ft2 <- NULL 

# There are some rows with NA for temp--not sure why. Need to remove these. 
hobo <- na.omit(hobo)

# Make hobo$datetime a valid timeseries
hobo$datetime <- mdy_hm(hobo$datetime)

# Shorten data to experimental START (6/25/22) and STOP (8/4/22)
hobo2 <- hobo[hobo$datetime > "2022-06-25" & hobo$datetime < "2022-08-04",]

# Run some stats
hobo_stats <- hobo2 %>%
  group_by(treatment, position) %>%
  summarize(mean=mean(temp_c), sd=sd(temp_c), n=n(), se=(sd/sqrt(n)))

hobo2$position2 <- ifelse(grepl("cup", hobo2$position), "Cup", "Sump")

hobo2$treatment2 <- ifelse(grepl("18", hobo2$treatment), "18°C",
                        ifelse(grepl("22", hobo2$treatment), "22°C",
                               ifelse(grepl("28", hobo2$treatment), "28°C", "32°C")))

hobo2$treatment2 <- as.character(hobo2$treatment2)

hobo2$treatment2 <- factor(hobo2$treatment2, levels = "32°C","28°C","22°C", "18°C")

# Basic graph of temps in sumps and cups for duration of experiment: 
ggplot(hobo2, aes(x=datetime, y=temp_c, color=treatment2))+
  geom_line(size=0.4)+
  theme_classic(base_size = 12)+
  geom_hline(yintercept = 18, linetype=2)+
  geom_hline(yintercept = 22, linetype=2)+
  geom_hline(yintercept = 28, linetype=2)+
  geom_hline(yintercept = 32, linetype=2)+
  facet_wrap(~position2, ncol=1)+
  ylab("Temperature (°C)")+
  xlab("Date")+
  labs(color = "Treatment")+
  scale_color_brewer(type = "div", palette = 5, direction=-1)

                     
# Graph of mean temperatures for duration of experiment (SE is very low! error bars are not visible)
pd=position_dodge(width=0.1)

ggplot(hobo_stats, aes(x=treatment, y=mean, color=treatment, shape=position, group=position))+
  geom_point(size=3.5, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic()+
  geom_hline(yintercept = 18, linetype=2)+
  geom_hline(yintercept = 22, linetype=2)+
  geom_hline(yintercept = 28, linetype=2)+
  geom_hline(yintercept = 32, linetype=2)+
  ylab("Mean Temp (°C)")+
  xlab("Treatment")

#### PAR (Grow Lights)

# Read in data
par <- read.csv("par.csv")

par <- par[c(1:9),]

keycol <- "date"
valuecol <- "position"
gathercols <- c("c1_r1", "c1_r2","c1_r3","c1_r4","c2_r1","c2_r2","c2_r3","c2_r4","c3_r1","c3_r2","c3_r3","c3_r4","c4_r1","c4_r2","c4_r3","c4_r4", "c5_r1", "c5_r2","c5_r3","c5_r4","c6_r1","c6_r2","c6_r3","c6_r4","c7_r1","c7_r2","c7_r3","c7_r4","c8_r1","c8_r2","c8_r3","c8_r4")

par_long <- gather(par, keycol, valuecol, gathercols)
par_long$position <- par_long$keycol
par_long$par <- par_long$valuecol
par_long <- par_long[,-c(2:3)]

par2 <- read.csv("par_2.csv")
par2 <- par2[c(1:288),c(1:4)]

par_stats2 <- par2 %>%
  group_by(treatment)%>%
  summarize(mean=mean(par), sd=sd(par), n=n(), se=(sd/sqrt(n)))

ggplot(par_stats2, aes(x=treatment, y=mean))+
  geom_point(size=0.5, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic()

par2$date <- as.factor(par2$date)
par2$treatment <- as.factor(par2$treatment)
par2$position <- as.factor(par2$position)

anova_par <- aov(par ~ treatment*date, data = par2)
summary(anova_par)
tukey1<-TukeyHSD(anova_par)
tukey1

# Run some stats... 
par_stats <- par_long %>%
  group_by(position)%>%
  summarize(mean=mean(par), sd=sd(par), n=n(), se=(sd/sqrt(n)))

# Graph of PAR for each position per date 
ggplot(par_stats, aes(x=position, y=mean))+
  geom_point(size=0.5, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic()

par_stats2 <- par_long %>%
  group_by(date)%>%
  summarize(mean=mean(par), sd=sd(par), n=n(), se=(sd/sqrt(n)))

# Graph of PAR for each sump per date 
ggplot(par_stats2, aes(x=date, y=mean))+
  geom_point(size=3, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic(base_size = 12)+
  ylab("Mean PAR (umol/m2/s)")+
  xlab("Date")

par_stats3 <- par %>%
  group_by(sumps)%>%
  summarize(mean=mean(par), sd=sd(par), n=n(), se=(sd/sqrt(n)))

par_stats4 <- par %>%
  group_by(date)%>%
  summarize(mean=mean(par), sd=sd(par), n=n(), se=(sd/sqrt(n)))

mean(par$par)

ggplot(par_stats4, aes(x=date, y=mean))+
  geom_point(size=1, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic()+
  ylab("Mean PAR (umol/m2/s)")+
  xlab("Date")+
  geom_hline(yintercept = 132.375, linetype=2)

#### Flow Rates (Cups)

# Read in data 
conditions <- read.csv("conditions.csv")

# Make sump a character
conditions$sump <- as.character(conditions$sump)

# Make new column for mean flow 
conditions$meanflow <- (conditions$flow_l1 + conditions$flow_l2 + conditions$flow_r1 + conditions$flow_r2) / 4

# Need to change units of mean flow sec/100mL to mL/sec
conditions$meanflow2 <- conditions$meanflow / 100

# Run some stats... 
flow_stats <- conditions %>%
  group_by(date)%>%
  summarize(mean=mean(meanflow2), sd=sd(meanflow2), n=n(), se=(sd/sqrt(n)))

mean(flow_stats$mean)

flow_stats2 <- conditions %>%
  group_by(sump)%>%
  summarize(mean=mean(meanflow2), sd=sd(meanflow2), n=n(), se=(sd/sqrt(n)))

# Make a flow graph 
ggplot(flow_stats2, aes(x=sump, y=mean))+
  geom_point(size=3, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic(base_size = 12)+
  ylab("Mean Flow (ml/sec)")+
  xlab("Date")

#### Dissolved Oxygen (Sumps)

do_stats <- conditions %>%
  group_by(date)%>%
  summarize(mean=mean(do), sd=sd(do), n=n(), se=(sd/sqrt(n)))

mean(do_stats$mean)

# Make a flow graph 
ggplot(do_stats, aes(x=date, y=mean))+
  geom_point(size=3, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic(base_size = 12)+
  ylab("Mean % DO")+
  xlab("Date")

do_stats2 <- conditions %>%
  group_by(sump)%>%
  summarize(mean=mean(do), sd=sd(do), n=n(), se=(sd/sqrt(n)))

mean(do_stats2$mean)

# Make a flow graph 
ggplot(do_stats2, aes(x=sump, y=mean))+
  geom_point(size=3, position=pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1,position=pd)+
  theme_classic(base_size = 12)+
  ylab("Mean % DO")+
  xlab("Treatment")
