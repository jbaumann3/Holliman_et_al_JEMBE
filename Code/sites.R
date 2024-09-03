# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# Collection Sites

# Last Updated: 6/4/23

### MAP

# Set working directory
setwd("~/Desktop/Honors Project/Deva_astrangia_thermal_phys_2022/deva_manuscript_jembe/Data")

# Load libraries
library (tidyverse)
library(lubridate)
library(ggmap)
library(maps)
#library(maptools)
library(scales)
library(mapproj)
library(sp)

# Load in data
gps <- read.csv("gps.csv")
gps <- gps[c(1:3)]

us <- map_data("world", "US", xlim=c(-77,-71),ylim=c(35,41))

pd<- position_dodge(0.03)

gps$site <- factor(gps$site, levels = c("MA", "RI", "NC"))

map<-ggplot()+
  geom_polygon(data=us, aes(x=long, y=lat, group= group), fill='gray70', color='black')+
  coord_fixed(xlim=c(-78,-68),ylim=c(34,44), ratio = 1)+
  theme_bw(base_size = 12)+
  geom_point(data=gps, aes(x=long, y=lat, color=site), position=pd, size=5)+
  ylab("Latitude")+
  xlab("Longitude")+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  labs(color="Site")+
  #theme(axis.text=element_text(size=15))+
  #theme(axis.title=element_text(size=15))+
  #theme(legend.title=element_text(size=15))+
  #theme(legend.text =element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
map

map2<-ggplot()+
  geom_polygon(data=us, aes(x=long, y=lat, group= group), fill='gray70', color='black')+
  coord_fixed(xlim = c(-90,-65), ylim = c(25,50), ratio = 1)+ 
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(35,43, 1))
map2

### MAKING A NOAA BUOY FIGURE 

# Set working directory
setwd("~/Desktop/JEMBE Paper Data Analysis/Data/NOAA Buoy Data")

# Read in data
MAbuoy <- read.csv("NOAA_MA_2021.csv")
MAbuoy$site <- "MA"
MAbuoy <- MAbuoy[-1,]

RIbuoy <- read.csv("NOAA_RI_2021.csv")
RIbuoy$site <- "RI"
RIbuoy <- RIbuoy[-1,]

NCbuoy <- read.csv("NOAA_NC_2021.csv")
NCbuoy$site <- "NC"
NCbuoy <- NCbuoy[-1,]

MAbuoy20 <- read.csv("NOAA_MA_2020.csv")
MAbuoy20$site <- "MA"
MAbuoy20 <- MAbuoy20[-1,]

RIbuoy20 <- read.csv("NOAA_RI_2020.csv")
RIbuoy20$site <- "RI"
RIbuoy20 <- RIbuoy20[-1,]

NCbuoy20 <- read.csv("NOAA_NC_2020.csv")
NCbuoy20$site <- "NC"
NCbuoy20 <- NCbuoy20[-1,]

# Merge dfs
noaa <- do.call("rbind", list(MAbuoy, RIbuoy, NCbuoy, MAbuoy20, RIbuoy20, NCbuoy20))

# Looks like there's some error with the RI dataset 
noaa <- noaa[!noaa$WTMP>50,]

# Make a column that is a valid timeseries
noaa$time <- paste(noaa$X.YY,"-",noaa$MM,"-",noaa$DD,"-",noaa$hh,"-",noaa$mm)
noaa$time <- ymd_hm(noaa$time)

noaa$MM <- as.integer(noaa$MM)

noaa$MM <- month.name[noaa$MM]

# The original dataset is way too big to graph, so run summary stats
noaa$WTMP <- as.numeric(noaa$WTMP)
#noaa$MM <- as.integer(noaa$MM)

noaa_stats <- noaa %>%
  group_by(MM, site) %>%
  dplyr::summarize(mean=mean(WTMP), sd=sd(WTMP), n=n(), se=(sd/sqrt(n)))

noaa_stats$MM <- factor(noaa_stats$MM, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
noaa_stats$site <- factor(noaa_stats$site, levels = c("MA", "RI", "NC"))

# Graph water temp over time 
pd=position_dodge(width=0.5)

ggplot(noaa_stats, aes(x=MM, y=mean, color=site))+
  geom_errorbar(aes(x=MM, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  geom_point(position=pd, size=3)+
  geom_line(aes(group=site), position=pd)+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  theme_classic(base_size = 12)+
  ylab("Mean SST (°C)")+
  xlab("")+
  labs(color="Site")
