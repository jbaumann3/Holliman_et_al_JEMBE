# DEVA HOLLIMAN
# ASTRANGIA PROJECT
# DATA ANALYSIS
# Fv/Fm

# Last Updated: 6/4/23

#################### GETTING STARTED ####################

# Set working directory
setwd("~/Desktop/JEMBE Paper Data Analysis/Data")

# Load libraries
library(dplyr)
library(tidyverse)
library(broom)
library(data.table)
library(emmeans)
library(patchwork)
library(ggeffects)
library(lme4)
library(knitr)
library(ggh4x) #to add colors to facet_grid labels

# Read in data 
pam <- read.csv("pam.csv")
pam <- pam[c(1:568), c(1:14)]

# Make some new columns
pam$meanfvfm <- (pam$Fv.Fm_1 + pam$Fv.Fm_2 + pam$Fv.Fm_3)/3

pam$population <- ifelse(grepl("NC", pam$coral_ID), "NC", 
                         ifelse(grepl("MA", pam$coral_ID), "MA", "RI"))



##########################################################
# # Run some stats 
# pam_stats <- pam %>%
#   group_by(week, population, treatment_temp) %>%
#   dplyr::summarize(mean=mean(meanfvfm), sd=sd(meanfvfm), n=n(), se=(sd/sqrt(n)))
# 
# pam_stats$population <- as.factor(pam_stats$population)
# pam_stats$treatment_temp <- as.factor(pam_stats$treatment_temp)
# 
# # Do some organizational things 
# pd=position_dodge(width=0.4)
# pam_stats$treatment_temp <- as.character(pam_stats$treatment_temp)
# pam_stats$week <- as.character(pam_stats$week)
# 
# # Rearrange levels
# pam_stats2 <- pam_stats                             
# pam_stats2$population <- factor(pam_stats2$population, levels = c("MA", "RI", "NC"))
# 
# pam_stats2$real_temp <- ifelse(grepl("18", pam_stats2$treatment_temp), "18",
#                                ifelse(grepl("22", pam_stats2$treatment_temp), "22",
#                                       ifelse(grepl("28", pam_stats2$treatment_temp), "28", "32")))
# 
# # Graph...
# ggplot(pam_stats2, aes(x=week, y=mean, color=real_temp))+
#   geom_point(data=pam_stats2, aes(x=week, y=mean, color=real_temp), size=4, alpha=0.2, position=pd)+
#   geom_point(size=4, position=pd)+
#   geom_errorbar(aes(x=week, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
#   #geom_line(aes(group=treatment_temp), position=pd)+
#   geom_smooth(aes(group=treatment_temp), method='lm', alpha=0.1)+
#   theme_classic(base_size = 13)+
#   ylab("Fv/Fm")+
#   xlab("Week")+
#   facet_wrap(population ~ .)+
#   scale_color_manual(values=c("purple4","cyan4", "springgreen4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
#   labs(color = "Treatment (°C)")
# 
# pam_stats2$treatment_temp2 <- ifelse(grepl("18", pam_stats2$real_temp), "18.7", 
#                           ifelse(grepl("22", pam_stats2$real_temp), "22.4", 
#                                  ifelse(grepl("28", pam_stats2$real_temp), "28.0", "31.5")))
# 
# pam_stats2$treatment_temp2 <- factor(pam_stats2$treatment_temp2, levels = c("18.7", "22.4", "28.0", "31.5"))
# 
# 
# ggplot(pam_stats2, aes(x=week, y=mean, color=population))+
#   geom_point(data=pam_stats2, aes(x=week, y=mean, color=population), size=3, alpha=0.2, position=pd)+
#   geom_point(size=3, position=pd)+
#   geom_errorbar(aes(x=week, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
#   #geom_line(aes(group=population), position=pd)+
#   geom_smooth(aes(group=population), method='lm', alpha=0.1)+
#   theme_classic(base_size = 12)+
#   ylab("Fv/Fm")+
#   xlab("Week")+
#   facet_wrap(treatment_temp2 ~ .)+
#   scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
#   labs(color = "Population")
# 
# # Run Stats for Fv/Fm over time
# pam$treatment_temp=as.factor(pam$treatment_temp)
# pam$week=as.factor(pam$week)
# 
# anova_pam <- aov(meanfvfm ~ treatment_temp*population*week, data = pam)
# summary(anova_pam)
# tukey<-TukeyHSD(anova_pam)
# tukey
# 
# anova_pam1 <- aov(meanfvfm ~ treatment_temp*population, data = pam)
# summary(anova_pam1)
# tukey1<-TukeyHSD(anova_pam1)
# tukey1
# 
# anova_pam2 <- aov(meanfvfm ~ population*week, data = pam)
# summary(anova_pam2)
# tukey2<-TukeyHSD(anova_pam2)
# tukey2
# 
# anova_pam3 <- aov(meanfvfm ~ treatment_temp*week, data = pam)
# summary(anova_pam3)
# tukey3<-TukeyHSD(anova_pam3)
# tukey3
# 
# pam_stats3 <- pam %>%
#   group_by(week) %>%
#   dplyr::summarize(mean=mean(meanfvfm), sd=sd(meanfvfm), n=n(), se=(sd/sqrt(n)))
# 
# 
#   transform(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))
# 
# 
# ylab("Fv/Fm")+
#   xlab("Treatment (°C)")+
#   labs(color="Population")+
#   scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
#   scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))
# 
# anova_pam <- aov(mean ~ treatment_temp*population, data = pam_stats)
# summary(anova_pam)
# tukey1<-TukeyHSD(anova_pam)
# tukey1
# 
# pam_stats$week <- as.factor(pam_stats$week)
# 
# anova_pam_3 <- aov(mean ~ treatment_temp*population*week, data = pam_stats)
# summary(anova_pam_3)
# TukeyHSD(anova_pam_3)
# 
# 
# pam$week <- as.factor(pam$week)
# pam$treatment_temp <- as.factor(pam$treatment_temp)
# pam$population <- as.factor(pam$population)
# 
# anova_pam <- aov(meanfvfm ~ treatment_temp*population*week, data = pam)
# summary(anova_pam)
# tukey1<-TukeyHSD(anova_pam)
# tukey1
# 
# anova_pam_3 <- aov(meanfvfm ~ treatment_temp*population*week, data = pam)
# summary(anova_pam_3)
# TukeyHSD(anova_pam_3)
# 
# ###################### JUST LOOK AT WEEK 6
# 
# pam_end <- pam[c(421:488,490:505),]
# 
# end_stats <- pam_end %>%
#   group_by(population, treatment_temp) %>%
#   dplyr::summarize(mean=mean(meanfvfm), sd=sd(meanfvfm), n=n(), se=(sd/sqrt(n)))
# 
# end_stats$population <- factor(end_stats$population, levels = c("MA", "RI", "NC"))
# end_stats$treatment_temp <- as.factor(end_stats$treatment_temp)
# 
# ggplot(end_stats, aes(x=population, y=mean, color=treatment_temp)) + 
#   geom_point(size=4, position=pd)+
#   geom_point(data=pam_end, aes(x=population, y=meanfvfm, color=treatment_temp), size=4, alpha=0.2, position=pd)+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
#   theme_classic(base_size = 15)+
#   ylab("Fv/Fm")+
#   xlab("Population")+
#   labs(color="Treatment (°C)")+
#   scale_color_manual(labels=c("T0 Control" = "Control", "18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"), values=c("purple4","cyan4", "springgreen4", "goldenrod"))
# 
# pam_end$population <- factor(pam_end$population, levels = c("MA", "RI", "NC"))
# pam_end$treatment_temp <- as.factor(pam_end$treatment_temp)
# 
# ggplot(end_stats, aes(x=treatment_temp, y=mean, color=population)) + 
#   geom_point(size=3, position=pd)+
#   #geom_point(data=pam_end, aes(x=treatment_temp, y=meanfvfm, color=population), size=3, alpha=0.2, position=pd)+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
#   theme_classic(base_size = 12)+
#   ylab("Fv/Fm")+
#   xlab("Treatment (°C)")+
#   labs(color="Population")+
#   scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
#   scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))
# 
# pam_end_anova <- aov(meanfvfm ~ treatment_temp*population, data = pam_end)
# summary(pam_end_anova)
# TukeyHSD(pam_end_anova)

######################

#Linear Regression
pam$treatment_temp=as.factor(pam$treatment_temp)
pam$population <- factor(pam$population, levels = c("MA", "RI", "NC"))


lm1<- lm(meanfvfm ~ week*treatment_temp*population, data=pam)
summary(lm1) 
anova(lm1)



test1<- lm1 %>%
  augment(se_fit=TRUE,  interval = "confidence",
          conf.level = 0.95)
test1



lmgraph1<- lm1 %>%
  augment(se_fit=TRUE, interval = "confidence", conf.level = 0.95) %>%
  ggplot(aes(x=week, y=meanfvfm, color= treatment_temp))+
  #geom_ribbon(aes(ymin=.lower, ymax=.upper),color='grey',alpha=0.1)+
  geom_point()+
  geom_line(aes(y=.fitted), linewidth=1)+
  theme_classic()+
  theme_classic(base_size = 13)+
  ylab("Fv/Fm")+
  xlab("Week")+
  facet_grid(population ~ .)+
  scale_color_manual(values=c("#2166ac","#92c5de", "#f4a582", "#b2182b"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  labs(color = "Treatment (°C)")
lmgraph1





# same graph using ggpredict to get error ribbon

# Extract the prediction data frame
predlm<-ggpredict(lm1, terms = c('week','treatment_temp','population'))
predlm

predlm2<-predict_response(lm1, terms = c('week','treatment_temp','population'))
predlm2

print(predlm2, collapse_table = TRUE, collapse_ci = TRUE)

#pam$week<-as.numeric(pam$week)

strip<-strip_themed(text_y=element_text(color='white'),background_y= list(element_rect(fill=alpha("purple4",1.0)),
                                       element_rect(fill="cyan4"),
                                       element_rect(fill="goldenrod")))



newfvfmgraph<-ggplot(predlm2, aes(x = x, y = predicted, colour = group)) +
  geom_jitter(data=pam, aes(x=week, y=meanfvfm, color=treatment_temp), width=0.1,alpha=0.3)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high),alpha=0.05)+
  geom_line(linewidth=1.3) +
  facet_grid2(facet~., strip=strip)+
  theme_classic(base_size = 13)+
  ylab("Fv/Fm")+
  xlab("Week")+
  scale_color_manual(values=c("#2166ac","#92c5de", "#f4a582", "#b2182b"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  labs(color = "Treatment (°C)")+
  theme(text=element_text(size=16))

newfvfmgraph


  


#95% CI coef plot
coefs<-tidy(lm1, quick=FALSE)
coefs

ci<-data.table(confint(lm1), keep.rownames='term')
ci

cidf<-cbind(coefs,ci)
cidf

cidf<-cidf[,-6]

cidf<- cidf %>%
  rename("lower"="2.5 %",
         "upper"="97.5 %")

cidf
cidf$term=as.factor(cidf$term)

ggplot(data=cidf, aes(x=estimate, y=term))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_point(size=3)+
  geom_errorbarh(aes(xmax=lower, xmin=upper),height=0.2)+
  theme_classic()

#slope comparison -> from here: https://stats.stackexchange.com/questions/33013/what-test-can-i-use-to-compare-slopes-from-two-or-more-regression-models

lm1
anova(lm1)
lm1$coefficients

slopecomp<-lstrends(lm1, c("population","treatment_temp"), var="week")
slopecomp

pairs<-as.data.frame(pairs(slopecomp))

sigpairs<- pairs %>%
  filter(p.value <0.05)

sigpairs

#### Playing with color options 

devastest<-ggplot(pam_stats2, aes(x=week, y=mean, color=real_temp))+
  geom_point(data=pam_stats2, aes(x=week, y=mean, color=real_temp), size=4, alpha=0.2, position=pd)+
  geom_point(size=4, position=pd)+
  geom_errorbar(aes(x=week, ymin=mean-se, ymax=mean+se), width=0.2, position=pd)+
  #geom_line(aes(group=treatment_temp), position=pd)+
  geom_smooth(aes(group=treatment_temp), method='lm', alpha=0.1)+
  theme_classic(base_size = 13)+
  ylab("Fv/Fm")+
  xlab("Week")+
  facet_grid(population ~ .)+
  scale_color_manual(values=c("#2166ac","#92c5de", "#f4a582", "#b2182b"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  labs(color = "Treatment (°C)")
devastest

newfvfmgraph/devastest



