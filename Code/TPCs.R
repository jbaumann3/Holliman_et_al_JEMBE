# DEVA HOLLIMAN
# HONORS PROJECT 2022-2023
# DATA ANALYSIS
# TPCs

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
resp <- read.csv("respirometry.csv")
resp <- resp[,c(1:24)]

sym <- read.csv("sym_count.csv")
sym$total_sym <- sym$total_sym_number*10

sym$sym_sci <- sym$total_sym / 1000000
sym_exp <- sym[c(25:108), c(1,2,7,17)]
# delete MAH (outlier) from sym
sym_exp <- sym_exp[-c(56),]
sym_T0 <- sym[c(1:24),c(1,2,7,17)]

# Make some new columns
resp$population <- ifelse(grepl("NC", resp$coral_ID), "NC", 
                          ifelse(grepl("MA", resp$coral_ID), "MA", 
                                 ifelse(grepl("RI", resp$coral_ID), "RI", "blank")))

resp$coral_volume <- resp$empty_chamber_volume - resp$full_chamber_volume 
resp$coral_volume[is.na(resp$coral_volume)] <- 0

resp$mean_bw <- (resp$BW_1 + resp$BW_2 + resp$BW_3)/3

# Separate dfs by timepoint, delete repeated and unnecessary rows, and outlier MAH 
resp_final_treatments <- resp[c(102:161, 163:197, 199:216, 221:222, 227),]
resp_final_ambient <- resp[c(186:197, 204:223, 225:293),]

############## DATA WRANGLING FOR DF RESP_FINAL_TREATMENTS ##############

# Note that there are 33 blanks in this df
treatment_stats <- resp_final_treatments %>%
  group_by(population, treatment_temp)%>%
  summarize(mean=mean(treatment_temp), sd=sd(treatment_temp), n=n(), se=(sd/sqrt(n)))

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

resp_final_treatments$resp_rate_bw <- ifelse(grepl(0, resp_final_treatments$coral_volume), (resp_final_treatments$dark_DO_produced_umol * (1/resp_final_treatments$empty_chamber_volume) * (1/resp_final_treatments$dark_time_elapsed_hr)) ,
                                          (resp_final_treatments$dark_DO_produced_umol * (1/resp_final_treatments$mean_bw) * (1/resp_final_treatments$full_chamber_volume) * (1/resp_final_treatments$dark_time_elapsed_hr)))

resp_final_treatments <- merge(resp_final_treatments, sym_exp, by='coral_ID', all.x = TRUE, all.y=TRUE)
resp_final_treatments <- resp_final_treatments[,-c(31)]
resp_final_treatments$treatment_temp <- resp_final_treatments$treatment_temp.x
resp_final_treatments <- resp_final_treatments[-c(43,111),]

resp_final_treatments$resp_rate_sa <- ifelse(grepl(0, resp_final_treatments$coral_volume), (resp_final_treatments$dark_DO_produced_umol * (1/resp_final_treatments$empty_chamber_volume) * (1/resp_final_treatments$dark_time_elapsed_hr)) ,
                                             (resp_final_treatments$dark_DO_produced_umol * (1/resp_final_treatments$surface_area) * (1/resp_final_treatments$full_chamber_volume) * (1/resp_final_treatments$dark_time_elapsed_hr)))

# Fix some quick mistakes:
resp_final_treatments$light_time_start[is.na(resp_final_treatments$light_time_start)] <- 0
resp_final_treatments$light_time_stop[is.na(resp_final_treatments$light_time_stop)] <- 4000

# Calculate resp_final_treatments photosynthesis rate                                       
resp_final_treatments$light_time_elapsed_hr <- (resp_final_treatments$light_time_stop - resp_final_treatments$light_time_start) / 3600
resp_final_treatments$light_DO_produced_umol <- (resp_final_treatments$light_DO_final - resp_final_treatments$light_DO_initial) * (1000/31.999)

resp_final_treatments$photo_rate_bw <- ifelse(grepl(0, resp_final_treatments$coral_volume), (resp_final_treatments$light_DO_produced_umol * (1/resp_final_treatments$empty_chamber_volume) * (1/resp_final_treatments$light_time_elapsed_hr)) ,
                                           (resp_final_treatments$light_DO_produced_umol * (1/resp_final_treatments$mean_bw) * (1/resp_final_treatments$full_chamber_volume) * (1/resp_final_treatments$light_time_elapsed_hr)))

resp_final_treatments$photo_rate_sa <- ifelse(grepl(0, resp_final_treatments$coral_volume), (resp_final_treatments$light_DO_produced_umol * (1/resp_final_treatments$empty_chamber_volume) * (1/resp_final_treatments$light_time_elapsed_hr)) ,
                                              (resp_final_treatments$light_DO_produced_umol * (1/resp_final_treatments$surface_area) * (1/resp_final_treatments$full_chamber_volume) * (1/resp_final_treatments$light_time_elapsed_hr)))

# Standardize by subtracting resp_rate and photo_rate from blanks from each trial

#step 1: make df for corals and df for blanks
rft_blanks<-subset(resp_final_treatments, coral_ID=='blank')
head(rft_blanks)
ncol(rft_blanks)
rft_corals<-subset(resp_final_treatments, coral_ID!='blank')
head(rft_corals)

#step 2: Trim out the stuff we don't need for blanks
rft_blanks2<-rft_blanks[,c(4,30,34,37,38)]
head(rft_blanks2)

#step 2.5: rename columns
rft_blanks2<- rft_blanks2 %>%
  rename('blank_resp_rate_bw'='resp_rate_bw', 'blank_photo_rate_bw'='photo_rate_bw','blank_resp_rate_sa'='resp_rate_sa', 'blank_photo_rate_sa'='photo_rate_sa')
head(rft_blanks2)

#step 3: Mean blank values per trial
mean_blanks<- rft_blanks2 %>%
  group_by(trial) %>%
  summarize(mean_blank_R_bw= mean (blank_resp_rate_bw), mean_blank_P_bw=mean(blank_photo_rate_bw), mean_blank_R_sa= mean (blank_resp_rate_sa), mean_blank_P_sa=mean(blank_photo_rate_sa))
mean_blanks

#step 3: merge blank P and R to corals
rft_2<-merge(rft_corals, mean_blanks, by='trial')
head(rft_2)

#step 4: create corrected P/R columns
rft_2$p_correct_bw <- rft_2$photo_rate_bw - rft_2$mean_blank_P_bw
rft_2$r_correct_bw <- rft_2$resp_rate_bw - rft_2$mean_blank_R_bw
rft_2$p_correct_sa <- rft_2$photo_rate_sa - rft_2$mean_blank_P_sa
rft_2$r_correct_sa <- rft_2$resp_rate_sa - rft_2$mean_blank_R_sa

rft_2$r_correct_bw<-rft_2$r_correct_bw*(-1)

rft_2$r_correct_bw[rft_2$r_correct_bw < 0] <- 0

rft_2$p_gross <- rft_2$p_correct_bw + rft_2$r_correct_bw

# Run some stats
rfstats1_bw <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(r_correct_bw), sd=sd(r_correct_bw), n=n(), se=(sd/sqrt(n)))

rfstats1_bw$light_treatment <- "dark"
rfstats1_bw$sym_standard <- "no"

rfstats2_bw <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(p_correct_bw), sd=sd(p_correct_bw), n=n(), se=(sd/sqrt(n)))

rfstats2_bw$light_treatment <- "light"
rfstats2_bw$sym_standard <- "no"

rfstats1_sa <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(r_correct_sa), sd=sd(r_correct_sa), n=n(), se=(sd/sqrt(n)))

rfstats1_sa$light_treatment <- "dark"
rfstats1_sa$sym_standard <- "no"

rfstats2_sa <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(p_correct_sa), sd=sd(p_correct_sa), n=n(), se=(sd/sqrt(n)))

rfstats2_sa$light_treatment <- "light"
rfstats2_sa$sym_standard <- "no"

# Now standarize rates of respiration by symbiont density (NOTE: missing RIT and MAH)
# Note that this is P per 1,000,000 cells!
rft_2$p_correct_bw_sym <- (rft_2$p_correct_bw / rft_2$sym_sci)
rft_2$p_correct_sa_sym <- (rft_2$p_correct_sa / rft_2$sym_sci)

# Run some stats
rfstats2_bw_sym <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(p_correct_bw_sym), sd=sd(p_correct_bw_sym), n=n(), se=(sd/sqrt(n)))

rfstats2_bw_sym$light_treatment <- "light"
rfstats2_bw_sym$sym_standard <- "yes"

rfstats2_sa_sym <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(p_correct_sa_sym), sd=sd(p_correct_sa_sym), n=n(), se=(sd/sqrt(n)))

rfstats2_sa_sym$light_treatment <- "light"
rfstats2_sa_sym$sym_standard <- "yes"

# Combine tables rfstats_bw and rfstats_sa
rfstats_bw <- rbind(rfstats1_bw, rfstats2_bw, rfstats2_bw_sym)
rfstats_sa <- rbind(rfstats1_sa, rfstats2_sa, rfstats2_sa_sym)

# Make treatment_temp a chr
rfstats_bw$treatment_temp <- as.character(rfstats_bw$treatment_temp)
rfstats_sa$treatment_temp <- as.character(rfstats_sa$treatment_temp)

###################### NO TPCS JUST GRAPH ######################

pd=position_dodge(width=0.5)

rft_2$treatment_temp_int<-rft_2$treatment_temp.x

rft_2$population <- factor(rft_2$population, levels = c("MA", "RI", "NC"))

## graph gross photosynthesis

rft_gross <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(p_gross), sd=sd(p_gross), n=n(), se=(sd/sqrt(n)))

rft_gross$treatment_temp<-as.factor(rft_gross$treatment_temp)

## graph gross photosynthesis
ggplot() +
  geom_point(aes(x=treatment_temp, y=mean, color=population), rft_gross, size = 3, position = pd) +
  theme_classic(base_size = 11)+
  labs(x = 'Treatment (ºC)',
       y = 'Gross Photosynthesis (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  #geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(rft_gross, mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2, position = pd)+
  labs(color="Population", fill="Population")+
  ggtitle("Treatment Temperature")+
  scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))+
  ylim(0,1.5)

## another version of resp

rft_resp <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(r_correct_bw), sd=sd(r_correct_bw), n=n(), se=(sd/sqrt(n)))

rft_resp$treatment_temp<-as.factor(rft_resp$treatment_temp)

rft_2$treatment_temp<-as.factor(rft_2$treatment_temp)

ggplot() +
  geom_point(aes(x=treatment_temp, y=mean, color=population), rft_resp, size = 3, position = pd) +
  #geom_point(aes(x=treatment_temp, y=r_correct_bw, color=population), rft_2, size=3, position=pd, alpha=0.2)+
  theme_classic(base_size = 11)+
  labs(x = 'Treatment (ºC)',
       y = 'Gross Respiration (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  #geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(rft_resp, mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2, position = pd)+
  labs(color="Population", fill="Population")+
  ggtitle("Treatment Temperature")+
  scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))+
  ylim(0,1.5)

#### P:R 

rft_2$pr <- (rft_2$p_gross / (-1*(rft_2$r_correct_bw)))

rft_2_ish <- rft_2[-c(71),]

rft_pr <- rft_2_ish %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(pr), sd=sd(pr), n=n(), se=(sd/sqrt(n)))

ggplot() +
  geom_point(aes(x=treatment_temp, y=mean, color=population), rft_pr, size = 3, position = pd) +
  #geom_point(aes(x=treatment_temp, y=pr, color=population), rft_2_ish, size = 3, alpha=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(x = 'Treatment (ºC)',
       y = 'P:R')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(rft_pr, mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2, position = pd)+
  labs(color="Population", fill="Population")+
  ggtitle("Treatment Temperature")+
  scale_x_discrete(labels=c("18"="18.7","22" = "22.4","28" = "28.0", "32" = "31.5"))

view(rft_2[,c(34)])

deva <- rft_2[,c(1:3,25,44,47,52)]

####### THERMAL PERFORMANCE CURVES #######

# 1. Make two new dfs
    #   rft_tpc (data per coral) --> this is rft_2 in Justin's code
    #   mean_rft_tpc (means and error) --> this is EITHER error_data OR respdf

rft_tpc <- rft_2[,c(2,5,25,34,43:49)]

rft_tpc$treatment_temp <- rft_tpc$treatment_temp.x
rft_tpc <- rft_tpc[,-c(2,4)]


############## Some ANOVAS #######

rft_tpc$treatment_temp <- as.factor(rft_tpc$treatment_temp)
rft_tpc$treatment_temp <- as.factor(rft_tpc$treatment_temp)
rft_tpc$population <- as.factor(rft_tpc$population)
rft_tpc$population <- as.factor(rft_tpc$population)

anova_resp <- aov(r_correct_bw ~ treatment_temp*population, data = rft_tpc)
summary(anova_resp)
tukey_resp<-TukeyHSD(anova_resp)
tukey_resp

anova_photo <- aov(p_correct_bw ~ treatment_temp*population, data = rft_tpc)
summary(anova_photo)
tukey_photo<-TukeyHSD(anova_photo)
tukey_photo

anova_photo_gross <- aov(p_gross ~ treatment_temp*population, data = rft_tpc)
summary(anova_photo_gross)
tukey_photo_gross<-TukeyHSD(anova_photo_gross)
tukey_photo_gross

###################### EVERYTHING ABOVE THIS LINE IS PERFECT

### Net P and Gross R (bw)

# 2. Group data 

rft_mean_bw<-rft_tpc %>%
  group_by(population, treatment_temp) %>%
  summarize(meanp=mean(p_correct_bw), sdp=sd(p_correct_bw), n=n(), sep=sdp/sqrt(n),meanr=mean(r_correct_bw), sdr=sd(r_correct_bw), n=n(), ser=sdr/sqrt(n))

# 3. Model Selection

mean_rft_tpc$population <- as.factor(mean_rft_tpc$population)
mean_rft_tpc$measurement <- ifelse(grepl("dark", mean_rft_tpc$light_treatment), "Respiration", "Photosynthesis")
mean_rft_tpc$measurement <- as.factor(mean_rft_tpc$measurement)

mean_rft_tpc_bw <- mean_rft_tpc[c(1:24),]
mean_rft_tpc_bw$temp <- mean_rft_tpc_bw$treatment_temp
mean_rft_tpc_bw$rate <- mean_rft_tpc_bw$mean
mean_rft_tpc_bw <- mean_rft_tpc_bw[,-c(2:3)]

mean_rft_tpc_bw$temp <- as.integer(mean_rft_tpc_bw$temp)

# Need to put data in highly specific format
mean_df <- mean_rft_tpc_bw[,c(1,8:10)]
mean_df$incubation <- 'treatment'

mean_df_resp<-subset(mean_df, measurement=='Respiration')
mean_df_photo<-subset(mean_df, measurement=='Photosynthesis')


### HERE
mean_df_resp$rate <- -1*(mean_df_resp$rate)

mean_df <- rbind(mean_df_resp, mean_df_photo)

bw_rft_mean_fits <- mean_df %>%
  group_by(population, measurement, incubation) %>%
  nest(data = c(temp, rate)) %>%
  mutate(thomas1 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a,b,c,topt),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))

# 4. Model Predictions

# stack models
d_stack <- select(bw_rft_mean_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', thomas1)
d_stack

# get predictions using augment
newdata <- tibble(temp = seq(min(mean_df$temp), max(mean_df$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

d_preds

# 5. Make an initial plot, check to make sure everything looks good: 

ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(color = population)) +
  geom_point(aes(temp, rate, color=population), mean_df) +
  theme_classic(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_aaas()+
  facet_grid(model_name~measurement*incubation, scales='free')

# 6. Calculate estimated parameters

params <- d_stack %>%
  mutate(., params = map(fit, calc_params)) %>%
  select(-fit) %>%
  unnest(params)

params

# 7. Bootstrapping

#extract coefficients and refitting models using minpack.lm::nlsLM()

#get coefs
fitsR<-mutate(bw_rft_mean_fits, coefs=map(thomas1, coef))
head(fitsR)

# fit with nlsLM instead
fitsR2boots <- mutate(fitsR, nls_fit = map2(data, coefs, ~nlsLM(rate ~ thomas_2012(temp = temp, a, b, c, topt),
                                                                data = .x,
                                                                start = .y,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'))))

head(fitsR2boots)

# create empty list column
fitsR2boots <- mutate(fitsR2boots, bootstrap = list(rep(NA, n())))

# run for loop to bootstrap each refitted model
for(i in 1:nrow(fitsR2boots)){
  temp_data <- fitsR2boots$data[[i]]
  temp_fit <- nlsLM(rate ~ thomas_2012(temp = temp, a, b, c, topt),
                    data = temp_data,
                    start = fitsR2boots$coefs[[i]],
                    lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'thomas_2012'),
                    upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'thomas_2012'))
  
  boot <- Boot(temp_fit, method = 'residual')
  fitsR2boots$bootstrap[[i]] <- boot
  rm(list = c('temp_fit', 'temp_data', 'boot'))
}

fitsR2boots

# get the raw values of each bootstrap
fitsR2boots <- mutate(fitsR2boots, output_boot = map(bootstrap, function(x) x$t))


# calculate predictions 
fitsR2boots <- mutate(fitsR2boots, preds = map2(output_boot, data, function(x, y){
  temp <- as.data.frame(x) %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(y$temp), max(y$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = thomas_2012(temp = temp, a, b, c, topt))
  return(temp)
}))

fitsR2boots

#select and calc 95% CIs of predictions
boot_conf_preds <- select(fitsR2boots, population, preds) %>%
  unnest(preds) %>%
  group_by(population, incubation, measurement, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

boot_conf_preds

# 8. Wrangle the error bars 

bw_error_data <- mean_rft_tpc[c(1:24),]

bw_error_resp<-subset(bw_error_data, measurement=='Respiration')
bw_error_photo<-subset(bw_error_data, measurement=='Photosynthesis')

bw_error_resp$mean <- -1*(bw_error_resp$mean)

error <- rbind(bw_error_resp, bw_error_photo)

error$population <- factor(error$population, levels = c("MA", "RI", "NC"))
mean_df$population <- factor(mean_df$population, levels = c("MA", "RI", "NC"))
boot_conf_preds$population <- factor(boot_conf_preds$population, levels = c("MA", "RI", "NC"))
d_preds$population <- factor(d_preds$population, levels = c("MA", "RI", "NC"))

error$treatment_temp <- as.numeric(error$treatment_temp)

# 9. Make some pretty graphs!

pd=position_dodge(width=0.3)

# Units are Gross P/R are umol DO/g coral/L water/hr

# Graph with P/R side by side
ggplot() +
  geom_line(aes(temp, .fitted, color=population), d_preds) +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill=population), boot_conf_preds, alpha = 0.3) +
  geom_point(aes(temp, rate, color=population), mean_df, size = 2) +
  theme_classic(base_size = 12)+
  labs(x = 'Treatment (ºC)',
       y = 'Metabolic Rate')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"))+
  facet_wrap(~measurement)+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(error, mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2)+
  labs(color="Population", fill="Population")

#d_preds$temp <- as.factor(d_preds$temp)
#boot_conf_preds$temp <- as.factor(boot_conf_preds$temp)
#mean_df$temp <- as.factor(mean_df$temp)

# Graph with one metabolic rate (in this case, respiration)
ggplot() +
  geom_line(aes(temp, .fitted, color=population), subset(d_preds, measurement %in% "Respiration")) +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill=population), subset(boot_conf_preds, measurement %in% "Respiration"), alpha = 0.3) +
  geom_point(aes(temp, rate, color=population), subset(mean_df, measurement %in% "Respiration"), size = 3) +
  theme_classic(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Gross Respiration (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(subset(error, measurement %in% "Respiration"), mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2)+
  labs(color="Population", fill="Population")

## graph just photosynthesis
ggplot() +
  geom_line(aes(temp, .fitted, color=population), subset(d_preds, measurement %in% "Photosynthesis")) +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill=population), subset(boot_conf_preds, measurement %in% "Photosynthesis"), alpha = 0.3) +
  geom_point(aes(temp, rate, color=population), subset(mean_df, measurement %in% "Photosynthesis"), size = 3) +
  theme_classic(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Net Photosynthesis (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"))+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(subset(error, measurement %in% "Photosynthesis"), mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2)+
  labs(color="Population", fill="Population")

###################### NO TPCS JUST GRAPH ######################

pd=position_dodge(width=1.2)

rft_2$treatment_temp_int<-rft_2$treatment_temp.x

rft_2$r_correct_bw_neg <- rft_2$r_correct_bw * -1

rft_2$population <- factor(rft_2$population, levels = c("MA", "RI", "NC"))

# Graph with one metabolic rate (in this case, respiration)
ggplot() +
  geom_point(aes(temp, rate, color=population), subset(mean_df, measurement %in% "Respiration"), size = 3, position = pd) +
  geom_point(aes(x=treatment_temp_int, y=r_correct_bw_neg, color=population), rft_2, size = 3, alpha=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Gross Respiration (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(subset(error, measurement %in% "Respiration"), mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2, position = pd)+
  labs(color="Population", fill="Population")

## graph just photosynthesis
ggplot() +
  geom_point(aes(temp, rate, color=population), subset(mean_df, measurement %in% "Photosynthesis"), size = 3, position = pd) +
  geom_point(aes(x=treatment_temp_int, y=p_correct_bw, color=population), rft_2, size = 3, alpha=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Net Photosynthesis (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(subset(error, measurement %in% "Photosynthesis"), mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2, position = pd)+
  labs(color="Population", fill="Population")

### gross photo

rft_gross <- rft_2 %>%
  group_by(population, treatment_temp) %>%
  summarize(mean=mean(p_gross), sd=sd(p_gross), n=n(), se=(sd/sqrt(n)))

rft_gross$treatment_temp<-as.character(rft_gross$treatment_temp)
rft_gross$treatment_temp<-as.integer(rft_gross$treatment_temp)

## graph gross photosynthesis
ggplot() +
  geom_point(aes(x=treatment_temp, y=mean, color=population), rft_gross, size = 3, position = pd) +
  #geom_point(aes(x=treatment_temp_int, y=p_gross, color=population), rft_2, size = 3, alpha=0.2, position=pd)+
  theme_classic(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Gross Photosynthesis (umol DO/g/L/hr)')+
  scale_color_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  scale_fill_manual(values=c("purple4","cyan4", "goldenrod"), labels=c("18" = "18.7", "22" = "22.4","28" = "28.0", "32" = "31.5"))+
  geom_hline(yintercept=0, alpha=0.2, linetype='dashed')+
  geom_errorbar(rft_gross, mapping=aes(x=treatment_temp, ymin=mean-se, ymax=mean+se, color=population), width=0.2, position = pd)+
  labs(color="Population", fill="Population")+
  ggtitle("Treatment Temperatures")
