


library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)


# edit filepath here
setwd("work/epfl-altshuler/Microcosm")




# Load MetaData ---------------------------------------------------------------


meta <- read.csv("DATA/gas_analysis/gas_flux_meta.csv", header = T)

# convert to time format
meta[c("start", "end")] <- lapply(meta[c("start", "end")], as.ITime, format = "%H:%M") 

# convert to US style (to match licor date)
meta$date <- ifelse(nchar(meta$date) != 10, paste0("0", meta$dat), meta$dat) # format date to XX/XX/XXXX
meta$Date <- paste(substr(meta$date, 4, 6), substr(meta$date, 1, 3), 2025, sep = "")

meta$Subsite <- paste(meta$treat, meta$bottle, meta$day, sep = "_")

meta <- meta[meta$Subsite != "fn_27_20", ] # remove error readings

# edit times more accurately for these subsites
meta[meta$Subsite == "a_17_83", ]$end <- meta[meta$Subsite == "a_17_83", ]$end - 30
meta[meta$Subsite == "f_3_2", ]$end <- meta[meta$Subsite == "f_3_2", ]$end - 30
meta[meta$Subsite == "f_3_42", ]$start <- meta[meta$Subsite == "f_3_42", ]$start + 30
meta[meta$Subsite == "a_15_97", ]$start <- meta[meta$Subsite == "a_15_97", ]$start + 30


# manual changes to the meta dataframe excel
# a_20_27 # manually changed hour that was incorrect
# f_2_62 had to also manually change as hour was incorrect
# f*_8*_4 also wrong hour
# f_11_42 changed end by 10 min (probably typo?)
# a_15_49 formatted :
# a_22_27 added 3 min to end
# f_8*_27 added 3 min to end
# f_3_0 changed hour
# fn_25_2 changed end by 20 minutes
# fn_26_35 changed end by 3 min
# f_5_35 changed by 10 min
# f_5_6 changed hour 
# f_2_0 changed hour


# CH4 lm -----------------------------------------------------------

# extract file names from analyser folders
CH4_files <- data.frame(CH4 = list.files("DATA/gas_analysis/gas_data/", pattern = "^Ultra100782", recursive = T))

CH4_files$Date <- substr(CH4_files$CH4, 1, 6)


# for each file format and get rep slopes
site_CH4_all <- data.frame()

for(i in 1:length(CH4_files$CH4)){

  # load individual files
  CH4 <- read.csv(paste("DATA/gas_analysis/gas_data/", CH4_files[i, ]$CH4, sep = ""), sep = ",") # select file here
  
  # format date and time
  CH4$Date <- substr(CH4$Time.Stamp, 1, 10) # format date
  CH4$Time <- as.ITime(substr(CH4$Time.Stamp, 12, 19))
  CH4$Seconds <- c(1:length(CH4$Time))
  
  # associated site data
  site_CH4 <- subset(meta, Date == CH4[i,]$Date) # matching date
  
  # time of site measurement
  start_CH4 <- as.ITime(min(na.omit(site_CH4$start)))
  end_CH4 <- as.ITime(max(na.omit(c(site_CH4$end)))) 
  
  # subset site gas data
  site_gas_CH4 <- CH4[CH4$Time >= start_CH4 & CH4$Time <= end_CH4,]
  
  # convert time to seconds from start of measurement
  site_gas_CH4$Seconds <- c(1:length(site_gas_CH4$Time.Stamp))
  
  # calculate slope for each subsite tech rep
  site_gas_CH4$Rep <- NA
  
  # get slope for each rep
  for(i in 1:length(site_CH4$start)){
    rep_start <- site_CH4[i, ]$start + 45 # adjust so that it is more accurate
    rep_end <- site_CH4[i, ]$end
    print(i)
    
    # assign subsite to gas data
    index <- which(site_gas_CH4$Time >= rep_start & site_gas_CH4$Time <= rep_end) # get index of first matching start time
    site_gas_CH4[index, ]$Rep <- site_CH4[i, ]$Subsite
    
    # calculate subsite slope
    model <- lm(CH4..ppm.~Seconds, data = site_gas_CH4[site_gas_CH4$Rep == site_CH4[i, ]$Subsite,])
    slope <- coef(summary(model))[2]  # coefficents, slope is second estimate value
    site_CH4[i, "slope_CH4"] <- slope
  }
  
  
  site_CH4_all <- rbind(site_CH4_all, site_CH4) # save results from one day to total dataframe

  plot_file <- ggplot(site_gas_CH4, aes(x = Seconds, y = CH4..ppm.)) + geom_point(alpha  = 0.1) + geom_smooth(data = na.omit(site_gas_CH4), aes(colour = Rep), method = "lm", fill = NA, fullrange = FALSE) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~Rep, scales = "free") + theme_bw()
  
  # save each plot
  ggsave(plot_file, file = paste0("ch4_plots/", gsub("/", "", site_CH4[1, ]$Date), ".png", sep = ""), width = 15, height = 8)

}








# CO2 N2O lm -------------------------------------------------------


# extract file names from analyser folders
CO2_files <- data.frame(CO2 = list.files("DATA/gas_analysis/gas_data/", pattern = "^Ultra100880", recursive = T))

CO2_files$Date <- substr(CO2_files$CO2, 1, 6)


# edit metadata
meta_co2 <- meta


# edit more accurate measurement times
meta_co2[meta_co2$Subsite == "an_28_2", ]$start <- meta_co2[meta_co2$Subsite == "an_28_2", ]$start + 90
meta_co2[meta_co2$Subsite == "a_14_35", ]$start <- meta_co2[meta_co2$Subsite == "a_14_35", ]$start + 45
meta_co2[meta_co2$Subsite == "a_16_83", ]$start <- meta_co2[meta_co2$Subsite == "a_16_83", ]$start + 15
meta_co2[meta_co2$Subsite == "f_2_83", ]$start <- meta_co2[meta_co2$Subsite == "f_2_83", ]$start + 30
meta_co2[meta_co2$Subsite == "f_7_83", ]$start <- meta_co2[meta_co2$Subsite == "f_7_83", ]$start + 15
meta_co2[meta_co2$Subsite == "f_11_97", ]$start <- meta_co2[meta_co2$Subsite == "f_11_97", ]$start + 30
meta_co2[meta_co2$Subsite == "f_2_97", ]$start <- meta_co2[meta_co2$Subsite == "f_2_97", ]$start + 30

# remove error measurement
meta_co2 <- meta_co2[meta_co2$Subsite != "fn_25_0", ]


# for each file format and get rep slopes
site_CO2_all <- data.frame()

for(i in 1:length(CO2_files$CO2)){
  
  # load individual files
  CO2 <- read.csv(paste("DATA/gas_analysis/", CO2_files[i, ]$CO2, sep = ""), sep = ",") # select file here
  
  # format date and time
  CO2$Date <- substr(CO2$Time.Stamp, 1, 10) # format date
  CO2$Time <- as.ITime(substr(CO2$Time.Stamp, 12, 19))
  CO2$Seconds <- c(1:length(CO2$Time))
  
  # associated site data
  site_CO2 <- subset(meta_co2, Date == CO2[i,]$Date) # matching date
  
  # time of site measurement
  start_CO2 <- as.ITime(min(na.omit(site_CO2$start)))
  end_CO2 <- as.ITime(max(na.omit(c(site_CO2$end)))) 
  
  # subset site gas data
  site_gas_CO2 <- CO2[CO2$Time >= start_CO2 & CO2$Time <= end_CO2,]
  
  # convert time to seconds from start of measurement
  site_gas_CO2$Seconds <- c(1:length(site_gas_CO2$Time.Stamp))
  
  # calculate slope for each subsite
  site_gas_CO2$Rep <- NA
  
  for(i in 1:length(site_CO2$start)){
    rep_start <- site_CO2[i, ]$start + 45 # adjust so that it is more accurate
    rep_end <- site_CO2[i, ]$end
    print(i)
    
    # assign subsite to gas data
    index <- which(site_gas_CO2$Time >= rep_start & site_gas_CO2$Time <= rep_end) # get index of first matching start time
    site_gas_CO2[index, ]$Rep <- site_CO2[i, ]$Subsite
    
    # calculate subsite slope
    model <- lm(CO2..ppm.~Seconds, data = site_gas_CO2[site_gas_CO2$Rep == site_CO2[i, ]$Subsite,])
    slope <- coef(summary(model))[2]  # coefficients, slope is second estimate value
    site_CO2[i, "slope_CO2"] <- slope
    
    # calculate subsite slope n20
    model2 <- lm(N2O..ppm.~Seconds, data = site_gas_CO2[site_gas_CO2$Rep == site_CO2[i, ]$Subsite,])
    slope2 <- coef(summary(model2))[2]  # coefficients, slope is second estimate value
    site_CO2[i, "slope_N2O"] <- slope2
    
  }
  
  
  site_CO2_all <- rbind(site_CO2_all, site_CO2) # save results from one day to total dataframe
  
  
  # CO2 plots
  plot_file <- ggplot(site_gas_CO2, aes(x = Seconds, y = CO2..ppm.)) + geom_point(alpha  = 0.1) + geom_smooth(data = na.omit(site_gas_CO2), aes(colour = Rep), method = "lm", fill = NA, fullrange = FALSE) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~Rep, scales = "free") + theme_bw()
  
  ggsave(plot_file, file = paste0("co2_plots/", gsub("/", "", site_CO2[1, ]$Date), ".png", sep = ""), width = 15, height = 8)
  
  
  # n2o plots
  plot_file2 <- ggplot(site_gas_CO2, aes(x = Seconds, y = N2O..ppm.)) + geom_point(alpha  = 0.1) + geom_smooth(data = na.omit(site_gas_CO2), aes(colour = Rep), method = "lm", fill = NA, fullrange = FALSE) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~Rep, scales = "free") + theme_bw()
  
  ggsave(plot_file2, file = paste0("n2o_plots/", gsub("/", "", site_CO2[1, ]$Date), ".png", sep = ""), width = 15, height = 8)
  
}





# Gas flux calculations ---------------------------------------------------


# load soil weights
soil_weight <- read.csv("DATA/gas_analysis/soil_weights.csv", header = T)

# merge with slope values
site_CO2_all[!(site_CO2_all$Subsite %in% site_CH4_all$Subsite),] # all co2 subsites are in ch4
gas_values <- left_join(site_CH4_all, site_CO2_all[c("Subsite", "slope_CO2", "slope_N2O")], by = "Subsite")

gas_values <- merge(gas_values, soil_weight[c("bottle", "soil_dry_weight")], by = "bottle")




# gas formula
# following gas flux calculations: https://info.gasmet.com/hubfs/Greenhouse%20Gas%20Flux%20Measurement%20Guide%20Gasmet%202020.pdf
gas_flux <- function(condition, air_temp, weight, slope){
  # weight in g
  
  # Volume
  bottle <- ifelse(grepl("f", condition), 
                   0.00023, # 230g for flooded soils
                   0.00025) # 250g for aerobic soils
                   
  # analyser volume 60cm3 in m3
  analyser <- 60/ 1000000 
  
  # displacement by soil volume (m3) assuming 1g = 1cm3
  soil_vol <- weight /1000000
  
  v <- bottle + analyser - soil_vol # total volume (m3)
  
  va <- 0.0224 # molar volume for ideal gas at 273 = 22.4 x 10-3 m 3 /mol
  T1 <- air_temp + 273.15 # air temperature + K (K)
  T2 <- 273.15 # standard temperature (K)
  
  # soil flux <- slope * (v / (va * (T2/T1))) * (1 / A) 
  flux <- slope * (v / (va * (T2/T1))) * (1 / weight) # umol g-1 s-1
  
  return(flux)
}



# calculate fluxes
gas_values$CH4_flux <- gas_flux(gas_values$treat, gas_values$tempC, gas_values$soil_dry_weight, gas_values$slope_CH4)
gas_values$CO2_flux <- gas_flux(gas_values$treat, gas_values$tempC, gas_values$soil_dry_weight, gas_values$slope_CO2)
gas_values$N2O_flux <- gas_flux(gas_values$treat, gas_values$tempC, gas_values$soil_dry_weight, gas_values$slope_N2O)

# save
write.csv(gas_values[, c(1,3,4, 11:18)], file = "DATA/gas_analysis/gas_data/nccr_gas_fluxes.csv", row.names = F)




# Plot treatment fluxes ---------------------------------------------------


gas_values <- read.csv("DATA/gas_analysis/gas_data/nccr_gas_fluxes.csv", header = T)

# CH4
# outlier f_5_42 is removed from anaylsis for now

# note: 16 bottle f_16_35 is meant to be aerobic, relabelled f_9_25
gas_values[gas_values$Subsite == "f_16_35", ]$bottle <- "9"
gas_values[gas_values$Subsite == "f_16_35", ]$Subsite <- "f_9_35"

# note: f* is 8* relabel 
gas_values[gas_values$bottle == "8*", ]$treat <- "f*"
gas_values[gas_values$bottle == "8*", ]$Subsite <- paste(gas_values[gas_values$bottle == "8*", ]$treat,
                                                         gas_values[gas_values$bottle == "8*", ]$bottle, 
                                                         gas_values[gas_values$bottle == "8*", ]$day,
                                                         sep = "_")

gas_values <- gas_values[gas_values$bottle != "8*", ] # remove this bottle from analysis

gas_values$treat <- factor(gas_values$treat, levels = c("a", "f", "an", "fn", "f*"))
ggplot(gas_values[gas_values$Subsite != "f_5_42", ], aes(x = day, y = CH4_flux, colour = bottle)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = bottle), alpha = 0.5, linewidth = 0.8) +
  facet_wrap(~treat, scales = "free", ncol = 2) + ylab("CH4 umol g-1 s-1") + xlab("Time (Days)")
#ggsave("ch4_all.png", plot = last_plot(), height = 4, width = 6)

#check per bottle
ggplot(gas_values, aes(x = day, y = CH4_flux, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = bottle), alpha = 0.5, linewidth = 1) +
  facet_wrap(~bottle, scales = "free_y", ncol = 4) + ylab("CH4 umol g-1 s-1") + xlab("Time (Days)")


ch4_means <- gas_values[gas_values$Subsite != "f_5_42", ] %>% group_by(treat, day) %>% summarise(mean_ch4 = mean(CH4_flux), sd_ch4 = sd(CH4_flux))

ggplot(ch4_means, aes(x = day, y = mean_ch4, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_ch4 + sd_ch4, ymin = mean_ch4 - sd_ch4), width = 0.4) +
  facet_wrap(~treat, scales = "free", ncol = 2) + ylab("CH4 umol g-1 s-1") + xlab("Time (Days)")
#ggsave("ch4_means.png", plot = last_plot(), height = 4, width = 6)

ch4_plot <- ggplot(ch4_means[ch4_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_ch4, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_ch4 + sd_ch4, ymin = mean_ch4 - sd_ch4), width = 0.4) +
  ylab("CH4 umol g-1 s-1") + xlab("") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment")





# CO2
# the outlier before is fine here
ggplot(na.omit(gas_values), aes(x = day, y = CO2_flux, colour = bottle)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = bottle), alpha = 0.5, linewidth = 0.8) +
  facet_wrap(~treat, scales = "free", ncol = 2) + ylab("CO2 umol g-1 s-1") + xlab("Time (Days)")
#ggsave("co2_all.png", plot = last_plot(), height = 4, width = 6)

ggplot(na.omit(gas_values), aes(x = day, y = CO2_flux, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = bottle), alpha = 0.5, linewidth = 1) +
  facet_wrap(~bottle, scales = "free", ncol = 4) + ylab("CO2 umol g-1 s-1") + xlab("Time (Days)")


co2_means <- na.omit(gas_values) %>% group_by(treat, day) %>% summarise(mean_co2 = mean(CO2_flux), sd_co2 = sd(CO2_flux))

ggplot(co2_means, aes(x = day, y = mean_co2, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_co2 + sd_co2, ymin = mean_co2 - sd_co2), width = 0.4) +
  facet_wrap(~treat, scales = "free", ncol = 2) + ylab("CO2 umol g-1 s-1") + xlab("Time (Days)")
#ggsave("co2_means.png", plot = last_plot(), height = 4, width = 6)

co2_plot <- ggplot(co2_means[co2_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_co2, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_co2 + sd_co2, ymin = mean_co2 - sd_co2), width = 0.4) +
  ylab("CO2 umol g-1 s-1") + xlab("") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment")




# N2O
# the outlier before is fine here
ggplot(na.omit(gas_values), aes(x = day, y = N2O_flux, colour = bottle)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = bottle), alpha = 0.5, linewidth = 0.8) +
  facet_wrap(~treat, scales = "free", ncol = 2) + ylab("N2O umol g-1 s-1") + xlab("Time (Days)")
#ggsave("n2o_all.png", plot = last_plot(), height = 4, width = 6)

n2o_means <- na.omit(gas_values) %>% group_by(treat, day) %>% summarise(mean_n2o = mean(N2O_flux), sd_n2o = sd(N2O_flux))

ggplot(n2o_means, aes(x = day, y = mean_n2o, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_n2o + sd_n2o, ymin = mean_n2o - sd_n2o), width = 0.4) +
  facet_wrap(~treat, scales = "free", ncol = 2) + ylab("N2O umol g-1 s-1") + xlab("Time (Days)")
#ggsave("n2o_means.png", plot = last_plot(), height = 4, width = 6)

n2o_plot <- ggplot(n2o_means[n2o_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_n2o, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_n2o + sd_n2o, ymin = mean_n2o - sd_n2o), width = 0.4) +
  ylab("N2O umol g-1 s-1") + xlab("Time (Days)") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment")


# without stat tests
ggarrange(ch4_plot, co2_plot, n2o_plot, ncol = 1, common.legend = T, legend = "right")








# Statistical tests -------------------------------------------------------

library(rstatix)
# friedman (non normal ranked anova)
gas_values_c <- gas_values[gas_values$treat %in% c("f", "a"), ] # remove controls

# friedman.test(CO2_flux ~ treat | day, gas_values_c)
# failed due to uneven sampling design (sacrifices)
# friedman.test(CH4_flux ~ treat | day, gas_values)

# wilcox test each timepoint
days <- unique(gas_values$day)
days <- days[!(days %in% c(2,124))] # remove day 2 comparison

pval_all <- data.frame()
for(day in days){
  print(day)
  co2_p <- wilcox.test(na.omit(gas_values[gas_values$treat == "f" & gas_values$day == day, ]$CO2_flux),
              na.omit(gas_values[gas_values$treat == "a" & gas_values$day == day, ]$CO2_flux))
  
  n2o_p <- wilcox.test(na.omit(gas_values[gas_values$treat == "f" & gas_values$day == day, ]$N2O_flux),
                       na.omit(gas_values[gas_values$treat == "a" & gas_values$day == day, ]$N2O_flux))
  
  ch4_p <- wilcox.test(gas_values[gas_values$treat == "f" & gas_values$day == day & gas_values$Subsite != "f_5_42", ]$CH4_flux,
                       gas_values[gas_values$treat == "a" & gas_values$day == day & gas_values$Subsite != "f_5_42", ]$CH4_flux)
  
  pvals <- data.frame(day = day, CH4 = ch4_p$p.value, CO2 = co2_p$p.value, N2O = n2o_p$p.value)
  pval_all <- rbind(pval_all, pvals)
}


# adjust for multiple testing
pval_all$CH4_adj <- p.adjust(pval_all$CH4, method = "fdr")
pval_all$CO2_adj <- p.adjust(pval_all$CO2, method = "fdr")
pval_all$N2O_adj <- p.adjust(pval_all$N2O, method = "fdr")

# add symbol
labs.function = function(x){
  case_when(x > 0.051 ~ "",
            x <= 0.051 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")}


pval_all$CH4_sym <- labs.function(pval_all$CH4_adj)
pval_all$CO2_sym <- labs.function(pval_all$CO2_adj)
pval_all$N2O_sym <- labs.function(pval_all$N2O_adj)




# ch4 plot
ch4_max <- ch4_means %>% group_by(day) %>% summarise(max_CH4 = max(mean_ch4)) # get max values for y position
pval_all <- merge(pval_all, ch4_max, by = "day" ) # merge for y position

# try applying log scale
ch4_plot <- ggplot(ch4_means[ch4_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_ch4, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_ch4 + sd_ch4, ymin = mean_ch4 - sd_ch4), width = 0.4) +
  ylab("CH4 umol g-1 s-1") + xlab("") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment") +
  geom_text(data = pval_all, aes(x = day, y = pval_all$max_CH4 + 0.0000002, label = CH4_sym), inherit.aes = F, size = 3) # + scale_y_log10()


# set groups that need to be scaled by day
ch4_means$scale_ch4 <- ""
ch4_means[ch4_means$day > 60, ]$scale_ch4 <- "zoomed"

pval_all$scale_ch4 <- ""
pval_all[pval_all$day > 60, ]$scale_ch4 <- "zoomed"

ggplot(ch4_means[ch4_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_ch4, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_ch4 + sd_ch4, ymin = mean_ch4 - sd_ch4), width = 0.4) +
  ylab("CH4 umol g-1 s-1") + xlab("") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment") +
  geom_text(data = pval_all, aes(x = day, y = pval_all$max_CH4 + 0.00000002, label = CH4_sym), inherit.aes = F, size = 3) +
  facet_wrap(~scale_ch4, scales = "free") + theme(strip.text = element_blank())



# co2 plot
co2_max <- co2_means %>% group_by(day) %>% summarise(max_CO2 = max(mean_co2))
pval_all <- merge(pval_all, co2_max, by = "day" ) # merge for y position

co2_plot <- ggplot(co2_means[co2_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_co2, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_co2 + sd_co2, ymin = mean_co2 - sd_co2), width = 0.4) +
  ylab("CO2 umol g-1 s-1") + xlab("") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment") +
  geom_text(data = pval_all, aes(x = day, y = pval_all$max_CO2 + 0.00005, label = CO2_sym), inherit.aes = F, size = 3) # + scale_y_log10()




# n2o plot
n2o_max <- n2o_means %>% group_by(day) %>% summarise(max_N2O = max(mean_n2o))
pval_all <- merge(pval_all, n2o_max, by = "day" ) # merge for y position

n2o_plot <- ggplot(n2o_means[n2o_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_n2o, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_n2o + sd_n2o, ymin = mean_n2o - sd_n2o), width = 0.4) +
  ylab("N2O umol g-1 s-1") + xlab("Time (Days)") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment") +
  geom_text(data = pval_all, aes(x = day, y = pval_all$max_N2O + 0.0000001, label = N2O_sym), inherit.aes = F, size = 3) 


# set groups that need to be scaled by day
n2o_means$scale_n2o <- ""
n2o_means[n2o_means$day > 43, ]$scale_n2o <- "zoomed"

pval_all$scale_n2o <- ""
pval_all[pval_all$day > 43, ]$scale_n2o <- "zoomed"

n2o_plot <- ggplot(n2o_means[n2o_means$treat %in% c("a", "f"), ], aes(x = day, y = mean_n2o, colour = treat)) + geom_point(alpha = 0.5) + 
  theme_bw() + geom_line(aes(group = treat), alpha = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymax = mean_n2o + sd_n2o, ymin = mean_n2o - sd_n2o), width = 0.4) +
  ylab("N2O umol g-1 s-1") + xlab("Time (Days)") + scale_colour_discrete(labels = c("Aerobic", "Flooded"), name = "Treatment") +
  geom_text(data = pval_all, aes(x = day, y = pval_all$max_N2O + 0.000000001, label = N2O_sym), inherit.aes = F, size = 3) +
  facet_wrap(~scale_n2o, scales = "free") + theme(strip.text = element_blank())



# log10 abs transformation doesnt work when crossing 0
# facet wrap and zoom look weird, and main trends are better as non scaled

# plot together
# looks better without scaling
ggarrange(ch4_plot, co2_plot, n2o_plot, ncol = 1, common.legend = T, legend = "right")
#ggsave("flux_plot_all.png", plot = last_plot(), height = 8, width = 8)





















