# data analysis RQ4
# Gerard H. ros & Hongzheng
# 13 May 2024

# require packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)

# read in the database
d1 <- as.data.table(read_xlsx('data/RQ4_database.xlsx',sheet='RQ4_database'))

print(data)
data$kpi_mean<-as.numeric(data$kpi_mean)
# convert the variables to a factor with the desired order
data$end_product_type<- factor(data$end_product_type, levels = c("solid_fraction","liquid_fraction","microalgae","struvite","permeate","concentrate","ammonium_salts","ammonium_sulphate","ammonium_nitrate"))
data$experiment_scale<- factor(data$experiment_scale, levels = c("full_scale","farm","pilot","batch","laboratory"))
data$animal_control<- factor(data$animal_control, levels = c("mix","pig","poultry","cattle","dairy"))
data$man_code<-factor(data$man_code,levels = c("NRECO","MTREAT"))
data$kpi_code<-factor(data$kpi_code,levels=c("TN","TAN","TP","NH3"))
data$man_treatment<-factor(data$man_treatment,levels = c("ammonia scrubbing","stripping and scrubbing","struvite precipitation","centrifugation","screw press","acidification","ultrafiltration","microfiltration","reverse osmosis","porous polypropylene","pressure filtration","microalgae","gas permeable","wetland","NDN"))
# what are the unique KPIs included
table(data$kpi_code)

summary(data)
#remove the missing values
#data<-na.omit(data$kpi_mean)
#print(data)
# convert dataframe to data.table
data <- as.data.table(data)
# subset the file for NRECO and MTREAT
d1 <- data[man_code =='NRECO']
# boxplot for the end products
ggplot(d1,aes(x= end_product_type,y=kpi_mean,fill=kpi_code)) + 
  geom_boxplot() + 
  theme_bw() + xlab('end product') + ylab('recovery (%)')+
  scale_x_discrete(limits = c("solid_fraction","liquid_fraction","microalgae","struvite","permeate","concentrate","ammonium_salts","ammonium_sulphate","ammonium_nitrate")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
# boxplot for the scales
ggplot(d1,aes(x= end_product_type,y=kpi_mean,fill=experiment_scale)) + 
  geom_boxplot() + 
  theme_bw() + xlab('end product') + ylab('recovery (%)')+
  scale_x_discrete(limits = c("solid_fraction","liquid_fraction","microalgae","struvite","permeate","concentrate","ammonium_salts","ammonium_sulphate","ammonium_nitrate")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())

## regression model for the recovery rate of TN
d2 <- d1[kpi_code =='TN']
#remove the missing values
#d2<-na.omit(d2$kpi_mean)

# Load necessary packages
library(dplyr)

# Check unique values for experiment_scale
cat("\nUnique values in experiment_scale before conversion:\n")
print(unique(d2$experiment_scale))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    experiment_scale = factor(experiment_scale, levels = c("full_scale", "farm", "pilot", "batch", "laboratory"))
  )

# Check the unique values after conversion
cat("\nUnique values in experiment_scale after conversion:\n")
print(unique(d2$experiment_scale))

# Check for levels that are not present in the data
cat("\nLevels in experiment_scale:\n")
print(levels(d2$experiment_scale))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ end_product_type * experiment_scale, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))

## Regression model for the recovery rate of TAN
d2 <- d1[kpi_code =='TAN']

# Check unique values for experiment_scale
cat("\nUnique values in experiment_scale before conversion:\n")
print(unique(d2$experiment_scale))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    experiment_scale = factor(experiment_scale, levels = c("full_scale", "farm", "pilot", "batch", "laboratory"))
  )

# Check the unique values after conversion
cat("\nUnique values in experiment_scale after conversion:\n")
print(unique(d2$experiment_scale))

# Check for levels that are not present in the data
cat("\nLevels in experiment_scale:\n")
print(levels(d2$experiment_scale))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ end_product_type * experiment_scale, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))

## Regression model for the recovery rate of TP
d2 <- d1[kpi_code =='TP']

# Check unique values for experiment_scale
cat("\nUnique values in experiment_scale before conversion:\n")
print(unique(d2$experiment_scale))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    experiment_scale = factor(experiment_scale, levels = c("full_scale", "farm", "pilot", "batch", "laboratory"))
  )

# Check the unique values after conversion
cat("\nUnique values in experiment_scale after conversion:\n")
print(unique(d2$experiment_scale))

# Check for levels that are not present in the data
cat("\nLevels in experiment_scale:\n")
print(levels(d2$experiment_scale))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ end_product_type * experiment_scale, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))







# subset the file for MTREAT
d1 <- data[man_code =='MTREAT']
# boxplot for the end products
ggplot(d1,aes(x= man_treatment,y=kpi_mean,fill=kpi_code)) + 
  geom_boxplot() + 
  theme_bw() + xlab('manure processing') + ylab('removal rate (%)')+
  scale_x_discrete(limits = c("ammonia scrubbing","stripping and scrubbing","struvite precipitation","centrifugation","screw press","acidification","ultrafiltration","microfiltration","reverse osmosis","porous polypropylene","pressure filtration","microalgae","gas permeable","wetland","NDN")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
# boxplot for the scales
ggplot(d1,aes(x= man_treatment,y=kpi_mean,fill=experiment_scale)) + 
  geom_boxplot() + 
  theme_bw() + xlab('manure processing') + ylab('removal rate (%)')+
  scale_x_discrete(limits = c("ammonia scrubbing","stripping and scrubbing","struvite precipitation","centrifugation","screw press","acidification","ultrafiltration","microfiltration","reverse osmosis","porous polypropylene","pressure filtration","microalgae","gas permeable","wetland","NDN")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())

library(dplyr)

# subset the file for removal rate of TN
d2 <- d1[kpi_code =='TN']

# Check unique values for animal control
cat("\nUnique values in man_treatment before conversion:\n")
print(unique(d2$man_treatment))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    animal_control = factor(man_treatment, levels = c("struvite precipitation","stripping and scrubbing", "centrifugation","ammonia scrubbing","screw press","ultrafiltration","microfiltration","reverse osmosis","microalgae"))
  )

# Check the unique values after conversion
cat("\nUnique values in man_treatment after conversion:\n")
print(unique(d2$man_treatment))

# Check for levels that are not present in the data
cat("\nLevels in man_treatment:\n")
print(levels(d2$man_treatment))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ experiment_scale * man_treatment, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))


## subset the file for removal rate of TAN
d2 <- d1[kpi_code =='TAN']

# Check unique values for animal control
cat("\nUnique values in man_treatment before conversion:\n")
print(unique(d2$man_treatment))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    animal_control = factor(man_treatment, levels = c("struvite precipitation","stripping and scrubbing", "centrifugation","ammonia scrubbing","screw press","ultrafiltration","microfiltration","reverse osmosis","microalgae"))
  )

# Check the unique values after conversion
cat("\nUnique values in man_treatment after conversion:\n")
print(unique(d2$man_treatment))

# Check for levels that are not present in the data
cat("\nLevels in man_treatment:\n")
print(levels(d2$man_treatment))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ experiment_scale * man_treatment, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))

## subset the file for removal rate of NH3
d2 <- d1[kpi_code =='NH3']

# Check unique values for animal control
cat("\nUnique values in man_treatment before conversion:\n")
print(unique(d2$man_treatment))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    animal_control = factor(man_treatment, levels = c("struvite precipitation","stripping and scrubbing", "centrifugation","ammonia scrubbing","screw press","ultrafiltration","microfiltration","reverse osmosis","microalgae"))
  )

# Check the unique values after conversion
cat("\nUnique values in man_treatment after conversion:\n")
print(unique(d2$man_treatment))

# Check for levels that are not present in the data
cat("\nLevels in man_treatment:\n")
print(levels(d2$man_treatment))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ experiment_scale * man_treatment, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))


## subset the file for removal rate of TP
d2 <- d1[kpi_code =='TP']

# Check unique values for animal control
cat("\nUnique values in man_treatment before conversion:\n")
print(unique(d2$man_treatment))

# Convert columns to factors with specified levels
d2 <- d2 %>%
  mutate(
    end_product_type = as.factor(end_product_type),
    animal_control = factor(man_treatment, levels = c("struvite precipitation","stripping and scrubbing", "centrifugation","ammonia scrubbing","screw press","ultrafiltration","microfiltration","reverse osmosis","microalgae"))
  )

# Check the unique values after conversion
cat("\nUnique values in man_treatment after conversion:\n")
print(unique(d2$man_treatment))

# Check for levels that are not present in the data
cat("\nLevels in man_treatment:\n")
print(levels(d2$man_treatment))

# Drop unused levels
d2 <- droplevels(d2)

# Verify the data after dropping unused levels
cat("\nSummary of d2 after dropping unused levels:\n")
print(summary(d2))

# Check if experiment_scale still has fewer than two levels
if (length(unique(d2$experiment_scale)) < 2) {
  stop("experiment_scale must have at least two levels")
}

# Fit the linear model again
m1 <- lm(kpi_mean ~ experiment_scale * man_treatment, data = d2)
cat("\nSummary of the linear model:\n")
print(summary(m1))
