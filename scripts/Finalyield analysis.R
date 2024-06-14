# data analysis RQ1
# Maria Anna Antonovardaki & Gerard Ros
# 22 May 2024

#required packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)

# read in the database (NOTE that I updated some SE and SD values)
d1a <- as.data.table(read_xlsx('C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Finalyield_database.xlsx',sheet = 1))
d1b <- as.data.table(read_xlsx('C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Finalyield_database.xlsx',sheet = 2))
d1c <- as.data.table(read_xlsx('C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Finalyield_database.xlsx',sheet = 3))

# combine the three sheets
d1 <- rbind(d1a,d1b,d1c)

# update column names
setnames(d1,
         old = c('clay (%)','Decimal degree lat','Decimal degree lon'),
         new = c('clay','lat_deg','lon_deg'))

# save the unique lon-lat for covariate extraction
# the covariates are collected in script covariate extraction

  # subset only relevant columns, and write csv file, only done once, therefore in FALSE ifelse statement
  if(FALSE){
    d1.lonlat <- d1[,.(lat_deg,lon_deg)]
    d1.lonlat <- unique(d1.lonlat[!is.na(lat_deg)])
    fwrite(d1.lonlat,'data/yield_lonlat_unique.csv')
    
  } else {
    
    # after collection site properties, read these covariates in
    d1.cov <- fread('C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/240523_covariates_yield.csv')
  }
  
  # combine both
  d1 <- merge(d1,d1.cov,by.x = c('lat_deg','lon_deg'),by.y = c('lat','lon'),all.x = TRUE)

# estimate uncertainty when SE and SD are missing
  
  # make kpi_treat_sd from character to numeric
  d1$kpi_treat_sd <- as.numeric(d1$kpi_treat_sd)
  
  # remove SD values exceeding the mean value (checked by MA)
  d1[kpi_treat_sd > kpi_treat, kpi_treat_sd := NA_real_]
  d1[kpi_contr_sd > kpi_control, kpi_contr_sd := NA_real_]
  d1[kpi_treat_se > kpi_treat, kpi_treat_se := NA_real_]
  d1[kpi_contr_se > kpi_control, kpi_contr_se := NA_real_]
  
  # first fill-up dataset for missing input data 
  d1[is.na(kpi_contr_sd) & !is.na(kpi_contr_se), kpi_contr_sd := kpi_contr_se * sqrt(replication)]
  d1[is.na(kpi_treat_sd) & !is.na(kpi_treat_se), kpi_treat_sd := kpi_treat_se * sqrt(replication)]
  d1[!is.na(kpi_contr_sd) & is.na(kpi_contr_se), kpi_contr_se := kpi_contr_sd / sqrt(replication)]
  d1[!is.na(kpi_treat_sd) & is.na(kpi_treat_se), kpi_treat_se := kpi_treat_sd / sqrt(replication)]

  # add coefficient of variation 
  d1[,cv_treat := kpi_treat_sd / kpi_treat]
  d1[,cv_control := kpi_contr_sd / kpi_control]

  # add mean CV for full dataset
  d1[,cv_treat_mean := mean(cv_treat,na.rm=TRUE)]
  d1[,cv_control_mean := mean(cv_control,na.rm=TRUE)]

  # estimate missing SD
  d1[is.na(kpi_treat_sd), kpi_treat_sd := cv_treat_mean * 1.25 * kpi_treat]
  d1[is.na(kpi_contr_sd), kpi_contr_sd := cv_control_mean * 1.25 * kpi_control]

  # remove columns not needed any more
  d1[,c('cv_treat','cv_control','cv_treat_mean','cv_control_mean') := NULL]
  d1[,c('lat_deg','lon_deg','lat','lon','reference','GEnZname','kpi_treat_se','kpi_contr_se') := NULL]
  
  # use only topsoil properties ISRIC and remove subsoil properties
  cols <- colnames(d1)[grepl('_5_15$|15_30$',colnames(d1))]
  d1[,c(cols):=NULL]
  
  # update dataset (MA check that not all data entries are numeric)
  d1[, mat := as.numeric(mat)]
  d1[, map := as.numeric(map)]
  d1[, soc := as.numeric(soc)]
  d1[, duration := as.numeric(gsub('growing seasons|years|year| +','',duration))]
  d1[, year := NULL]
  d1[grepl('<',pH), ph := as.numeric(gsub('>|<','',pH))-0.5]
  d1[grepl('>',pH), ph := as.numeric(gsub('>|<','',pH))+0.5]
  d1[is.na(ph),ph := as.numeric(pH)]  
  d1[,pH := NULL]
  eclay <- function(x){mean(as.numeric(unlist(strsplit(gsub('<|>','',x),'-|~'))))}
  d1[, clay := eclay(clay),by=dataset_ID] # note that ranges are converted into a mean value
  
  # regroup crop type to minimize options
  d1[,crop_type := tolower(crop_type)]
  d1[grepl('barley|wheat|oat|cereal|rye|spelt',crop_type),ctype := 'cereal']
  d1[is.na(ctype) & grepl('vegetable|carrot|onion|melon|tomato|pepper|cabbage',crop_type),ctype := 'vegetable']
  d1[is.na(ctype) & grepl('bean|pea|rape|oil',crop_type),ctype := 'beans_oilcrop']
  d1[is.na(ctype) & grepl('potato|sugar beet',crop_type),ctype := 'arable']
  d1[is.na(ctype) & grepl('maize',crop_type),ctype := 'maize']
  d1[is.na(ctype) & grepl('rice',crop_type),ctype := 'rice']
  d1[is.na(ctype),ctype := 'other']
  
  # regroup fertilizer type
  d1[, fertilizer_type := tolower(fertilizer_type)]
  d1[is.na(fertilizer_type), fertilizer_type := 'unknown']
  d1[grepl('^mineral$|inorganic',fertilizer_type), fertilizer_type := 'inorganic']
  d1[grepl(',|integrat',fertilizer_type), fertilizer_type := 'combined']
  
  # regroup nutricode
  d1[grepl('N, K, P|N/P/K|N, P, K',nutri_code), nutri_code := 'NPK']
  d1[grepl('N, K',nutri_code), nutri_code := 'NK']
  d1[grepl('N, P',nutri_code), nutri_code := 'NP']
  d1[is.na(nutri_code), nutri_code := 'unknown']
  
  # add unknown
  d1[is.na(cover_crop), cover_crop := 'unknown']
  d1[is.na(crop_rotation), crop_rotation := 'unknown']
  d1[is.na(crop_residue), crop_residue := 'unknown']
  d1[is.na(tillage) | tillage == 'yes', tillage := 'CT']
  d1[grepl(',',tillage), tillage := 'CT']
  d1[tillage == 'no', tillage :='NT']
  
  # estimate missing NPKMg doses (MA, please update the input for strange combinations)
  d1[,nutri_dose_N := as.numeric(nutri_dose_N)]
  d1[,nutri_dose_P := as.numeric(nutri_dose_P)]
  d1[,nutri_dose_K := as.numeric(nutri_dose_K)]
  d1[,nutri_dose_Mg := as.numeric(nutri_dose_Mg)]
  d1[is.na(nutri_dose_N),nutri_dose_N := median(d1$nutri_dose_N,na.rm=T)]
  d1[is.na(nutri_dose_P),nutri_dose_P := median(d1$nutri_dose_P,na.rm=T)]
  d1[is.na(nutri_dose_K),nutri_dose_K := median(d1$nutri_dose_K,na.rm=T)]
  d1[is.na(nutri_dose_Mg),nutri_dose_Mg := median(d1$nutri_dose_Mg,na.rm=T)]
  
  # replace missing properties
  d1[is.na(ph), ph := median(d1$ph,na.rm=T)]
  d1[is.na(duration), duration := median(d1$duration,na.rm=T)]
  d1[is.na(mat), mat := median(d1$mat,na.rm=T)]
  d1[is.na(map), map := median(d1$map,na.rm=T)]
  d1[is.na(soc), soc := median(d1$soc,na.rm=T)]
  
  # replace missing replication
  d1[is.na(replication), replication := 2]
  
# make plot sampling location
  
  require(sf)

  # load in the base map
  world <- map_data("world")
  ggplot() +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "#999999", fill = "#CCCCCC", linewidth = 0.1) +
    geom_point(data = d1.cov,aes(lon, lat), alpha = 1, size = 2,color='black') +
    theme_bw()+
    ylim(30,70)+
    xlim(-20,40)+
    coord_sf(crs=4326)
  
# --- Analysis for KPI crop yield -----

  # remove observations without a control or treatment value
  d2 <- d1[!(is.na(kpi_treat) | is.na(kpi_control))]

  # check the data.table on missing values
  summary(d2)

  # add response ratio as effect size
  d2 <- escalc(measure = "ROM", data = d2, 
               m1i = kpi_treat , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control, sd2i= kpi_contr_sd, n2i = replication)

  # convert back to data.table
  d2 <- as.data.table(d2)

  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]

  # make first plot of individual observations
  ggplot(data = d2,aes(x = id,y = yi)) + 
    geom_point() + 
    geom_line() +
    geom_errorbar(aes(x = id,ymin = yi-vi,ymax = yi + vi)) +
    theme_bw() + xlab('study id') + ylab('lnRR, change in crop yield') + 
    ggtitle('Observed changes in crop yield 
            due to crop, soil and 
            fertilisation practices')
  
  d2 <- d2[abs(yi) <= 2]

  # estimate mean response across all studies
m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)

  # see the summary of the model
  summary(m1)
 
# analyse whether the mean response differs per factor
fcols <- c('crop_type','crop_residue','cover_crop','crop_rotation', 'tillage', 
           'fertilizer_type','ph','man_code','fertilizer_type',
           'nutri_dose_N','nutri_dose_P','nutri_dose_K','nutri_dose_Mg',
           'bdod_mean_0_5','cec_mean_0_5','clay_mean_0_5','ntot_mean_0_5',
           'phw_mean_0_5','soc_mean_0_5',
           'tmp_mean','pet_mean','tmp_sd','pet_sd',
           'ctype')

# make an empty list to store output
out.dgr <- list()

# We analyse in a for loop the impact of individual factors explaining variation in response
# this to avoid repetition (otherwise we have to build similar code for each factor seperately)
for(i in fcols){
  
  # make a temporary variable 'evar'
  d2[,evar := get(i)]

  # is the variable categorial or not
  vartype = is.character(d2[,get(i)])

  # do the meta-analysis for categorial variable
  if(vartype==TRUE){
    
    # do the meta-analysis for categorial variable
    m2 <- metafor::rma.mv(yi,vi, mods = ~factor(evar),data=d2,
                          random= list(~ 1|study_ID), method="REML",sparse = TRUE)
    
    # collect model coefficients
    m2.sum <- summary(m2)
    m2.out <- as.data.table(coefficients(m2.sum))
    m2.out[,factor := i]
    m2.out[,var := gsub('factor\\(evar\\)','',rownames(m2.sum$b))]
    
  } else{
    
  # do the meta-analysis for a numerical variable
  m2 <- metafor::rma.mv(yi,vi, mods = ~evar,data=d2,
                        random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # collect model coefficients
  m2.sum <- summary(m2)
  m2.out <- as.data.table(coefficients(m2.sum))
  m2.out[,factor := i]
  m2.out[,var := c('intercept',i)]
  
  }
  
  # add to list
  out.dgr[[i]] <- copy(m2.out)
  
  # print message
  print(paste0('model build for parameter: ',i))
}

# convert list to one data.table
out.dgr <- rbindlist(out.dgr)

# remove cropname (too mucht options, use only grouped one)
out.dgr <- out.dgr[factor != 'crop_type']
# make barplot (figure formatting: need to be done)
#out.dgr[, pfactor := as.factor(factor)]
#ggplot(data=out.dgr[var!='intercept']) +
  #geom_bar(aes(x=estimate,y= var),stat="identity") + 
  #geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
  #ggtitle('Effect of main factors on Crop yield') +
  #theme(legend.position = 'bottom')

#Create a barplot for the management practices and their estimates
library(dplyr)

# Filter the data

filtered_data <- out.dgr %>%
  filter(factor == 'man_code' & !grepl('intercept|intrcpt', var))

# Plot the filtered data of management practices
ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = var), stat = "identity") + 
  geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  theme_bw() +
  ggtitle('Effect of management practices 
          on Crop yield') +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 10, hjust = 0),  # Adjust title size and align to left
    axis.title.x = element_text(size = 12),         # Adjust x-axis title size
    axis.title.y = element_text(size = 12),         # Adjust y-axis title size
    axis.text.x = element_text(size = 8),           # Adjust x-axis text size
    axis.text.y = element_text(size = 6, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )
#Create a barplot for the fertilisation doses and their estimates
# Filter the data
filtered_data <- out.dgr %>%
  filter(factor %in% c('nutri_dose_N', 'nutri_dose_P', 'nutri_dose_Mg','nutri_dose_K') & !grepl('intercept|intrcpt', var))


# Add asterisks for significance
filtered_data <- filtered_data %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Plot the filtered data of fertilization dose type
ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = var), stat = "identity") + 
  geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  geom_text(aes(x = estimate, y = var, label = significance), 
            vjust = -0.5, hjust = -0.2, size = 5, color = "red") +
  theme_bw() +
  ggtitle('Effect of fertilization dose on Crop yield') +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 10, hjust = 0),  # Adjust title size and align to left
    axis.title.x = element_text(size = 12),           # Adjust x-axis title size
    axis.title.y = element_text(size = 12),           # Adjust y-axis title size
    axis.text.x = element_text(size = 8),             # Adjust x-axis text size
    axis.text.y = element_text(size = 6, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )

#Create a barplot for the site properties and their estimates
# Filter the data
  filtered_data <- out.dgr %>%
  filter(factor %in% c('ph','bdod_mean_0_5', 'cec_mean_0_5','clay_mean_0_5',
                       'ntot_mean_0_5','soc_mean_0_5', 'phw_mean_0_5', 'tmp_mean','pet_mean') & !grepl('intercept|intrcpt', var))


# Add asterisks for significance
filtered_data <- filtered_data %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Plot the filtered data of site properties type
ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = var), stat = "identity") + 
  geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  geom_text(aes(x = estimate, y = var, label = significance), 
            vjust = -0.5, hjust = -0.2, size = 5, color = "red") +
  theme_bw() +
  ggtitle('Effect of site properties on Crop yield') +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 10, hjust = 0),  # Adjust title size and align to left
    axis.title.x = element_text(size = 12),           # Adjust x-axis title size
    axis.title.y = element_text(size = 12),           # Adjust y-axis title size
    axis.text.x = element_text(size = 8),             # Adjust x-axis text size
    axis.text.y = element_text(size = 6, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )

##Barplot for crop type
# Filter the data
filtered_data <- out.dgr %>%
  filter(factor == 'ctype' & !grepl('intercept|intrcpt', var))

# Add asterisks for significance
filtered_data <- filtered_data %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Plot the filtered data of crop type
ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = var), stat = "identity") + 
  geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  geom_text(aes(x = estimate, y = var, label = significance), 
            vjust = -0.5, hjust = -0.2, size = 5, color = "red") +
  theme_bw() +
  ggtitle('Effect of crop type on Crop yield') +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 10, hjust = 0),  # Adjust title size and align to left
    axis.title.x = element_text(size = 12),           # Adjust x-axis title size
    axis.title.y = element_text(size = 12),           # Adjust y-axis title size
    axis.text.x = element_text(size = 8),             # Adjust x-axis text size
    axis.text.y = element_text(size = 6, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )

# do one meta-regression model with multiple factors

# make a function to extract relevant model statistics
estats <- function(model_new,model_base){
  out <- data.table(AIC = model_new$fit.stats[4,2],
                    ll = model_new$fit.stats[1,2],
                    ll_impr = round(100 * (1-model_new$fit.stats[1,2]/model_base$fit.stats[1,2]),2),
                    r2_impr = round(100*max(0,(sum(model_base$sigma2)-sum(model_new$sigma2))/sum(model_base$sigma2)),2),
                    pval = round(anova(model_new,model_base)$pval,3))
  return(out)
}

# make first an empty model
m3.empty <- metafor::rma.mv(yi,vi, data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
#make the different practice groups
combination_practices <- c("CC + M","CC + MCR + RR", "CC + MCR + RR + CT", "CC + MCR + RR + RT", "CT + CC + M", "RT + NF", "FT + NF","IR + FT","RR + RT","IMP + NF","MCR + NF", "DD + ML", data = d2)
tillage_practices <- c("NT", "ST", "RT", data = d2)
fertilisation_practices <- c("DF", "NF", "MF", data = d2)
irrigation_practices <- c("IR", data = d2)
mulching_practices <- "ML"
residue_practices <- c("RR", data = d2)
biochar_practices <- c("IMP", data = d2)
crop_practices <- c("CC", "INC", "MCR", data = d2)
drilling_practices <- c("DD", data = d2)

#create new variable indicating the groups
d2$man_code_new <- ifelse(d2$man_code %in% combination_practices, "Combination practices", 
                          ifelse(d2$man_code %in% tillage_practices, "Tillage Practices", 
                                 ifelse(d2$man_code %in% fertilisation_practices, "Fertilisation practices",
                                        ifelse(d2$man_code %in% irrigation_practices, "Irrigation practices",
                                               
                                                      ifelse(d2$man_code %in% residue_practices, "Residue practices",
                                                             ifelse(d2$man_code %in% biochar_practices, "Biochar practices",
                                                                    ifelse(d2$man_code %in% crop_practices, "Crop practices",
                                                                           ifelse(d2$man_code %in% drilling_practices, "Drilling practices", "Other"))))))))

# Print the updated data to verify the new column
print(d2)

# Function to standardize numeric columns
library(dplyr)
standardize <- function(x) {
  return ((x - mean(x)) / sd(x))
}
# Apply the standardization function to numeric columns
d2_standardized <- d2 %>%
  mutate(across(where(is.numeric), standardize))

# Display the standardized data frame
print(d2_standardized)

# Fit the model using the new grouping variable
m3.full <- rma.mv(yi, vi, 
                  mods = ~ man_code_new + ctype + crop_residue + cover_crop + crop_rotation + 
                    tillage + fertilizer_type + nutri_dose_N + nutri_dose_P + nutri_dose_K +
                    nutri_dose_Mg + ntot_mean_0_5 +
                    ph + tmp_mean + pet_mean,
                  data = d2, 
                  random = list(~ 1 | study_ID), 
                  method = "REML", 
                  sparse = TRUE)
# analyse summary stats
summary(m3.full)


#different grouping and final refined model
#make the different practice groups
combination_practices <- c("CC + MCR + RR", "CC + M", "CC + MCR + RR + CT", "CC + MCR + RR + RT", "CT + CC + M", "RT + NF", "FT + NF","IR + FT","RR + RT","IMP + NF","MCR + NF", "DD + ML", data = d2)
soil_practices <- c("NT", "ST", "RT","ML", "IMP","RR","DD", data = d2)
fertilisation_practices <- c("DF", "NF", "MF", data = d2)
irrigation_practices <- c("IR", data = d2)
crop_practices <- c("CC", "INC", "MCR", data = d2)


#create new variable indicating the groups
d2$man_code_new2 <- ifelse(d2$man_code %in% combination_practices, "Combination practices", 
                          ifelse(d2$man_code %in% soil_practices, "Soil Practices", 
                                 ifelse(d2$man_code %in% fertilisation_practices, "Fertilisation practices",
                                        ifelse(d2$man_code %in% irrigation_practices, "Irrigation practices", "Other"))))

# Print the updated data to verify the new column
print(d2)

# Function to standardize numeric columns
library(dplyr)
standardize <- function(x) {
  return ((x - mean(x)) / sd(x))
}
# Apply the standardization function to numeric columns
d2_standardized <- d2 %>%
  mutate(across(where(is.numeric), standardize))

# Display the standardized data frame
print(d2_standardized)

m3.full1 <- rma.mv(yi, vi, 
                  mods = ~ man_code_new2 + ctype + cover_crop + crop_rotation + 
                   nutri_dose_N + nutri_dose_K +
                     +ntot_mean_0_5 +
                    ph + tmp_mean + pet_mean,
                  data = d2, 
                  random = list(~ 1 | study_ID), 
                  method = "REML", 
                  sparse = TRUE)

# analyse summary stats
summary(m3.full1)


# I DO NOT UNDERSTAND THIS PART refine the model (if desired)
# steps to do: add variables plus interactions (* for all interactions plus main factors, : for interactions only), check pvalue
# if within a variable one subgroup is significant and others not, then adjust the groupings
# so, you have the options:
# only a model with additive factors => use the plus sign in the mods argument
# if you have a model with main factors AND interactions => use the "*" sign in the mods argument
# if you have a model with ONLY interactions => use then the ":" sign
#d2[,sup_cat := fifelse(grepl('micro',crop_type),'x','other')]
#m3.full <- metafor::rma.mv(yi,vi, 
 #                          mods = ~ supplemental_rate + stage : sup_cat,
 #                          data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)

# It doesn't work probably because of above collect stats of the model
estats(m3.full,m3.full1)

# It doesn't work probabl ybecause of above show anova whether full model is signifiantly different from empty model
anova(m3.full,m3.full1,refit = T)