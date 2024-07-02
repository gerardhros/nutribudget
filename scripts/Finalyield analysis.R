# data analysis RQ1
# Maria Anna Antonovardaki & Gerard Ros
# 22 May 2024

#required packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)
require(broom)
require(dplyr)
require(sf)

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
# load in the base map
  world <- map_data("world")
  ggplot() +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "#999999", fill = "#CCCCCC", linewidth = 0.1) +
    geom_point(data = d1.cov,aes(lon, lat), alpha = 1, size = 2,shape =21, fill = 'green', color='black') +
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
    #scale_color_gradient(low="red", high="green") +   # Sequential color scheme
    geom_point() + 
    geom_line() +
    geom_errorbar(aes(x = id,ymin = yi-vi,ymax = yi + vi)) +
    theme_bw() + xlab('study id') + ylab('lnRR, change in crop yield') + 
    ggtitle('Observed changes in crop yield due to crop, soil and fertilisation practices') +
    theme(plot.title = element_text(size = rel(0.8)),
          axis.title.x = element_text(size = rel(0.8)),  # Adjust x-axis label size
          axis.title.y = element_text(size = rel(0.8)))
  
  d2 <- d2[abs(yi) <= 2]

  # random effect model estimate mean response across all studies
m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)

  # see the summary of the model
  summary(m1)
 
# mian factor analysis analyse whether the mean response differs per factor
fcols <- c('crop_type','crop_residue','cover_crop','crop_rotation', 'tillage', 
           'fertilizer_type','man_code','fertilizer_type',
           'nutri_dose_N','nutri_dose_P','nutri_dose_K','nutri_dose_Mg',
           'bdod_mean_0_5','cec_mean_0_5','clay_mean_0_5','ntot_mean_0_5',
           'phw_mean_0_5','soc_mean_0_5',
           'tmp_mean','pet_mean','tmp_sd','pet_sd','pre_mean','pre_sd',
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
    m2 <- metafor::rma.mv(yi,vi, mods = ~factor(evar)-1,data=d2,
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

##Create a barplot for the management practices and their estimates
# Filter the data

filtered_data <- out.dgr %>%
  filter(factor == 'man_code' & !grepl('intercept|intrcpt', var))

# Add asterisks for significance
filtered_data <- filtered_data %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ))
practice_order <- filtered_data$var

# Plot the filtered data of management practices
p1 <- ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = var), stat = "identity", color='lightblue', fill = 'lightblue') + 
  geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  theme_bw() +
  geom_text(aes(x = estimate, y = var, label = significance), vjust = +0.4, hjust = -0.2, size = 5, color = "red") +  # Add the asterisks for significance
  ggtitle(      'Effect of management practices on crop yield') + 
  labs(x = 'lnRR response of crop yield', y = '') +  # Update x and y axis titles
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),  # Adjust title size and align to left
    axis.title.x = element_text(size = 10),           # Adjust x-axis title size
    axis.title.y = element_text(size = 12),           # Adjust y-axis title size
    axis.text.x = element_text(size = 8),             # Adjust x-axis text size
    axis.text.y = element_text(size = 8, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )

#Export and save the plot
# Save the plot to the specified directory
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p1, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

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
p2 <- ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = var), stat = "identity", color = "lightblue", fill = "lightblue") + 
  geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  geom_text(aes(x = estimate, y = var, label = significance), 
            vjust = -0.5, hjust = -0.2, size = 6, color = "red") +
  theme_bw() +
  ggtitle('Effect of fertilization doses on management induced changes in crop yield') +
  labs(x = 'lnRR response of crop yield', y = '') +  # Update x and y axis titles
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 12, hjust = 0.5),  # Adjust title size and align to left
    axis.title.x = element_text(size = 10),           # Adjust x-axis title size
    axis.title.y = element_text(size = 12),           # Adjust y-axis title size
    axis.text.x = element_text(size = 8),             # Adjust x-axis text size
    axis.text.y = element_text(size = 9, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )
#Export and save the plot
# Save the plot to the specified directory
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p2, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

#Create a barplot for the site properties and their estimates
# Filter the data
filtered_data <- out.dgr %>%
  filter(factor %in% c('crop_residue', 'cover_crop', 'crop_rotation', 'tillage', 'fertilizer_type', 
                        'bdod_mean_0_5', 'cec_mean_0_5', 'clay_mean_0_5', 'ntot_mean_0_5',
                       'soc_mean_0_5', 'phw_mean_0_5', 'tmp_mean', 'pre_mean', 'pet_mean', 'ctype') & 
           !grepl('intercept|intrcpt', var))

# Add descriptive labels
filtered_data <- filtered_data %>%
  mutate(label = case_when(
    factor == 'crop_residue' & var == 'no' ~ 'crop residue (no)',
    factor == 'crop_residue' & var == 'yes' ~ 'crop residue (yes)',
    factor == 'crop_residue' & var == 'unknown' ~ 'crop residue (unknown)',
    factor == 'cover_crop' & var == 'no' ~ 'cover crop (no)',
    factor == 'cover_crop' & var == 'yes' ~ 'cover crop (yes)',
    factor == 'cover_crop' & var == 'unknown' ~ 'cover crop (unknown)',
    factor == 'crop_rotation' & var == 'no' ~ 'crop rotation (no)',
    factor == 'crop_rotation' & var == 'yes' ~ 'crop rotation (yes)',
    factor == 'crop_rotation' & var == 'unknown' ~ 'crop rotation (unknown)',
    factor == 'tillage' & var == 'CT' ~ 'tillage (CT)',
    factor == 'tillage' & var == 'NT' ~ 'tillage (NT)',
    factor == 'tillage' & var == 'RT' ~ 'tillage (RT)',
    factor == 'fertilizer_type' & var == 'combined' ~ 'fertilizer type (combined)',
    factor == 'fertilizer_type' & var == 'inorganic' ~ 'fertilizer type (inorganic)',
    factor == 'fertilizer_type' & var == 'organic' ~ 'fertilizer type (organic)',
    factor == 'fertilizer_type' & var == 'unknown' ~ 'fertilizer type (unknown)',
    factor == 'bdod_mean_0_5' ~ 'bulk density',
    factor == 'cec_mean_0_5' ~ 'cation exchange capacity',
    factor == 'clay_mean_0_5' ~ 'clay content',
    factor == 'ntot_mean_0_5' ~ 'total nitrogen',
    factor == 'soc_mean_0_5' ~ 'soil organic carbon',
    factor == 'phw_mean_0_5' ~ 'water pH',
    factor == 'tmp_mean' ~ 'temperature mean',
    factor == 'pre_mean' ~ 'precipitation mean',
    factor == 'pet_mean' ~ 'evapotranspiration',
    factor == 'ctype' & var == 'arable' ~ 'arable',
    factor == 'ctype' & var == 'beans_oilcrop' ~ 'beans & oilcrop',
    factor == 'ctype' & var == 'cereal' ~ 'cereals',
    factor == 'ctype' & var == 'maize' ~ 'maize',
    factor == 'ctype' & var == 'other' ~ 'other',
    factor == 'ctype' & var == 'rice' ~ 'rice',
    factor == 'ctype' & var == 'vegetable' ~ 'vegetables'
  ))

# Add asterisks for significance
filtered_data <- filtered_data %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Define the order for the labels
ordered_labels <- c(
  'crop residue (no)', 'crop residue (yes)', 'crop residue (unknown)',
  'cover crop (no)', 'cover crop (yes)', 'cover crop (unknown)',
  'crop rotation (no)', 'crop rotation (yes)', 'crop rotation (unknown)',
  'tillage (CT)', 'tillage (NT)', 'tillage (RT)',
  'fertilizer type (combined)', 'fertilizer type (inorganic)', 'fertilizer type (organic)', 'fertilizer type (unknown)',
 'bulk density', 'cation exchange capacity', 'clay content', 'total nitrogen',
  'soil organic carbon', 'water pH', 'temperature mean', 'precipitation mean', 'evapotranspiration',
  'arable', 'beans & oilcrop', 'cereals', 'maize', 'other', 'rice', 'vegetables'
)

# Convert the label column to a factor with the specified levels
filtered_data <- filtered_data %>%
  mutate(label = factor(label, levels = ordered_labels))

# Plot the filtered data with descriptive labels
p3 <- ggplot(data = filtered_data) +
  geom_bar(aes(x = estimate, y = label), stat = "identity", color = 'lightblue', fill = 'lightblue') + 
  geom_errorbar(aes(y = label, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
  theme_bw() +
  geom_text(aes(x = estimate, y = label, label = significance), vjust = +0.4, hjust = -0.2, size = 4.5, color = "red") + 
  ggtitle('Effect of site properties on crop yield') + 
  labs(x = 'lnRR response of crop yield', y = '') +  # Update x and y axis titles
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),  # Adjust title size and align to left
    axis.title.x = element_text(size = 10),           # Adjust x-axis title size
    axis.title.y = element_text(size = 12),           # Adjust y-axis title size
    axis.text.x = element_text(size = 8),             # Adjust x-axis text size
    axis.text.y = element_text(size = 9, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
  )
#Export and save the plot
# Save the plot to the specified directory
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p3, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)


## do one meta-regression model with multiple factors

# mixed effect model make a function to extract relevant model statistics
estats <- function(model_new,model_base){
out <- data.table(AIC = model_new$fit.stats[4,2],
ll = model_new$fit.stats[1,2],
ll_impr = round(100 * (1-model_new$fit.stats[1,2]/model_base$fit.stats[1,2]),2),
r2_impr = round(100*max(0,(sum(model_base$sigma2)-sum(model_new$sigma2))/sum(model_base$sigma2)),2),
pval = round(anova(model_new,model_base)$pval,3))
return(out)
}

## make first an empty model
m3.empty <- metafor::rma.mv(yi,vi, data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
##make the different practice groups
combination_practices <- c("CC + M","CC + MCR + RR", "CC + MCR + RR + CT", "CC + MCR + RR + RT", "CT + CC + M", "RT + NF", "FT + NF","IR + FT","RR + RT","IMP + NF","MCR + NF", "DD + ML", "CC + RR + CT", "DD + MCR")
tillage_practices <- c("NT", "ST", "RT")
fertilisation_practices <- c("DF", "NF", "MF", "OF")
irrigation_practices <- c("IR")
mulching_practices <- c("ML")
residue_practices <- c("RR")
biochar_practices <- c("IMP")
crop_practices <- c("CC", "INC", "MCR", "L")
drilling_practices <- c("DD")

#create new variable indicating the groups
d2$man_code_new <- ifelse(d2$man_code %in% combination_practices, "Combination practices", ifelse(d2$man_code %in% tillage_practices, "Tillage Practices", 
ifelse(d2$man_code %in% irrigation_practices, "Irrigation practices",
ifelse(d2$man_code %in% fertilisation_practices, "Fertilisation practices",
ifelse(d2$man_code %in% residue_practices, "Residue practices",
ifelse(d2$man_code %in% biochar_practices, "Biochar practices",
ifelse(d2$man_code %in% crop_practices, "Crop practices",
ifelse(d2$man_code %in% mulching_practices, "Mulching practices", 
       "Drilling practices"))))))))

# Print the updated data to verify the new column
print(d2)

# a non linear response of crop yield due to n dose
d2$n_dose <- (d2$nutri_dose_N)^2


# Function to standardize numeric columns
standardize <- function(x) {
  return ((x - mean(x)) / sd(x))
}

cols <- c('n_dose', 'nutri_dose_P', 'nutri_dose_K',
          'nutri_dose_Mg', 'ntot_mean_0_5', 'tmp_mean', 'pre_mean', 'pet_mean', 'soc_mean_0_5', 'clay_mean_0_5', 'bdod_mean_0_5', 'cec_mean_0_5', 'phw_mean_0_5' )
d2[,c(cols) := lapply(.SD,standardize),.SDcols=cols]


# Fit the model using the new grouping variable
m3.full <- rma.mv(yi, vi, 
                  mods = ~ man_code_new + ctype + crop_residue + cover_crop + crop_rotation + tillage + fertilizer_type + n_dose + nutri_dose_P + nutri_dose_K +
nutri_dose_Mg + ntot_mean_0_5 +
 tmp_mean + pre_mean + pet_mean +  soc_mean_0_5 + clay_mean_0_5 + bdod_mean_0_5 + cec_mean_0_5 + phw_mean_0_5 -1,
                  data = d2, 
                  random = list(~ 1 | study_ID), 
                  method = "REML", 
                  sparse = TRUE)
# analyse summary stats
summary(m3.full)



##Barplot of practices & site properties for the full model
# Extract and tidy the results
m3_tidy <- broom::tidy(m3.full)

# Define new column values
new_column_values <- c("management_practice","management_practice", "management_practice", "management_practice","management_practice","management_practice",
                       "management_practice","management_practice","management_practice",
                       "crop_type","crop_type","crop_type","crop_type",
                       "crop_type","crop_type","crop_residue", "crop_residue", 
                       "cover_crop", "cover_crop", "crop_rotation", 
                       "crop_rotation", "tillage", "tillage", "fertilizer_type ",
                       "fertilizer_type", "fertilizer_type","nutri_does", 
                       "nutri_dose", "nutri_dose", "nutri_dose", "soil_property", 
                       "soil_property", "weather_condition", "weather_condition", "weather_condition", "soil_property", "soil_property","soil_property", "soil_property")

m3_tidy$colorfactor <- new_column_values

# Rename factor levels using recode with named arguments
m3_tidy <- m3_tidy %>%
  mutate(factor = recode(term,
                         `man_code_newCombination practices` = "Combination",
                         `man_code_newBiochar practices` = "Biochar",  
                         `man_code_newCrop practices` = "Crop",
                         `man_code_newDrilling practices` = "Drilling",
                         `man_code_newFertilisation practices` = "Fertilisation",
                         `man_code_newIrrigation practices` = "Irrigation",
                         `man_code_newMulching practices` = "Mulching",
                         `man_code_newResidue practices` = "Residue",
                         `man_code_newTillage Practices` = "Tillage",
                         `man_code_newOther` = "Other practices",
                         `ctypebeans_oilcrop` = "Beans & oilcrop",
                         `ctypecereal` = "Cereals",
                         `ctypemaize` = "Maize",
                         `ctypeother` = "Other",
                         `ctyperice` = "Rice",
                         `ctypevegetable` = "Vegetables",
                         `crop_residueunknown` = "Crop residue (U)", 
                         `crop_residueyes` = "Crop residue",
                         `cover_cropunknown` = "Cover crop (U)",
                         `cover_cropyes` = "Cover crop", 
                         `crop_rotationunknown` = "Crop rotation (U)",
                         `crop_rotationyes` = "Crop rotation", 
                         `tillageNT` = "No tillage", 
                         `tillageRT` = "Reduced tillage", 
                         `fertilizer_typeinorganic ` = "Inorganic fertiliser",
                         `fertilizer_typeorganic` = "Organic fertiliser",
                         `fertilizer_typeunknown` = "Fertiliser (U)",
                         `n_dose` = "N dose",
                         `nutri_dose_P` = "P dose", 
                         `nutri_dose_K` = "K dose", 
                         `nutri_dose_Mg` = "Mg dose", 
                         `ntot_mean_0_5` = "Total nitrogen", 
                         `tmp_mean` = "Temp", 
                         `pre_mean` = "Precipitation",
                         `pet_mean` = "Evapotranspiration",
                         `soc_mean_0_5`  = "Soil organic C", 
                         `clay_mean_0_5` = "Clay",
                         `bdod_mean_0_5` = "Bulk density", 
                         `cec_mean_0_5` = "Cation exchange capacity", 
                         `phw_mean_0_5` = "pH Water"))

# Add a column for significance stars
m3_tidy$significance <- ifelse(m3_tidy$p.value < 0.001, "***", 
ifelse(m3_tidy$p.value < 0.01, "**", 
ifelse(m3_tidy$p.value < 0.05, "*", 
ifelse(m3_tidy$p.value < 0.1, ".", ""))))

# Define the order of the terms
m3_tidy$factor <- factor(m3_tidy$factor, levels = unique(m3_tidy$factor))

# Split data into management practices and site properties
management_practices <- m3_tidy[m3_tidy$factor %in% c("Combination", "Crop", "Drilling", "Fertilisation","Irrigation", "Mulching", "Residue", "Tillage", "Biochar"), ]

site_properties <- m3_tidy[!m3_tidy$factor %in% c("Intrecept", "Combination", "Crop", "Drilling", "Fertilisation", "Irrigation","Mulching", "Residue", "Tillage", "Other practices", "Biochar"), ]

# Plot for management practices with Standard error (SE)
p4 <- ggplot(management_practices, aes(x = factor, y = estimate)) + 
  geom_bar(stat = "identity", fill = "lightblue", color = NA, alpha = 0.7) +  # Set fill to lightblue and remove black border
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
            vjust = ifelse(management_practices$estimate > 0, -0.5, 1.5), size = 4) + 
  theme_minimal() + 
  labs(x = 'Crop, soil and fertilisation practices', y = 'Parameter estimate', 
       title = 'Meta-regression model estimation on crop yield') + 
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1), 
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.ticks.length = unit(-0.5, "cm"), 
    panel.grid = element_blank(), 
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),  # Set plot background to white
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) + 
  ylim(-1, 1)

# Save the plot
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p4, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

# Plot for site properties with SE
p5 <- ggplot(site_properties, aes(x = factor, y = estimate)) + 
  geom_bar(stat = "identity", aes(fill = colorfactor), alpha = 0.7) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
            vjust = ifelse(site_properties$estimate > 0, -0.5, 1.5), size = 4) + 
  theme_minimal() + 
  labs(x = 'Site properties', y = 'Parameter estimate', 
       title = 'Meta-regression model estimation on crop yield') + 
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7.6), 
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.ticks.length = unit(-0.5, "cm"), 
    panel.grid = element_blank(), 
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),  # Set plot background to white
    plot.title = element_text(size = 14, hjust = 0.5)  # Center the title and adjust the size
  ) + 
  ylim(-0.7, 1.1)

# Save the plot
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p5, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

##different grouping and final refined model
#make the different practice groups
combination_practices1 <- c("CC + MCR + RR", "CC + M", "CC + MCR + RR + CT", "CC + MCR + RR + RT", "CT + CC + M", "RT + NF", "FT + NF","IR + FT","RR + RT","IMP + NF","MCR + NF", "DD + ML")
soil_practices1 <- c("NT", "ST", "RT","ML", "IMP","RR","DD")
fertilisation_practices1 <- c("DF", "NF", "MF", "OF")
irrigation_practices1 <- c("IR")
crop_practices1 <- c("CC", "INC", "MCR")


#create new variable indicating the groups, here the other practices are the crop practices
d2$man_code_new2 <- ifelse(d2$man_code %in% combination_practices1, "Combination practices",ifelse(d2$man_code %in% soil_practices1, "Soil Practices",
ifelse(d2$man_code %in% fertilisation_practices1, "Fertilisation practices",
ifelse(d2$man_code %in% crop_practices1, "Crop practices",
ifelse(d2$man_code %in% irrigation_practices1, "Irrigation practices", "Other")))))

# Print the updated data to verify the new column
print(d2)

# Function to standardize numeric columns
standardize <- function(x) {
  return ((x - mean(x)) / sd(x))
}

cols <- c('n_dose', 'nutri_dose_P', 'nutri_dose_K',
          'nutri_dose_Mg', 'ntot_mean_0_5', 'tmp_mean', 'pre_mean', 'pet_mean', 'soc_mean_0_5', 'clay_mean_0_5', 'bdod_mean_0_5', 'cec_mean_0_5', 'phw_mean_0_5' )
d2[,c(cols) := lapply(.SD,standardize),.SDcols=cols]

# Include additive values  between management practices and site properties
m3.full1 <- rma.mv(yi, vi, 
                  mods = ~ man_code_new2 + ctype + cover_crop + crop_rotation + 
                    + clay_mean_0_5 + bdod_mean_0_5 + cec_mean_0_5 + phw_mean_0_5 
                     + ntot_mean_0_5 + nutri_dose_P + pre_mean -1,
                  data = d2, 
                  random = list(~ 1 | study_ID), 
                  method = "REML", 
                  sparse = TRUE)

# analyse summary stats
summary(m3.full1)

##Barplot of practices & site properties for the refined model
# Extract and tidy the results
m3_tidy <- broom::tidy(m3.full1)

# Define new column values
new_column_values1 <- c("management_practice","management_practice", "management_practice", "management_practice","management_practice",
                       "crop_type","crop_type","crop_type","crop_type",
                       "crop_type","crop_type", 
                       "cover_crop", "cover_crop", "crop_rotation", 
                       "crop_rotation", "soil_property", 
                       "soil_property", "soil_property", "soil_property","soil_property", "nutrient_dose", "weather_condition")

m3_tidy$colorfactor1 <- new_column_values1


# Rename terms to more readable names
m3_tidy <- m3_tidy %>%
  mutate(factor = recode(term,
                         `man_code_new2Combination practices` = "Combination",
                         `man_code_new2Crop practices` = "Crop",
                         `man_code_new2Fertilisation practices` = "Fertilisation",
                         `man_code_new2Irrigation practices` = "Irrigation",
                         `man_code_new2Soil Practices` = "Soil",
                         `ctypebeans_oilcrop` = "Beans & oilcrop",
                         `ctypecereal` = "Cereals",
                         `ctypemaize` = "Maize",
                         `ctypeother` = "Other",
                         `ctyperice` = "Rice",
                         `ctypevegetable` = "Vegetables",
                         `cover_cropunknown` = "Cover crop (U)",
                         `cover_cropyes` = "Cover crop",
                         `crop_rotationunknown` = "Crop rotation (U)",
                         `crop_rotationyes` = "Crop rotation",
                         `ntot_mean_0_5` = "Total nitrogen",
                         `pre_mean` = "Precipitation",
                         `clay_mean_0_5` = "Clay",
                         `bdod_mean_0_5` = "Bulk density",
                         `cec_mean_0_5` = "Cation exchange capacity",
                         `phw_mean_0_5` = "pH Water"))

# Add a column for significance stars
m3_tidy$significance <- ifelse(m3_tidy$p.value < 0.001, "***", 
                                ifelse(m3_tidy$p.value < 0.01, "**", 
                                       ifelse(m3_tidy$p.value < 0.05, "*", 
                                              ifelse(m3_tidy$p.value < 0.1, ".", ""))))

# Define the order of the terms
m3_tidy$factor <- factor(m3_tidy$factor, levels = unique(m3_tidy$factor))

# Split data into management practices and site properties
management_practices <- m3_tidy %>%
filter(factor %in% c("Combination", "Crop", "Fertilisation", "Irrigation", "Soil"))

site_properties <- m3_tidy %>%
filter(!factor %in% c("Combination", "Crop", "Fertilisation", "Irrigation", "Soil"))

# Plot for management practices with SE
p6 <- ggplot(management_practices, aes(x = factor, y = estimate)) + 
  geom_bar(stat = "identity", fill = "lightblue", color = NA, alpha = 0.7) +  # Set fill to lightblue and remove black border
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
            vjust = ifelse(management_practices$estimate > 0, -0.5, 1.5), size = 4) + 
  theme_minimal() + 
  labs(x = 'Crop, soil and fertilisation practices', y = 'Parameter estimate', 
       title = 'Meta-regression model estimation on crop yield') + 
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1), 
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.ticks.length = unit(-0.5, "cm"), 
    panel.grid = element_blank(), 
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),  # Set plot background to white
    plot.title = element_text(size = 14, hjust = 0.5)  # Center the title and adjust the size
  ) + 
  ylim(-0.8, 0.8)

# Save the plot
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p6, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

# Plot for site properties with SE
p7 <- ggplot(site_properties, aes(x = factor, y = estimate)) + 
  geom_bar(stat = "identity", aes(fill = colorfactor1), alpha = 0.7) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
            vjust = ifelse(site_properties$estimate > 0, -0.5, 1.5), size = 4) + 
  theme_minimal() + 
  labs(x = 'Site properties', y = 'Parameter estimate', 
       title = 'Meta-regression model estimation on crop yield') + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.ticks.length = unit(-0.5, "cm"), 
        panel.grid = element_blank(), 
        legend.position = "none") + 
  ylim(-0.4, 0.45)
# Save the plot
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p7, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

p7 <- ggplot(site_properties, aes(x = factor, y = estimate)) + 
  geom_bar(stat = "identity", aes(fill = colorfactor1), alpha = 0.7) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
            vjust = ifelse(site_properties$estimate > 0, -0.5, 1.5), size = 4) + 
  theme_minimal() + 
  labs(x = 'Site properties', y = 'Parameter estimate', 
       title = 'Meta-regression model estimation on crop yield') + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.ticks.length = unit(-0.5, "cm"), 
        panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),  # Adjust title size, center, and color
        plot.background = element_rect(fill = "white", color = NA)) +  # Set plot background to white
  ylim(-0.4, 0.45)

# Save the plot
ggsave("C:/Users/anton040/OneDrive - Wageningen University & Research/PHD/Nutribudget/Task 1.2 Template/Final Excel for Crop yields/Final/plot.png", 
       plot = p7, 
       device = "png", 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)
# Estats of all models 
estats(m3.empty,m3.full1)
estats(m3.full,m3.full1)


# Empty model compared to the different full models
anova(m3.empty,m3.full,refit = T)
anova(m3.empty,m3.full1,refit = T)


# Comparison of first full model with the refined model
anova(m3.full,m3.full1,refit = T)