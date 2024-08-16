# data analysis RQ1
# Maria Anna Antonovardaki & Gerard Ros
# 22 May 2024

# clean environment
rm(list=ls())

#required packages
require(readxl);require(data.table);require(metafor)
require(ggplot2);require(sf);require(broom)

# --- load in databases and clean up data ----

  # read in the database (NOTE that I updated some SE and SD values)
  d1a <- as.data.table(read_xlsx('data/Finalyield_database.xlsx',sheet = 1))
  d1b <- as.data.table(read_xlsx('data/Finalyield_database.xlsx',sheet = 2))
  d1c <- as.data.table(read_xlsx('data/Finalyield_database.xlsx',sheet = 3))
  
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
    d1.cov <- fread('products/240523_covariates_yield.csv')
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
  if(FALSE){
    
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
  }
  
   
  
# --- Analysis for KPI crop yield -----

  # remove observations without a control or treatment value
  d2 <- d1[!(is.na(kpi_treat) | is.na(kpi_control))]

  # check the data.table on missing values
  # summary(d2)

  # add response ratio as effect size
  d2 <- escalc(measure = "ROM", data = d2, 
               m1i = kpi_treat , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control, sd2i= kpi_contr_sd, n2i = replication)

  # convert back to data.table
  d2 <- as.data.table(d2)
  
  if(FALSE){
    
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
  }
  
  
  d2 <- d2[abs(yi) <= 2]

  # random effect model estimate mean response across all studies
  m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)

  # see the summary of the model
  summary(m1)
 
# --- main factor analysis -----
  
  # main factor analysis analyse whether the mean response differs per factor
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

# --- plots main factor analysis ----
  
  if(FALSE){
    
  ## Create a barplot for the management practices and their estimates

  # Filter the data
  filtered_data <- out.dgr[factor == 'man_code' & !grepl('intercept|intrcpt', var)]

  # Add asterisks for significance
  filtered_data[pval > 0.05, significance := '']
  filtered_data[pval <= 0.05, significance := '*']
  filtered_data[pval <= 0.01, significance := '**']
  filtered_data[pval <= 0.001, significance := '***']
  
  
  # set order variable
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

  # Export and save the plot
  ggsave(filename = 'products/plot_yield_management_effect.png', plot = p1, width = 8,height = 6, units = "in", dpi = 300)

  #Create a barplot for the fertilisation doses and their estimates

  # Filter the data
  filtered_data <- out.dgr[grepl('^nutri_dose',factor) & !grepl('intercept|intrcpt', var)]

  # Add asterisks for significance
  filtered_data[pval > 0.05, significance := '']
  filtered_data[pval <= 0.05, significance := '*']
  filtered_data[pval <= 0.01, significance := '**']
  filtered_data[pval <= 0.001, significance := '***']
  
  
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
  ggsave(filename = "products/plot_yield_ndose_effect.png", plot = p2, width = 8, height = 6, units = "in",dpi = 300)

  #Create a barplot for the site properties and their estimates
  
  # Filter the data
  cols <- c('crop_residue', 'cover_crop', 'crop_rotation', 'tillage', 'fertilizer_type', 
            'bdod_mean_0_5', 'cec_mean_0_5', 'clay_mean_0_5', 'ntot_mean_0_5',
            'soc_mean_0_5', 'phw_mean_0_5', 'tmp_mean', 'pre_mean', 'pet_mean', 'ctype')
  
  filtered_data <- out.dgr[factor %in%  cols & !grepl('intercept|intrcpt', var)]

  # Add descriptive labels
  filtered_data[factor == 'crop_residue' & var == 'no',label := 'crop residue (no)']
  filtered_data[factor == 'crop_residue' & var == 'yes' ,label := 'crop residue (yes)']
  filtered_data[factor == 'crop_residue' & var == 'unknown' ,label := 'crop residue (unknown)']
  filtered_data[factor == 'cover_crop' & var == 'no' ,label := 'cover crop (no)']
  filtered_data[factor == 'cover_crop' & var == 'yes' ,label := 'cover crop (yes)']
  filtered_data[factor == 'cover_crop' & var == 'unknown' ,label := 'cover crop (unknown)']
  filtered_data[factor == 'crop_rotation' & var == 'no' ,label := 'crop rotation (no)']
  filtered_data[factor == 'crop_rotation' & var == 'yes' ,label := 'crop rotation (yes)']
  filtered_data[factor == 'crop_rotation' & var == 'unknown' ,label := 'crop rotation (unknown)']
  filtered_data[factor == 'tillage' & var == 'CT' ,label := 'tillage (CT)']
  filtered_data[factor == 'tillage' & var == 'NT' ,label := 'tillage (NT)']
  filtered_data[factor == 'tillage' & var == 'RT' ,label := 'tillage (RT)']
  filtered_data[factor == 'fertilizer_type' & var == 'combined' ,label := 'fertilizer type (combined)']
  filtered_data[factor == 'fertilizer_type' & var == 'inorganic' ,label := 'fertilizer type (inorganic)']
  filtered_data[factor == 'fertilizer_type' & var == 'organic' ,label := 'fertilizer type (organic)']
  filtered_data[factor == 'fertilizer_type' & var == 'unknown' ,label := 'fertilizer type (unknown)']
  filtered_data[factor == 'bdod_mean_0_5' ,label := 'bulk density']
  filtered_data[factor == 'cec_mean_0_5' ,label := 'cation exchange capacity']
  filtered_data[factor == 'clay_mean_0_5' ,label := 'clay content']
  filtered_data[factor == 'ntot_mean_0_5' ,label := 'total nitrogen']
  filtered_data[factor == 'soc_mean_0_5' ,label := 'soil organic carbon']
  filtered_data[factor == 'phw_mean_0_5' ,label := 'water pH']
  filtered_data[factor == 'tmp_mean' ,label := 'temperature mean']
  filtered_data[factor == 'pre_mean' ,label := 'precipitation mean']
  filtered_data[factor == 'pet_mean' ,label := 'evapotranspiration']
  filtered_data[factor == 'ctype' & var == 'arable' ,label := 'arable']
  filtered_data[factor == 'ctype' & var == 'beans_oilcrop' ,label := 'beans & oilcrop']
  filtered_data[factor == 'ctype' & var == 'cereal' ,label := 'cereals']
  filtered_data[factor == 'ctype' & var == 'maize' ,label := 'maize']
  filtered_data[factor == 'ctype' & var == 'other' ,label := 'other']
  filtered_data[factor == 'ctype' & var == 'rice' ,label := 'rice']
  filtered_data[factor == 'ctype' & var == 'vegetable' ,label := 'vegetables']
 
  # Add asterisks for significance
  filtered_data[pval > 0.05, significance := '']
  filtered_data[pval <= 0.05, significance := '*']
  filtered_data[pval <= 0.01, significance := '**']
  filtered_data[pval <= 0.001, significance := '***']
  
  
  # Define the order for the labels
  ordered_labels <- c('crop residue (no)', 'crop residue (yes)', 'crop residue (unknown)',
                    'cover crop (no)', 'cover crop (yes)', 'cover crop (unknown)',
                    'crop rotation (no)', 'crop rotation (yes)', 'crop rotation (unknown)',
                    'tillage (CT)', 'tillage (NT)', 'tillage (RT)',
                    'fertilizer type (combined)', 'fertilizer type (inorganic)', 'fertilizer type (organic)', 'fertilizer type (unknown)',
                   'bulk density', 'cation exchange capacity', 'clay content', 'total nitrogen',
                    'soil organic carbon', 'water pH', 'temperature mean', 'precipitation mean', 'evapotranspiration',
                    'arable', 'beans & oilcrop', 'cereals', 'maize', 'other', 'rice', 'vegetables'
                  )

  # Convert the label column to a factor with the specified levels
  filtered_data[,label := factor(label, levels = ordered_labels)]

  # Plot the filtered data with descriptive labels
  p12 <- ggplot(data = filtered_data) +
        geom_bar(aes(x = estimate, y = label), stat = "identity", color = 'lightblue', fill = 'lightblue') + 
        geom_errorbar(aes(y = label, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
        theme_bw() +
        geom_text(aes(x = estimate, y = label, label = significance), vjust = +0.4, hjust = -0.2, size = 4.5, color = "red") + 
        ggtitle('Effect of site properties on management induced changes in crop yield') + 
        labs(x = 'lnRR response of crop yield', y = '') +  # Update x and y axis titles
        theme(
          plot.title = element_text(size = 12, hjust = 0.5),  # Adjust title size and align to left
          axis.title.x = element_text(size = 10),           # Adjust x-axis title size
          axis.title.y = element_text(size = 12),           # Adjust y-axis title size
          axis.text.x = element_text(size = 8),             # Adjust x-axis text size
          axis.text.y = element_text(size = 9, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
        )
  # Export and save the plot
  ggsave(filename = "products/plot_yield_siteproperties_effect.png", plot = p12, width = 8, height = 6,units = "in",dpi = 300)
  
  }
  

# --- metaregression modelling ----
  
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

  #create new variable indicating the groups
  d2[grepl(' \\+ ',man_code),man_code_new := "Combination practices"]
  d2[man_code %in% c("NT", "ST", "RT"),man_code_new := "Tillage Practices"]
  d2[man_code %in% c("IR"),man_code_new := "Irrigation practices"]
  d2[man_code %in% c("DF", "NF", "MF", "OF"),man_code_new := "Fertilisation practices"]
  d2[man_code %in% c("RR"),man_code_new := "Residue practices"]
  d2[man_code %in% c("IMP"),man_code_new := "Biochar practice"]
  d2[man_code %in% c("CC", "INC", "MCR", "L"),man_code_new := "Crop practices"]
  d2[man_code %in%  c("ML"),man_code_new := "Mulching practices"]
  d2[is.na(man_code_new),man_code_new := "Drilling practices"]
  
  # a non linear response of crop yield due to n dose
  d2[,n_dose :=nutri_dose_N^2]

  # do you want to standardize (optional)
  if(FALSE){
    # Function to standardize numeric columns
    standardize <- function(x) {return ((x - mean(x)) / sd(x))}
    
    # standardize numerical columns
    cols <- c('n_dose', 'nutri_dose_P', 'nutri_dose_K','nutri_dose_Mg', 'ntot_mean_0_5', 'tmp_mean', 
              'pre_mean', 'pet_mean', 'soc_mean_0_5', 'clay_mean_0_5', 'bdod_mean_0_5', 'cec_mean_0_5', 'phw_mean_0_5' )
    d2[,c(cols) := lapply(.SD,standardize),.SDcols=cols]
  }
  

  # Fit the model using the new grouping variable
  m3.full <- rma.mv(yi, vi, 
                  mods = ~ man_code_new + ctype + crop_residue + cover_crop + crop_rotation + tillage + 
                    fertilizer_type + n_dose + nutri_dose_P + nutri_dose_K +
                    nutri_dose_Mg + ntot_mean_0_5 + tmp_mean + pre_mean + pet_mean +  soc_mean_0_5 + 
                    clay_mean_0_5 + bdod_mean_0_5 + cec_mean_0_5 + phw_mean_0_5 -1,
                  data = d2, 
                  random = list(~ 1 | study_ID), 
                  method = "REML", 
                  sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)

# --- plotting meta-regression results ----
  
  if(FALSE){
    
  ##Barplot of practices & site properties for the full model

  # Extract and tidy the results
  m3_tidy <- as.data.table(broom::tidy(m3.full))

  # Define new column values
  new_column_values <- c(rep("management_practice",9),rep("crop_type",6),"crop_residue", "crop_residue", 
                         "cover_crop", "cover_crop", "crop_rotation", 
                         "crop_rotation", "tillage", "tillage", rep("fertilizer_type ",3),
                         rep("nutri_dose", 4), "soil_property", 
                         "soil_property", rep("weather_condition", 3), rep("soil_property", 4))

  m3_tidy[,colorfactor := new_column_values]

  # Rename factor levels using recode with named arguments
  m3_tidy[grepl('man_code_new',term),factor := gsub('^man_code_new| practice| [P|p]ractices','',term)]
  m3_tidy[grepl('man_code_newOther',term),factor := 'Other practices']
  m3_tidy[grepl('^ctype',term) & grepl('oil',term), factor := 'Beans & oilcrop']
  m3_tidy[grepl('^ctype',term) & grepl('cereal',term), factor := 'Cereals']
  m3_tidy[grepl('^ctype',term) & grepl('maize',term), factor := 'Maize']
  m3_tidy[grepl('^ctype',term) & grepl('other',term), factor := 'Other']
  m3_tidy[grepl('^ctype',term) & grepl('rice',term), factor := 'Rice']
  m3_tidy[grepl('^ctype',term) & grepl('vege',term), factor := 'Vegetables']
  m3_tidy[grepl('^crop|^cover',term),factor := gsub('_',' ',term)]
  m3_tidy[grepl('^crop|^cover',term),factor := gsub('yes$','',factor)]
  m3_tidy[grepl('^crop|^cover',term),factor := gsub('unknown',' (U)',factor)]
  m3_tidy[grepl('^nutri',term), factor := paste0(gsub('nutri_dose_','',term),' dose')]
  m3_tidy[term == "tillageNT", factor := "No tillage"]
  m3_tidy[term == "tillageRT", factor := "Reduced tillage"] 
  m3_tidy[term == "fertilizer_typeinorganic", factor := "Inorganic fertiliser"]
  m3_tidy[term == "fertilizer_typeorganic", factor := "Organic fertiliser"]
  m3_tidy[term == "fertilizer_typeunknown", factor := "Fertiliser (U)"]
  m3_tidy[term == "n_dose" , factor := "N dose"]
  m3_tidy[term == "ntot_mean_0_5" , factor := "Total nitrogen"]
  m3_tidy[term == "tmp_mean" , factor := "Temp"]
  m3_tidy[term == "pre_mean" , factor := "Precipitation"]
  m3_tidy[term == "pet_mean" , factor := "Evapotranspiration"]
  m3_tidy[term == "soc_mean_0_5" , factor := "Soil organic C"]
  m3_tidy[term == "clay_mean_0_5" , factor := "Clay"]
  m3_tidy[term == "bdod_mean_0_5" , factor := "Bulk density"]
  m3_tidy[term == "cec_mean_0_5" , factor := "Cation exchange capacity"]
  m3_tidy[term == "phw_mean_0_5" , factor := "pH Water"]
  
# Add a column for significance stars
  m3_tidy[p.value > 0.1, significance := '']
  m3_tidy[p.value <= 0.1, significance := '.']
  m3_tidy[p.value <= 0.05, significance := '*']
  m3_tidy[p.value <= 0.01, significance := '**']
  m3_tidy[p.value <= 0.001, significance := '***']
  
  # Define the order of the terms
  m3_tidy[,factor := factor(factor, levels = unique(factor))]

  # Plot for site properties with SE
  p5 <- ggplot(m3_tidy, aes(x = factor, y = estimate)) + 
        geom_bar(stat = "identity", aes(fill = colorfactor), alpha = 0.7) + 
        geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                      width = 0.2, position = position_dodge(0.9)) + 
        geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
                  vjust = ifelse(m3_tidy$estimate > 0, -0.5, 1.5), size = 4) + 
        theme_minimal() + 
        labs(x = 'Management practices and site properties', y = 'Parameter estimate', 
             title = 'Meta-regression model estimation on crop yield') + 
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1, size = 7.6), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
          axis.ticks.length = unit(-0.5, "cm"), 
          panel.grid = element_blank(), 
          legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA),  # Set plot background to white
          plot.title = element_text(size = 14, hjust = 0.5)  # Center the title and adjust the size
        ) + 
        ylim(-1, 1.1)

  # Save the plot
  ggsave(filename = "products/plot_yield_metamodel1_effect.png", plot = p5, width = 8, height = 6, units = "in", dpi = 300)
  
  }
  
# --- final regression model ----
  
  ##different grouping and final refined model

  #create new variable indicating the groups
  d2[grepl(' \\+ ',man_code),man_code_new2 := "Combination practices"]
  d2[man_code %in% c("NT", "ST", "RT","ML", "IMP","RR","DD"),man_code_new2 := "Soil practices"]
  d2[man_code %in% c("DF", "NF", "MF", "OF"),man_code_new2 := "Fertilisation practices"]
  d2[man_code %in% c("IR"),man_code_new2 := "Irrigation practices"]
  d2[man_code %in% c("CC", "INC", "MCR"),man_code_new2 := "Crop practices"]
  d2[is.na(man_code_new),man_code_new2 := "Other"]

  
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

  # Estats of all models 
  estats(m3.empty,m3.full1)
  estats(m3.full,m3.full1)
  
  
  # Empty model compared to the different full models
  anova(m3.empty,m3.full,refit = T)
  anova(m3.empty,m3.full1,refit = T)
  
  # Comparison of first full model with the refined model
  anova(m3.full,m3.full1,refit = T)
  
# --- plots final regression model ----
  
  if(FALSE){
  
  # Extract and tidy the results
  m3_tidy <- as.data.table(broom::tidy(m3.full1))

  # Define new column values
  new_column_values1 <- c(rep("management_practice",5),rep("crop_type",6),
                         "cover_crop", "cover_crop", "crop_rotation", 
                         "crop_rotation", rep("soil_property", 5), "nutrient_dose", "weather_condition")
  
  m3_tidy[,colorfactor1 := new_column_values1]


  # Rename terms to more readable names
  m3_tidy[grepl('man_code_new',term),factor := gsub('^man_code_new2| practice| [P|p]ractices','',term)]
  m3_tidy[grepl('man_code_newOther',term),factor := 'Other practices']
  m3_tidy[grepl('^ctype',term) & grepl('oil',term), factor := 'Beans & oilcrop']
  m3_tidy[grepl('^ctype',term) & grepl('cereal',term), factor := 'Cereals']
  m3_tidy[grepl('^ctype',term) & grepl('maize',term), factor := 'Maize']
  m3_tidy[grepl('^ctype',term) & grepl('other',term), factor := 'Other']
  m3_tidy[grepl('^ctype',term) & grepl('rice',term), factor := 'Rice']
  m3_tidy[grepl('^ctype',term) & grepl('vege',term), factor := 'Vegetables']
  m3_tidy[grepl('^crop|^cover',term),factor := gsub('_',' ',term)]
  m3_tidy[grepl('^crop|^cover',term),factor := gsub('yes$','',factor)]
  m3_tidy[grepl('^crop|^cover',term),factor := gsub('unknown',' (U)',factor)]
  m3_tidy[grepl('^nutri',term), factor := paste0(gsub('nutri_dose_','',term),' dose')]
  m3_tidy[term == "tillageNT", factor := "No tillage"]
  m3_tidy[term == "tillageRT", factor := "Reduced tillage"] 
  m3_tidy[term == "fertilizer_typeinorganic", factor := "Inorganic fertiliser"]
  m3_tidy[term == "fertilizer_typeorganic", factor := "Organic fertiliser"]
  m3_tidy[term == "fertilizer_typeunknown", factor := "Fertiliser (U)"]
  m3_tidy[term == "n_dose" , factor := "N dose"]
  m3_tidy[term == "ntot_mean_0_5" , factor := "Total nitrogen"]
  m3_tidy[term == "tmp_mean" , factor := "Temp"]
  m3_tidy[term == "pre_mean" , factor := "Precipitation"]
  m3_tidy[term == "pet_mean" , factor := "Evapotranspiration"]
  m3_tidy[term == "soc_mean_0_5" , factor := "Soil organic C"]
  m3_tidy[term == "clay_mean_0_5" , factor := "Clay"]
  m3_tidy[term == "bdod_mean_0_5" , factor := "Bulk density"]
  m3_tidy[term == "cec_mean_0_5" , factor := "Cation exchange capacity"]
  m3_tidy[term == "phw_mean_0_5" , factor := "pH Water"]


  # Add a column for significance stars
  m3_tidy[p.value > 0.1, significance := '']
  m3_tidy[p.value <= 0.1, significance := '.']
  m3_tidy[p.value <= 0.05, significance := '*']
  m3_tidy[p.value <= 0.01, significance := '**']
  m3_tidy[p.value <= 0.001, significance := '***']
  
  # Define the order of the terms
  m3_tidy[,factor := factor(factor, levels = unique(factor))]

  p7 <- ggplot(m3_tidy, aes(x = factor, y = estimate)) + 
        geom_bar(stat = "identity", aes(fill = colorfactor1), alpha = 0.7) + 
        geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                      width = 0.2, position = position_dodge(0.9)) + 
        geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + std.error + 0.02, estimate - std.error - 0.02)), 
                  vjust = ifelse(m3_tidy$estimate > 0, -0.5, 1.5), size = 4) + 
        theme_minimal() + 
        labs(x = 'Management practices and site properties', y = 'Parameter estimate', 
             title = 'Meta-regression model estimation on crop yield') + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), 
              panel.border = element_rect(color = "black", fill = NA, size = 1), 
              axis.ticks.length = unit(-0.5, "cm"), 
              panel.grid = element_blank(), 
              legend.position = "none",
              plot.title = element_text(size = 14, hjust = 0.5, color = "black"),  # Adjust title size, center, and color
              plot.background = element_rect(fill = "white", color = NA)) +  # Set plot background to white
        ylim(-0.8, 0.8)
  
  # Save the plot
  ggsave(filename = 'products/plot_yield_metamodel2_effect.png',plot = p7, width = 8,height = 6, units = "in",dpi = 300)

  }


# --- application of the model on the INTEGRATOR dataset ----

  # rerun the above code to regenerate model m3.full1 but now without standardization
  # save the model and use this model for application
  saveRDS(m3.full1,file='products/mamodel_yield.rds')

  # read in the model
  m3.full1 <- readRDS('products/mamodel_yield.rds')
  
  # load external datasets not stored at github given its size
  floc <- 'D:/DATA/17 nutribudget/'
  d4 <- fread(paste0(floc,'db_final_europe.csv'))

  # load in covariates
  d4.cov <- fread('products/240723_covariates_ncu.csv')

  # merge both files(note that not all NCUs are in covariates, to be checked later)
  d4 <- merge(d4,d4.cov,by.x='ncu',by.y='gncu2010_ext',all.x=TRUE)
  
  # predict the yield response (lnRR) for each crop and site
  # note that this is done manually (from reading summary(m3.full1) to avoid all re-naming of variables)
  
  # model coefficients
  m3.coeff <- as.data.table(broom::tidy(m3.full1))
  
  # replacing missing inputs with median value
  d4[,parea.cr := pmax(0,parea.cr,na.rm=T)]
  d4[,parea.cc := pmax(0,parea.cc,na.rm=T)]
  d4[is.na(ntot_mean_0_5), ntot_mean_0_5 := median(d4$ntot_mean_0_5,na.rm=T)]
  d4[is.na(phw_mean_0_5), phw_mean_0_5 := median(d4$phw_mean_0_5,na.rm=T)]
  d4[is.na(cec_mean_0_5), cec_mean_0_5 := median(d4$cec_mean_0_5,na.rm=T)]
  d4[is.na(bdod_mean_0_5), bdod_mean_0_5 := median(d4$bdod_mean_0_5,na.rm=T)]
  d4[is.na(clay_mean_0_5), clay_mean_0_5 := median(d4$clay_mean_0_5,na.rm=T)]
  d4[is.na(pre_mean), pre_mean := median(d4$pre_mean,na.rm=T)]
  
  # estimate the baseline response given the soil properties
  d4[, lnRR := m3.coeff[grepl('pre',term),estimate] * pre_mean]
  d4[, lnRR := lnRR +  m3.coeff[grepl('nutri_dose_P',term),estimate] * mean(d1$nutri_dose_P)]
  d4[, lnRR := lnRR +  m3.coeff[grepl('ntot_mean',term),estimate] * ntot_mean_0_5]
  d4[, lnRR := lnRR +  m3.coeff[grepl('phw_mean',term),estimate] * phw_mean_0_5]
  d4[, lnRR := lnRR +  m3.coeff[grepl('cec_mean',term),estimate] * cec_mean_0_5]
  d4[, lnRR := lnRR +  m3.coeff[grepl('bdod_mean',term),estimate] * bdod_mean_0_5]
  d4[, lnRR := lnRR +  m3.coeff[grepl('clay_mean',term),estimate] * clay_mean_0_5]
  # for measures, if not applicable, select the "unknown" option of the model (might also be the "no")
  d4[, lnRR := lnRR + (parea.cr * -0.0851764663 + (area_ncu_ha - parea.cr) * 0.1750006009)/area_ncu_ha]
  d4[, lnRR := lnRR + (parea.cc *  0.2708926575 + (area_ncu_ha - parea.cc) * 0.1075491926)/area_ncu_ha]
  d4[crop_name == 'Rice', lnRR := lnRR + m3.coeff[grepl('ctyperice',term),estimate]]
  d4[grepl('wheat|barly|oats|cerea',tolower(crop_name)), lnRR := lnRR + m3.coeff[grepl('ctypecereal',term),estimate]]
  d4[grepl('puls|soya|rape|fodder|oil|rey',tolower(crop_name)), lnRR := lnRR + m3.coeff[grepl('ctypebeans_oilcrop',term),estimate]]
  d4[grepl('maize',tolower(crop_name)), lnRR := lnRR + m3.coeff[grepl('ctypemaize',term),estimate]]
  d4[grepl('vegetabl',tolower(crop_name)), lnRR := lnRR + m3.coeff[grepl('ctypevegetable',term),estimate]]
  d4[grepl('fruits|olive|fibre|wine|tobac|olive|industr',tolower(crop_name)), lnRR := lnRR + m3.coeff[grepl('ctypeother',term),estimate]]
  d4[grepl('potat|sugar beet',tolower(crop_name)), lnRR := lnRR + 0]
  d4[,mlnRR_soil := exp(lnRR + m3.coeff[grepl('man_code_new2Soil',term),estimate])-1]
  d4[,mlnRR_irr := exp(lnRR + m3.coeff[grepl('man_code_new2Irr',term),estimate])-1]
  d4[,mlnRR_fert := exp(lnRR + m3.coeff[grepl('man_code_new2Fert',term),estimate])-1]
  d4[,mlnRR_crop := exp(lnRR + m3.coeff[grepl('man_code_new2Crop',term),estimate])-1]
  d4[,mlnRR_combi := exp(lnRR + m3.coeff[grepl('man_code_new2Com',term),estimate])-1]
  
  # take weighted mean per ncu
  cols <- colnames(d4)[grepl('^mlnRR|^yield_',colnames(d4))]
  
  d5 <- d4[,lapply(.SD,function(x) weighted.mean(x,w=area_ncu)),.SDcols = cols,by='ncu']
  
# --- plot impact of measures on crop yield ----
  
  # laod in NUTS shapefile
  s.nuts <- st_read(paste0(floc,'eu_nuts.gpkg'),layer='eu_nuts')
  
  # get the raster to plot
  r1 <- terra::rast(paste0(floc,'gncu2010_ext.asc'))
  terra::crs(r1) <- 'epsg:3035'
  r1 <- terra::project(r1,'epsg:4326',method='near')
  
  # convert to data.frame
  r1.p <- as.data.frame(r1,xy=TRUE)
  r1.p <- as.data.table(r1.p)
  
  # join/merge d1 with r1.p
  r.ncu <- merge(r1.p, d5, by.x = 'gncu2010_ext', by.y = 'ncu')

  # plot impact of soil measures on crop yield
  pbreak <- c(-1000,0,20,40,60,80,1000)
  plabel <- c('<0','0-20','20-40','40-60','60-80','>80')
  p1 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(mlnRR_soil*100, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nyield (%)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of soil measures on crop yield") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p1, filename = paste0(floc,'up_yield_meas_soil.jpg'),width = 12,height=12,units='cm')  
  
  p2 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(mlnRR_fert *100, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nyield (%)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of fertilization (4R) on crop yield") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p2, filename = paste0(floc,'up_yield_meas_fert.jpg'),width = 12,height=12,units='cm') 
  
  p3 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(mlnRR_crop *100, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nyield (%)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of crop measures on crop yield") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p3, filename = paste0(floc,'up_yield_meas_crop.jpg'),width = 12,height=12,units='cm') 
  
  p4 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(mlnRR_combi *100, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nyield (%)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of combined measures on crop yield") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p4, filename = paste0(floc,'up_yield_meas_combi.jpg'),width = 12,height=12,units='cm') 
  
  p5 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= fifelse(yield_target <= yield_ref,'realised','not realised'))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Target yield\n achieved')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Crop target yield achieved, baseline") +
      coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p5, filename = paste0(floc,'up_yield_target_bau.jpg'),width = 12,height=12,units='cm') 
  
  p6 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= fifelse(yield_target <= yield_ref * (1+mlnRR_combi),'realised','not realised'))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Target yield\nachieved')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Crop target yield achieved, combi of measures") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p6, filename = paste0(floc,'up_yield_target_combi.jpg'),width = 12,height=12,units='cm') 
  
  
  data.table(quantile = paste0('Q',c(0.05,0.25,0.50,0.75,0.95)),
             soil = round(quantile(d5$mlnRR_soil*100,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             fert = round(quantile(d5$mlnRR_fert*100,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             crop = round(quantile(d5$mlnRR_crop   *100,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             combi = round(quantile(d5$mlnRR_combi*100,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2))
  
  # save predicted output
  d6 <- copy(d5)
  setnames(d6,tolower(gsub('mlnRR_','yield_',colnames(d6))))
  saveRDS(d6,paste0(floc,'ncu_yield_meas.rds'))
  