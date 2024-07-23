# data analysis RQ on pH
# Salim Belyazid & Gerard Ros
# 23 May 2024

#required packages
require(readxl);require(data.table);require(metafor);require(ggplot2); require(sf);require(terra)

# clear environment
rm(list=ls())

# --- part 1. loading and preparing datasets ----

  # read in the database
  d1 <- as.data.table(read_xlsx('data/ph_database.xlsx',sheet='database'))
  row.names(d1)<-NULL
  
  # update column names to simply internal references to the column names (so no spaces or brackets)
  # and set all column names to lower case
  #setnames(d1,old = c('CEC (meq/100g)'),new = c('cec'))
  setnames(d1,tolower(colnames(d1)))
  
  # remove text from latitude column, and make the column numeric (text is aumatically converted to NA)
  d1[,lat := as.numeric(lat)]
  
  # save the unique lon-lat as csv for covariate extraction
  # the covariates are separately collected in script covariate extraction
  # this you only need to do once (run the code within the FALSE statement)
  # then a csv file is saved, that will be read in al subseqent times you run this script
  if(FALSE){
    # subset only relevant columns
    d1.lonlat <- d1[,.(lat,lon)]
    # select only the cases where latitude is given
    d1.lonlat <- unique(d1.lonlat[!is.na(lat)])
    fwrite(d1.lonlat,'data/ph_lonlat_unique.csv')
    
  } else {
    
    # after collection site properties, read these covariates in
    d1.cov <- fread('products/240617_covariates_ph.csv')
  }
  
  # combine the extracted coariates with the original dataset
  d1 <- merge(d1,d1.cov,by = c('lat','lon'),all.x=TRUE)
  
  # set covariates to mean value when lon-lat is missing (action Salim: better to take most nearby proxy lon-lat)
  cols <- colnames(d1)[grepl('_mean|_sd',colnames(d1))]
  d1[,c(cols) := lapply(.SD,function(x) fifelse(is.na(x),median(x,na.rm=T),x)),.SDcols = cols]

  # estimate uncertainty when SE and SD are missing

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
  d1[,c('lat','lon','reference','GEnZname','year') := NULL]

  # use only topsoil properties ISRIC and remove subsoil properties
  cols <- colnames(d1)[grepl('_5_15$|15_30$',colnames(d1))]
  d1[,c(cols):=NULL]
  
  # add unique rowid
  d1[,id := 1:.N]
  
  # helper function to convert classes to the mean value
  hfunconvert <- function(x) {mean(as.numeric(unlist(strsplit(gsub('<|>','',x),'-|~'))))}
  
  # update dataset (MA check that not all data entries are numeric)
  cols <- c('duration','mat','map','soc','ph','clay')

  # apply the function to expected numeric variables  
  d1[, c(cols) := lapply(.SD,function(x) hfunconvert(x)),.SDcols = cols,by=id]
  d1[, c(cols) := lapply(.SD,function(x) as.numeric(x)),.SDcols = cols]

  # replace missing properties
  d1[is.na(ph), ph := median(d1$ph,na.rm=T)]
  d1[is.na(duration), duration := median(d1$duration,na.rm=T)]
  d1[is.na(mat), mat := median(d1$mat,na.rm=T)]
  d1[is.na(map), map := median(d1$map,na.rm=T)]
  d1[is.na(soc), soc := median(d1$soc,na.rm=T)]
  d1[is.na(clay), clay := median(d1$clay,na.rm=T)]

  # replace missing replication
  d1[is.na(replication), replication := 2]

# --- part 2. plot the datasets based on the covariance -----

  # make plot sampling location
  if(FALSE){
    
  # load in the base map
  world <- map_data("world")
  ggplot() +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "#999999", fill = "#CCCCCC", linewidth = 0.1) +
    geom_point(data = d1.cov,aes(lon, lat), alpha = 1, size = 2,shape = 21, 
               color = "red", fill = "black", stroke = 0.5) +
    theme_bw()+
    ylim(30,70)+
    xlim(-20,40)+
    coord_sf(crs=4326)
  }

# --- part 3. Analysis for KPI ph -----

  # remove observations without a control or treatment value
  d2 <- d1[!(is.na(kpi_treat) | is.na(kpi_control))]
  
  # check the data.table on missing values
  # note that bs adn cec has too limited numbers, so not use them
  summary(d2)

  # add standardized mean difference response ratio as effect size
  # other options are also possible, change measure argument.
  # see ?escalc for argument description and examples
  d2 <- escalc(measure = "SMD", data = d2, 
               m1i = kpi_treat , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control, sd2i= kpi_contr_sd, n2i = replication)
  
  # convert back to data.table
  d2 <- as.data.table(d2)

  # remove items with missing yi or vi
  d2 <- d2[!is.na(yi) | !is.na(vi)]
  
  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]

  # make first plot of individual observations
  if(FALSE){
    
    ggplot(data = d2,aes(x= id,y=yi)) + 
      geom_point() + 
      geom_line()+
      geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
      ylim(-5,5) +
      theme_bw() + xlab('study id') + ylab('SMD, change in soil pH') + 
      ggtitle('Observed changes in soil pH due to management practices')
  }
  
  # given the huge changes and deline, a standardized change bigger than 5 is unlikely
  d2 <- d2[abs(yi)<=5]

  # estimate mean response across all studies
  m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_id), 
                        method="REML",
                        sparse = TRUE)

  # see the summary of the model
  summary(m1)


  #1- management
  #fcols <- c('man_code')
             
  #2- climate
  #fcols <- c('map','mat','pet_mean','ph','soc','clay','ph','soc','clay','bdod_mean_0_5','cec_mean_0_5')
  
  #3- site conditions
  fcols <- c('crop_type','crop_rotation',
             'crop_residue','cover_crop','tillage',
             'fertilizer_type')

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
    m2 <- metafor::rma.mv(yi,vi, mods = ~factor(evar)-1, data=d2,  #-1 means we dont produce intercepts, without -1 the baseline is the avrage of all treatments 
                          random= list(~ 1|study_id), method="REML",sparse = TRUE)
    
    # collect model coefficients
    m2.sum <- summary(m2)
    m2.out <- as.data.table(coefficients(m2.sum))
    m2.out[,factor := i]
    m2.out[,var := gsub('factor\\(evar\\)','',rownames(m2.sum$b))]
    
  } else{
    
    # do the meta-analysis for a numerical variable
    m2 <- metafor::rma.mv(yi,vi, mods = ~evar, intercept=FALSE, data=d2,
                          random= list(~ 1|study_id), method="REML",sparse = TRUE)
    
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
  
  
  #extract the significance of the estimate
  out.dgr[pval <= 0.001, significance := '***']
  out.dgr[pval <= 0.01, significance := '**']
  out.dgr[pval <= 0.05, significance := '*']
  out.dgr[pval <= 0.10, significance := '.']
  out.dgr[pval > 0.10, significance := '']

  ## plot the var and factor from the model
  if(FALSE){
    
  p1 <- ggplot(data = out.dgr[var!='intercept']) +
        geom_bar(aes(x = estimate, y = var), stat = "identity", color='lightblue', fill = 'lightblue') + 
        geom_errorbar(aes(y = var, xmin = estimate - se, xmax = estimate + se), width = 0.4) + 
        theme_bw() +
        geom_text(aes(x = estimate, y = var, label = significance), vjust = +0.4, hjust = -0.2, size = 6, color = "red")  +                     #Add the asterisks for significance
      #  xlab('Effect estimate') + ylab('Experimental Management') +
      #  xlab('Effect estimate') + ylab('Climatic & soil conditions') +
        xlab('Effect estimate') + ylab('Site management') +
        theme(
          plot.title = element_text(size = 10, hjust = 0),  # Adjust title size and align to left
          axis.title.x = element_text(size = 12),         # Adjust x-axis title size
          axis.title.y = element_text(size = 12),         # Adjust y-axis title size
          axis.text.x = element_text(size = 12),           # Adjust x-axis text size
          axis.text.y = element_text(size = 12, margin = margin(t = 20, b = 10))  # Adjust y-axis text size
        )
  
  # Export and save the plot
  ggsave(filename = 'products/plot_ph_siteman_effect.png', plot = p1, width = 8,height = 6, units = "in", dpi = 300)
  }
  

# ---- part 4. do one meta-regression model with multiple factors -----

  # make a function to extract relevant model statistics
  estats <- function(model_new,model_base){
    out <- data.table(AIC = model_new$fit.stats[4,2],
                      ll = model_new$fit.stats[1,2],
                      ll_impr = round(100 * (1-model_new$fit.stats[1,2]/model_base$fit.stats[1,2]),2),
                      r2_impr = round(100*max(0,(sum(model_base$sigma2)-sum(model_new$sigma2))/sum(model_base$sigma2)),2),
                      pval = round(anova(model_new,model_base)$pval,3))
    return(out)
  }

  # Function to standardize numeric columns
  standardize <- function(x) {return ((x - mean(x)) / sd(x))}

  # standardize the numerical variables
  cols <- c('ph', 'clay', 'soc', 'cec_mean_0_5', 'bdod_mean_0_5', 
            'mat', 'map', 'pet_mean', 'crop_rotation')
  d2[,c(cols) := lapply(.SD,standardize),.SDcols=cols]

  # make first an empty model
  m3.empty <- metafor::rma.mv(yi,vi, data=d2,random= list(~ 1|study_id), method="REML",sparse = TRUE)

  # show summary
  summary(m3.empty)

  # make a full model with all main factors together
  # add here all the variables (as given in out.dgr) that are signficant. 
  # Insignifcant variables can be added but only if you think there is a mechanistic reason and you like to show/proof that there is limited/no effect
  # note that the total number of variables can not be too high given the limited number of observations.
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~man_code +
                               ph + clay + soc + cec_mean_0_5 + bdod_mean_0_5 + 
                               mat + map + pet_mean +
                               tillage + fertilizer_type + 
                               crop_type + crop_residue + cover_crop + crop_rotation-1,
                             data=d2,random= list(~ 1|study_id), method="REML",sparse = TRUE,control=list(rel.tol=1e-8))

  # analyse summary stats
  summary(m3.full)


  ##Barplot of practices & site properties for the full model

  if(FALSE){
    
    # Extract and tidy the results
    m3_tidy <- as.data.table(broom::tidy(m3.full))
    
    # Define new column values
    new_column_values <- c("intercept","management_practice","management_practice", "management_practice", "management_practice",
                           "management_practice","management_practice","management_practice","management_practice",
                           "soil_property","soil_property","soil_property","soil_property","soil_property",
                           "climatic_property","climatic_property","climatic_property",
                           "tillage","tillage","tillage",
                           "fertilizer_type","fertilizer_type","fertilizer_type",
                           "crop_type","crop_type","crop_type","crop_type","crop_type",
                           "crop_residue", 
                           "cover_crop", "cover_crop", "cover_crop",
                           "crop_rotation")
    
    m3_tidy[,colorfactor := new_column_values]
    
    # add significance
    m3_tidy[p.value <= 0.001, significance := '***']
    m3_tidy[p.value <= 0.01, significance := '**']
    m3_tidy[p.value <= 0.05, significance := '*']
    m3_tidy[p.value <= 0.10, significance := '.']
    m3_tidy[p.value > 0.10, significance := '']
    
    # add new column names
    m3_tidy[term == 'man_codeBC',factor := 'Biochar']
    m3_tidy[term == 'man_codeBC',factor:= "Biochar"]
    m3_tidy[term == 'man_codeCC',factor:= "Continuous cover"]  
    m3_tidy[term == 'man_codeCT',factor:= "Conventional tillage"]
    m3_tidy[term == 'man_codeLM',factor:= "Liming"]
    m3_tidy[term == 'man_codeNT',factor:= "No tillage Mngt"]
    m3_tidy[term == 'man_codeNF',factor:= "Mineral fertilization"]
    m3_tidy[term == 'man_codeOF',factor:= "Organic fertilization"]
    m3_tidy[term == 'man_codeRT',factor:= "Reduced tillage"]
    m3_tidy[term == 'man_codeZT',factor:= "Zonation tillage"]
    m3_tidy[term == 'ph',factor:= "Soil pH"]
    m3_tidy[term == 'clay',factor:= "Clay content"]
    m3_tidy[term == 'soc',factor:= "Soil organic carbon"]
    m3_tidy[term == 'cec_mean_0_5',factor:= "Cation exchnge capacity"]
    m3_tidy[term == 'bdod_mean_0_5',factor:= "Soil bulk density"]
    m3_tidy[term == 'mat',factor:= "Mean annual temperature"]
    m3_tidy[term == 'map',factor:= "Mean annual precipitation"]
    m3_tidy[term == 'pet_mean',factor:= "Mean annual evapotranspiration"]
    m3_tidy[term == 'tillageno tillage',factor:= "No tillage Cntrl"]
    m3_tidy[term == 'tillagereduced tillage',factor:= "Low tillage"]
    m3_tidy[term == 'tillagetillage unknown',factor:= "Unknown tillage"]
    m3_tidy[term == 'fertilizer_typeMineral Fertiliser',factor:= "Mineral fertilizer"]
    m3_tidy[term == 'fertilizer_typeNo fertilizer',factor:= "No fertilizer"]
    m3_tidy[term == 'fertilizer_typeOrganic fertiliser',factor:= "Organic fertilizer"]
    m3_tidy[term == 'crop_typefruit',factor:= "Fruit crop"]
    m3_tidy[term == 'crop_typegrass',factor:= "Grass crop"]
    m3_tidy[term == 'crop_typelegume',factor:= "Legume crop"]
    m3_tidy[term == 'crop_typemaize',factor:= "Maize crop"]
    m3_tidy[term == 'crop_typeroot crop',factor:= "Root crop"]
    m3_tidy[term == 'crop_residueResidue retained',factor:= "Residue retained"]
    m3_tidy[term == 'cover_cropContinuous cropping',factor:= "Continuous cropping"]
    m3_tidy[term == 'cover_cropCover crop unknow',factor:= "Cover crop unknown"]
    m3_tidy[term == 'cover_cropNo cover crop',factor:= "No cover crop"]
    m3_tidy[term == 'crop_rotation',factor:= "Crop rotation"]
    
    # Define the order of the terms
    m3_tidy[,factor := factor(factor, levels = unique(factor))]
    
    # split datasets into management practices and site properties
    dt.mp <- m3_tidy[grepl('^man',colorfactor)]
    dt.sp <- m3_tidy[!grepl('^man',colorfactor)]
    
    # Plot for management practices
    p4 <- ggplot(dt.mp, aes(x = factor, y = estimate)) + 
          geom_bar(stat = "identity", aes(fill = colorfactor), alpha = 0.7, fill = "lightblue") + 
          geom_errorbar(aes(y = estimate, ymin = estimate - std.error, ymax = estimate + std.error), width = 0.1) +
          geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + 0.001, estimate - 0.001)), 
                    vjust = ifelse(dt.mp$estimate > 0, -0.5, 1.5), size = 6) + 
          theme_minimal() + 
          labs(x = 'Management practices', y = 'Parameter estimate') + 
          theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
                axis.ticks.length = unit(-0.5, "cm"), 
                panel.grid = element_blank(), 
                legend.position = "none") + 
          ylim(-15, 15)
    
    ggsave(plot = p4, filename = 'products/plot_ph_management_effect.png', width = 8,height = 6, units = "in", dpi = 300)
  
    # Plot for site properties
    p5 <- ggplot(dt.sp, aes(x = factor, y = estimate)) + 
          geom_bar(stat = "identity", aes(fill = colorfactor), alpha = 0.7, fill = "lightblue") + 
          geom_errorbar(aes(y = estimate, ymin = estimate - std.error, ymax = estimate + std.error), width = 0.1) +
          geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + 0.001, estimate - 0.001)), 
                    vjust = ifelse(dt.sp$estimate > 0, -0.5, 1.5), size = 4) + 
          theme_minimal() + 
          labs(x = 'Site properties', y = 'Parameter estimate') + 
          theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
                axis.ticks.length = unit(-0.5, "cm"), 
                panel.grid = element_blank(), 
                legend.position = "none") + 
          ylim(-15, 15)
      
    ggsave(plot = p5, filename = 'products/plot_ph_siteprop_effect.png', width = 8,height = 6, units = "in", dpi = 300)
    
  }
  
# --- part 5. make full model ----

  # metaregression
  m4.full <- metafor::rma.mv(yi,vi, 
                           mods = ~man_code +  ph  + mat + tillage + crop_type + fertilizer_type -1,
                           data=d2,random= list(~ 1|study_id), method="REML",sparse = TRUE,control=list(rel.tol=1e-8))
  # add summary stats
  summary(m4.full)

  # Extract and tidy the results
  m4_tidy <- as.data.table(broom::tidy(m4.full))

  # make plot
  if(FALSE){
    
    # Define new column values
    m4_tidy[term == 'man_codeBC',factor := 'Biochar']
    m4_tidy[term == 'man_codeCC',factor:= "Continuous cover"]  
    m4_tidy[term == 'man_codeCT',factor:= "Conventional tillage"]
    m4_tidy[term == 'man_codeLM',factor:= "Liming"]
    m4_tidy[term == 'man_codeNT',factor:= "No tillage Mngt"]
    m4_tidy[term == 'man_codeNF',factor:= "Mineral fertilization"]
    m4_tidy[term == 'man_codeOF',factor:= "Organic fertilization"]
    m4_tidy[term == 'man_codeRT',factor:= "Reduced tillage"]
    m4_tidy[term == 'man_codeZT',factor:= "Zonation tillage"]
    m4_tidy[term == 'ph',factor:= "Soil pH"]
    m4_tidy[term == 'mat',factor:= "Mean annual temperature"]
    m4_tidy[term == 'tillageno tillage',factor:= "No tillage Cntrl"]
    m4_tidy[term == 'tillagereduced tillage',factor:= "Low tillage"]
    m4_tidy[term == 'tillagetillage unknown',factor:= "Unknown tillage"]
    m4_tidy[term == 'fertilizer_typeMineral Fertiliser',factor:= "Mineral fertilizer"]
    m4_tidy[term == 'fertilizer_typeNo fertilizer',factor:= "No fertilizer"]
    m4_tidy[term == 'fertilizer_typeOrganic fertiliser',factor:= "Organic fertilizer"]
    m4_tidy[term == 'crop_typefruit',factor:= "Fruit crop"]
    m4_tidy[term == 'crop_typegrass',factor:= "Grass crop"]
    m4_tidy[term == 'crop_typelegume',factor:= "Legume crop"]
    m4_tidy[term == 'crop_typemaize',factor:= "Maize crop"]
    m4_tidy[term == 'crop_typeroot crop',factor:= "Root crop"]
    m4_tidy[term == 'crop_residueResidue retained',factor:= "Residue retained"]
    
    # Define the order of the terms
    m4_tidy[,factor := factor(factor, levels = unique(factor))]
    
    # Add a column for significance stars
    m4_tidy[p.value <= 0.001, significance := '***']
    m4_tidy[p.value <= 0.01, significance := '**']
    m4_tidy[p.value <= 0.05, significance := '*']
    m4_tidy[p.value <= 0.10, significance := '.']
    m4_tidy[p.value > 0.10, significance := '']
    
    # remove two terms
    m4_tidy <- m4_tidy[!grepl('codeCC|codeNT',term)]
    
    # make plot
    p6 <- ggplot(m4_tidy, aes(x = factor, y = estimate)) + 
          geom_bar(stat = "identity", aes(fill = colorfactor), alpha = 0.7, fill = "lightblue") + 
          geom_errorbar(aes(y = estimate, ymin = estimate - std.error, ymax = estimate + std.error), width = 0.1) +
          geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + 0.001, estimate - 0.001)), 
                    vjust = ifelse(m4_tidy$estimate > 0, -0.5, 1.5), size = 4) + 
          theme_minimal() + 
          labs(x = 'Independent explanatory variables', y = 'Parameter estimate') + 
          theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
                axis.ticks.length = unit(-0.5, "cm"), 
                panel.grid = element_blank(), 
                legend.position = "none") + 
          ylim(-10, 10)
    
    ggsave(plot = p6, filename = 'products/plot_ph_finalmodel.png', width = 8,height = 6, units = "in", dpi = 300)
  }

# --- application of the model on the INTEGRATOR dataset ----
  
  # rerun the above code to regenerate model m3.full1 but now without standardization
  # save the model and use this model for application
  saveRDS(m4.full,file='products/mamodel_ph.rds')
  
  # read in the model
  m1 <- readRDS('products/mamodel_ph.rds')
  
  # load external datasets not stored at github given its size
  floc <- 'D:/DATA/17 nutribudget/'
  d4 <- fread(paste0(floc,'db_final_europe.csv'))
  
  # load in covariates
  d4.cov <- fread('products/240723_covariates_ncu.csv')
  
  # merge both files(note that not all NCUs are in covariates, to be checked later)
  d4 <- merge(d4,d4.cov,by.x='ncu',by.y='gncu2010_ext',all.x=TRUE)
  
  # predict the standardized pH change for each crop and site
  # note that this is done manually (from reading summary(m1) to avoid all re-naming of variables)
  
  # model coefficients
  m1.coeff <- as.data.table(broom::tidy(m1))
  
  # what is the mean and SD pH change
  SDp <- d2[, mean(sqrt((replication -1) * kpi_treat_sd^2 + (replication - 1)*kpi_contr_sd^2)/(2*replication - 2))]
  SMD <- d2[, mean((kpi_treat - kpi_control)/(sqrt((replication -1) * kpi_treat_sd^2 + (replication - 1)*kpi_contr_sd^2)/(2*replication - 2)))]
  
  # replacing missing inputs with median value
  d4[,parea.rtct := pmax(0,parea.rtct,na.rm=T)]
  d4[,parea.ntct := pmax(0,parea.ntct,na.rm=T)]
  d4[is.na(phw_mean_0_5), phw_mean_0_5 := median(d4$phw_mean_0_5,na.rm=T)]
  d4[is.na(tmp_mean), tmp_mean := median(d4$tmp_mean,na.rm=T)]
  
  # estimate the baseline response given the soil properties
  d4[, SMD := m1.coeff[grepl('mat',term),estimate] * tmp_mean]
  d4[, SMD := SMD + m1.coeff[grepl('^ph',term),estimate] * phw_mean_0_5 * 0.1]
  
  # add effects of crop type
  d4[grepl('grass',tolower(crop_name)), SMD := SMD + m1.coeff[grepl('crop_typegrass',term),estimate]]
  d4[grepl('fruit',tolower(crop_name)), SMD := SMD + m1.coeff[grepl('crop_typefruit',term),estimate]]
  d4[grepl('potat|sugar beet',tolower(crop_name)), SMD := SMD + m1.coeff[grepl('root crop',term),estimate]]
  d4[grepl('maize',tolower(crop_name)), SMD := SMD + m1.coeff[grepl('crop_typemaize',term),estimate]]
  d4[grepl('vegetabl',tolower(crop_name)), SMD := SMD + m1.coeff[grepl('crop_typelegume',term),estimate]]
  d4[grepl('puls|soya|rape|fodder|oil|rey|wheat|barly|oats|cerea|olive|fibre|wine|tobac|olive|industr',tolower(crop_name)), SMD := SMD + 0]
  
  # add effects of tillage and fertilizer site conditions
  d4[, SMD := SMD + pmax(0,parea.ntct - parea.rtct) *  m1.coeff[grepl('tillageno',term),estimate] / area_ncu_ha + 
                    parea.rtct * m1.coeff[grepl('tillagereduced',term),estimate]/ area_ncu_ha +
                    (area_ncu_ha - parea.ntct) * m1.coeff[grepl('tillage unknown',term),estimate] / area_ncu_ha]
  d4[, SMD := SMD + parea.nifnof * m1.coeff[grepl('fertilizer_typeNo',term),estimate]/ area_ncu_ha +
              (parea.mifhof + parea.hifhof) * m1.coeff[grepl('fertilizer_typeOrganic',term),estimate]/ area_ncu_ha+
              (parea.mifnof + parea.hifnof) * m1.coeff[grepl('fertilizer_typeMineral',term),estimate]/ area_ncu_ha]
  
  # effects of management
  d4[,SMD_BC := (SMD + m1.coeff[grepl('man_codeBC',term),estimate]) * SDp]
  d4[,SMD_CC := (SMD + m1.coeff[grepl('man_codeCC',term),estimate]) * SDp]
  d4[,SMD_CT := (SMD + m1.coeff[grepl('man_codeCT',term),estimate]) * SDp]
  d4[,SMD_LM := (SMD + m1.coeff[grepl('man_codeLM',term),estimate]) * SDp]
  d4[,SMD_NF := (SMD + m1.coeff[grepl('man_codeNF',term),estimate]) * SDp]
  d4[,SMD_OF := (SMD + m1.coeff[grepl('man_codeOF',term),estimate]) * SDp]
  d4[,SMD_RT := (SMD + m1.coeff[grepl('man_codeRT',term),estimate]) * SDp]
  d4[,SMD_ZT := (SMD + m1.coeff[grepl('man_codeZT',term),estimate]) * SDp]
  d4[, SMD_COMBI := (SMD + sum(sort(m1.coeff[grepl('man_codeRT|man_codeCC|man_codeOF|man_codeLM',term),estimate],decreasing = T) /c(1:4))) * SDp]
  
  # take weighted mean per ncu
  cols <- colnames(d4)[grepl('^SMD_|^ph$',colnames(d4))]
  
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
  pbreak <- c(-2,-1,0,0.5,1,1000)
  plabel <- c('< -1','-1 - 0','0 - 0.5','0.5 - 1.0','> 1')
  p1 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_BC, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of BC on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p1, filename = paste0(floc,'up_ph_meas_bc.jpg'),width = 12,height=12,units='cm')  
  
  p2 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_CC, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of cover crops on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p2, filename = paste0(floc,'up_ph_meas_cc.jpg'),width = 12,height=12,units='cm') 
  
  p3 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_CT, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of tillage on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p3, filename = paste0(floc,'up_ph_meas_ct.jpg'),width = 12,height=12,units='cm') 
  
  p4 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_LM, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of lime on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p4, filename = paste0(floc,'up_ph_meas_lm.jpg'),width = 12,height=12,units='cm') 
  
  p5 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_NF, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of inorganic fertilization on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p5, filename = paste0(floc,'up_ph_meas_nf.jpg'),width = 12,height=12,units='cm') 
  
  p6 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_OF, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of organic fertilization on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p6, filename = paste0(floc,'up_ph_meas_of.jpg'),width = 12,height=12,units='cm') 
  
  p7 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_RT, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of reduced tillage on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p7, filename = paste0(floc,'up_ph_meas_rt.jpg'),width = 12,height=12,units='cm') 
  
  p8 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_ZT, pbreak,labels = plabel))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Effect on\nsoil pH (-)')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Effect of zonation tillage on soil pH") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p8, filename = paste0(floc,'up_ph_meas_zt.jpg'),width = 12,height=12,units='cm') 
  
  p9 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= fifelse(ph >= 5.5,'realised','not realised'))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Target pH\n achieved')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Soil target pH achieved, baseline") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p9, filename = paste0(floc,'up_ph_target_bau.jpg'),width = 12,height=12,units='cm') 
  
  p10 <- ggplot() +
        geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
        geom_tile(data = r.ncu,aes(x=x,y=y,fill= fifelse(5.5 <= ph + SMD_COMBI,'realised','not realised'))) +
        scale_fill_viridis_d(direction=-1)+
        theme(legend.position.inside = c(0.1,0.8))+
        labs(fill = 'Target pH\nachieved')+
        xlab("Longitude") + ylab("Latitude") +
        ggtitle("Soil target pH achieved, combi of measures") +
        coord_sf(crs = 4326) + theme_bw()
  ggsave(plot = p10, filename = paste0(floc,'up_ph_target_combi.jpg'),width = 12,height=12,units='cm') 
  
  # give quantile
  data.table(quantile = paste0('Q',c(0.05,0.25,0.50,0.75,0.95)),
             BC = round(quantile(d5$SMD_BC,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             CC = round(quantile(d5$SMD_CC,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             CT = round(quantile(d5$SMD_CT,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             LM = round(quantile(d5$SMD_LM,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             NF = round(quantile(d5$SMD_NF,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             OF = round(quantile(d5$SMD_OF,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             RT = round(quantile(d5$SMD_RT,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             ZT = round(quantile(d5$SMD_ZT,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
             COMBI = round(quantile(d5$SMD_COMBI,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2))
  
  table(d5$ph <= 5.5)*100/nrow(d5)
  table(d5$ph + d5$SMD_COMBI<= 5.5)*100/nrow(d5)
  