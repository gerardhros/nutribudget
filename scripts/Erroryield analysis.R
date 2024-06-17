# data analysis RQ1
# Maria Anna Antonovardaki & Gerard Ros
# 22 May 2024

#required packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)

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
  
  # STARTING OF ERROR CONVERGENCE = 1 given the huge changes and deline with a factor 3 is unlikely
  d2 <- d2[abs(yi) <= 2]

  d2 <- d2[kpi_treat_sd>1 & kpi_contr_sd > 1]
  
  # estimate mean response across all studies
m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)

  # see the summary of the model
  summary(m1)
  #END OF ERROR CONVERGENCE
  
  #PART OF MARIANNA TO AVOID CONVERGENCE = 1
  # Transform Data and Exclude Outliers Based on Effect Size
  d2 <- d2[abs(yi) <= 2]
  
  # Check the distribution of sampling variances and exclude extreme values
  summary(d2$vi)
  hist(d2$vi, breaks=50, main="Distribution of Sampling Variances", xlab="Sampling Variance (vi)")
  
  # Define a threshold to exclude extreme variances (e.g., remove the top 1% of variances)
  threshold <- quantile(d2$vi, 0.99)
  d2 <- d2[d2$vi < threshold, ]
  
  # Apply log transformation to effect sizes and variances
  d2[, log_yi := log(abs(yi) + 1)]
  d2[, log_vi := log(vi + 1)]
  
  # Function to fit the model with different optimizers
  fit_model <- function(data) {
    optimizers <- c("nlminb", "optim", "bobyqa")
    optmethods <- c(NULL, "BFGS", NULL)
    
    for (i in 1:length(optimizers)) {
      optimizer <- optimizers[i]
      optmethod <- optmethods[i]
      
      m1 <- tryCatch({
        metafor::rma.mv(log_yi, log_vi, 
                        data=data,
                        random=list(~ 1|study_ID), 
                        method="REML",
                        sparse=TRUE,
                        verbose=TRUE,
                        control=list(optimizer=optimizer, optmethod=optmethod, 
                                     rel.tol=1e-6, stepadj=0.5, maxit=10000))
      }, warning=function(w) {
        message("Model encountered a warning with optimizer ", optimizer, ": ", w)
        return(NULL)
      }, error=function(e) {
        message("Model failed with optimizer ", optimizer, ": ", e)
        return(NULL)
      })
      
      if (!is.null(m1)) {
        return(m1)
      }
    }
    
    return(NULL)
  }
  
  # Fit the model and handle potential warnings/errors
  m1 <- fit_model(d2)
  
  if (!is.null(m1)) {
    print(summary(m1))
  } else {
    print("Model did not converge or failed with all optimizers.")
  }
#END FOR CONVERGENCE

# analyse whether the mean response differs per factor
fcols <- c('crop_type','crop_residue','cover_crop','crop_rotation', 'tillage', 
           'fertilizer_type','ph','man_code','fertilizer_type',
           'nutri_dose_N','nutri_dose_P','nutri_dose_K','nutri_dose_Mg',
           'bdod_mean_0_5','cec_mean_0_5','clay_mean_0_5','ntot_mean_0_5',
           'phw_mean_0_5','soc_mean_0_5',
           'tmp_mean','pet_mean','tmp_sd','pet_sd',
           'ctype')

# PROTOTYPE CODE GETTING AN ERROR that optimizer nlminb did not achieved convergence = 1
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
out.dgr[, pfactor := as.factor(factor)]
ggplot(data=out.dgr[var!='intercept']) +
  geom_bar(aes(x=estimate,y= var),stat="identity") + 
  geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
  ggtitle('effect of main factors on Crop yield') +
  theme(legend.position = 'bottom')

#END OF PROTOTYPE CODE GETTING AN ERROR that optimizer nlminb did not achieved convergence = 1

#Part FIRTY TRY/PART OF MARIANNA/trying to avoid the error about convergence = 1 
#create an empty list
out.dgr <- list()

#start the loop
for(i in fcols) {
  
  # make a temporary variable 'evar'
  d2[, evar := get(i)]
  # Remove rows with NA in evar
  d2 <- d2[!is.na(evar)]
  
  # Check for extreme variances and exclude them based on a quantile threshold
  threshold <- quantile(d2$vi, 0.99)
  d2 <- d2[d2$vi < threshold, ]
  
  # is the variable categorial or not
  vartype = is.character(d2[, get(i)])
  
  # do the meta-analysis for categorial variable
  if(vartype == TRUE) {
    
    # do the meta-analysis for categorial variable
    m2 <- tryCatch({
      metafor::rma.mv(yi, vi, mods = ~factor(evar), data=d2,
                      random= list(~ 1|study_ID), method="REML", sparse = TRUE, 
                      verbose=TRUE, control=list(rel.tol=1e-6, stepadj=0.5, maxit=10000))
    }, warning = function(w) {
      message("Model encountered a warning: ", w)
      return(NULL)
    }, error = function(e) {
      message("Model failed with error: ", e)
      return(NULL)
    })
    
    if (!is.null(m2)) {
      m2.sum <- summary(m2)
      m2.out <- as.data.table(m2.sum$coefficients, keep.rownames = TRUE)
      m2.out[, factor := i]
      m2.out[, var := gsub('factor\\(evar\\)', '', rn)]
      m2.out[, rn := NULL]
      out.dgr[[i]] <- m2.out
      print(paste0('Model built for parameter: ', i))
    } else {
      print(paste0('Model failed for parameter: ', i))
    }
    
  } else {
    
    # do the meta-analysis for a numerical variable
    m2 <- tryCatch({
      metafor::rma.mv(yi, vi, mods = ~evar, data=d2,
                      random= list(~ 1|study_ID), method="REML", sparse = TRUE, 
                      verbose=TRUE, control=list(rel.tol=1e-6, stepadj=0.5, maxit=10000))
    }, warning = function(w) {
      message("Model encountered a warning: ", w)
      return(NULL)
    }, error = function(e) {
      message("Model failed with error: ", e)
      return(NULL)
    })
    
    if (!is.null(m2)) {
      m2.sum <- summary(m2)
      m2.out <- as.data.table(m2.sum$coefficients, keep.rownames = TRUE)
      m2.out[, factor := i]
      m2.out[, var := c('intercept', i)]
      m2.out[, rn := NULL]
      out.dgr[[i]] <- m2.out
      print(paste0('Model built for parameter: ', i))
    } else {
      print(paste0('Model failed for parameter: ', i))
    }
    
  }
}

# convert list to one data.table
out.dgr <- rbindlist(out.dgr, use.names = TRUE, fill = TRUE)

# remove cropname (too much options, use only grouped one)
out.dgr <- out.dgr[factor != 'crop_type']

# make barplot (figure formatting: need to be done)
out.dgr[, pfactor := as.factor(factor)]
ggplot(data=out.dgr[var!='intercept']) +
  geom_bar(aes(x=estimate, y= var), stat="identity") + 
  geom_errorbar(aes(y=var, xmin = estimate - se, xmax = estimate + se), width=0.4) + 
  theme_bw() +
  ggtitle('Effect of main factors on Crop yield') +
  theme(legend.position = 'bottom')

#END OF 1st try/PART OF MARIANNA

#Start 2nd try/part of Marianna/changing the number of thershold

# Transform Data and Exclude Outliers Based on Effect Size
d2 <- d2[abs(yi) <= 2]

# Define a threshold to exclude extreme variances based on the histogram
# Here, I'm using 0.02 as an example threshold
variance_threshold <- 0.02
d2 <- d2[vi < variance_threshold, ]

# Apply log transformation to effect sizes and variances
d2[, log_yi := log(abs(yi) + 1)]
d2[, log_vi := log(vi + 1)]

# Function to fit the model with different optimizers
fit_model <- function(data) {
  optimizers <- c("nlminb", "optim", "bobyqa")
  optmethods <- c(NULL, "BFGS", NULL)
  
  for (i in 1:length(optimizers)) {
    optimizer <- optimizers[i]
    optmethod <- optmethods[i]
    
    m1 <- tryCatch({
      metafor::rma.mv(log_yi, log_vi, 
                      data=data,
                      random=list(~ 1|study_ID), 
                      method="REML",
                      sparse=TRUE,
                      verbose=TRUE,
                      control=list(optimizer=optimizer, optmethod=optmethod, 
                                   rel.tol=1e-6, stepadj=0.5, maxit=10000))
    }, warning=function(w) {
      message("Model encountered a warning with optimizer ", optimizer, ": ", w)
      return(NULL)
    }, error=function(e) {
      message("Model failed with optimizer ", optimizer, ": ", e)
      return(NULL)
    })
    
    if (!is.null(m1)) {
      return(m1)
    }
  }
  
  return(NULL)
}

# Fit the model and handle potential warnings/errors
m1 <- fit_model(d2)

if (!is.null(m1)) {
  print(summary(m1))
} else {
  print("Model did not converge or failed with all optimizers.")
}

#END OF SECOND TRY/PART OF MARIANNA

#START THIRD try/part of Marianna/changing the number of thershold after calculating with GetData GraphDigitizer

# Transform Data and Exclude Outliers Based on Effect Size
d2 <- d2[abs(yi) <= 2]

# Define a threshold to exclude extreme variances based on the histogram
# I'm using 0.014 to see if the threshold is the reason on why it doesn't work and somehow harmonise the variances
variance_threshold <- 0.014
d2 <- d2[vi < variance_threshold, ]

# Apply log transformation to stabilise the effect sizes and variances (It has also been checked without this part but still the code does not work)
d2[, log_yi := log(abs(yi) + 1)]
d2[, log_vi := log(vi + 1)]

# Function to fit the model with different optimizers
fit_model <- function(data) {
  optimizers <- c("nlminb", "optim", "bobyqa")
  optmethods <- c(NULL, "BFGS", NULL)
  
  for (i in 1:length(optimizers)) {
    optimizer <- optimizers[i]
    optmethod <- optmethods[i]
    
    m1 <- tryCatch({
      metafor::rma.mv(log_yi, log_vi, 
                      data=data,
                      random=list(~ 1|study_ID), 
                      method="REML",
                      sparse=TRUE,
                      verbose=TRUE,
                      control=list(optimizer=optimizer, optmethod=optmethod, 
                                   rel.tol=1e-6, stepadj=0.5, maxit=10000))
    }, warning=function(w) {
      message("Model encountered a warning with optimizer ", optimizer, ": ", w)
      return(NULL)
    }, error=function(e) {
      message("Model failed with optimizer ", optimizer, ": ", e)
      return(NULL)
    })
    
    if (!is.null(m1)) {
      return(m1)
    }
  }
  
  return(NULL)
}

# Fit the model and handle potential warnings/errors
m1 <- fit_model(d2)

if (!is.null(m1)) {
  print(summary(m1))
} else {
  print("Model did not converge or failed with all optimizers.")
}

#END OF THIRD TRY/PART OF MARIANNA


#START FOURTH TRY/PART OF MARIANNA/Trying to simplify the model

d2 <- d2[abs(yi) <= 2]

# Define a threshold to exclude extreme variances based on the histogram
variance_threshold <- 0.02
d2 <- d2[vi < variance_threshold, ]

# Apply log transformation to effect sizes and variances
d2[, log_yi := log(abs(yi) + 1)]
d2[, log_vi := log(vi + 1)]

# Check for data issues: missing values, data types, sample sizes
print(summary(d2))
print(str(d2))

# Ensure there are no missing values
d2 <- na.omit(d2)

# Function to fit the model with different optimizers
fit_model <- function(data) {
  optimizers <- c("nlminb", "optim", "bobyqa")
  optmethods <- c(NULL, "BFGS", NULL)
  
  for (i in 1:length(optimizers)) {
    optimizer <- optimizers[i]
    optmethod <- optmethods[i]
    
    m1 <- tryCatch({
      metafor::rma.mv(log_yi, log_vi, 
                      data=data,
                      random=list(~ 1|study_ID), 
                      method="REML",
                      sparse=TRUE,
                      verbose=TRUE,
                      control=list(optimizer=optimizer, optmethod=optmethod, 
                                   rel.tol=1e-6, stepadj=0.5, maxit=10000))
    }, warning=function(w) {
      message("Model encountered a warning with optimizer ", optimizer, ": ", w)
      return(NULL)
    }, error=function(e) {
      message("Model failed with optimizer ", optimizer, ": ", e)
      return(NULL)
    })
    
    if (!is.null(m1)) {
      return(m1)
    }
  }
  
  return(NULL)
}

# Fit a simplified model without random effects
fit_simplified_model <- function(data) {
  tryCatch({
    m2 <- metafor::rma(yi=log_yi, vi=log_vi, 
                       data=data, 
                       method="REML",
                       control=list(optimizer="nlminb", rel.tol=1e-6, stepadj=0.5, maxit=10000))
    return(m2)
  }, warning=function(w) {
    message("Simplified model encountered a warning: ", w)
    return(NULL)
  }, error=function(e) {
    message("Simplified model failed with error: ", e)
    return(NULL)
  })
}

# Try fitting the complex model
m1 <- fit_model(d2)

if (!is.null(m1)) {
  print(summary(m1))
} else {
  message("Complex model did not converge. Trying a simplified model...")
  # If the complex model fails, try a simplified model
  m2 <- fit_simplified_model(d2)
  if (!is.null(m2)) {
    print(summary(m2))
  } else {
    print("Simplified model did not converge or failed.")
  }
}
#END FOURTH TRY/PART OF MARIANNA

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

# make a full model with all main factors together    
m3.full <- metafor::rma.mv(yi,vi, 
                           mods = ~ctype + crop_residue + cover_crop + crop_rotation + 
                                   tillage + fertilizer_type + nutri_dose_N + ntot_mean_0_5 +
                                   ph + tmp_mean + pet_mean,
                           data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)

# analyse summary stats
summary(m3.full)


# I DO NOT UNDERSTAND THIS PART refine the model (if desired)
# steps to do: add variables plus interactions (* for all interactions plus main factors, : for interactions only), check pvalue
# if within a variable one subgroup is significant and others not, then adjust the groupings
# so, you have the options:
# only a model with additive factors => use the plus sign in the mods argument
# if you have a model with main factors AND interactions => use the "*" sign in the mods argument
# if you have a model with ONLY interactions => use then the ":" sign
d2[,sup_cat := fifelse(grepl('micro',crop_type),'x','other')]
m3.full <- metafor::rma.mv(yi,vi, 
                           mods = ~ supplemental_rate + stage : sup_cat,
                           data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)

# It doesn't work probably because of above collect stats of the model
estats(m3.full,m3.empty)

# It doesn't work probabl ybecause of above show anova whether full model is signifiantly different from empty model
anova(m3.full,m3.empty,refit = T)