# data analysis RQ1
# Maria Anna Antonovardaki & Gerard Ros
# 22 May 2024

#required packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)

# read in the database (NOTE that I updated some SE and SD values)
d1a <- as.data.table(read_xlsx('data/yield_database.xlsx',sheet=1))
d1b <- as.data.table(read_xlsx('data/yield_database.xlsx',sheet=2))
d1c <- as.data.table(read_xlsx('data/yield_database.xlsx',sheet=3))

# combine the three sheets
d1 <- rbind(d1a,d1b,d1c)

# update column names
setnames(d1,
         old = c('clay (%)','Decimal degree lat','Decimal degree lon'),
         new = c('clay','lat_deg','lon_deg'))

# save the unique lon-lat for covariate extraction
# the covariates are collected in script coariate extraction

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
  d1 <- merge(d1,d1.cov,by.x=c('lat_deg','lon_deg'),by.y=c('lat','lon'),all.x=TRUE)

# estimate uncertainty when SE and SD are missing
  
  # remove SD values exceeding the mean value (input error?, to be checked by MA)
  d1[kpi_treat_sd > kpi_treat, kpi_treat_sd := NA_real_]
  d1[kpi_contr_sd > kpi_control, kpi_contr_sd := NA_real_]
  d1[kpi_treat_se > kpi_treat, kpi_treat_se := NA_real_]
  d1[kpi_contr_se > kpi_control, kpi_contr_se := NA_real_]
  
  # first fill-up dataset for missing input data (i think forgotten by MA)
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
  d1[,c('lat_deg','lon_deg','lat','lon','reference') := NULL]
  
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
  d1[, clay2 := eclay(clay),by=dataset_ID] # note that this changes range into mean value
  
# --- Analysis for KPI crop yield

# remove observations without a control or treatment value
d2 <- d1[!(is.na(kpi_treat) | is.na(kpi_control))]

# check the data.table on missing values
summary(d2)


# duration have many missing values, so i replaced these values by the median value
d2[is.na(d2$duration), duration := median(d2$duration,na.rm =TRUE)]

#add the standart deviation when standart error is present
d2[, kpi_contr_sd1 := kpi_contr_se * sqrt(replication)]
d2[, kpi_treat_sd1 := kpi_treat_se * sqrt(replication)]

#fill in missing values in standart deviations by adding a median value
d2[is.na (d2$kpi_contr_sd1),kpi_contr_sd1 := median(d2$kpi_contr_sd1, 
                                                    na.rm = TRUE) ]
d2[is.na (d2$kpi_treat_sd1),kpi_treat_sd := median(d2$kpi_treat_sd1,
                                                   na.rm = TRUE) ]

# add response ratio as effect size
d2 <- escalc(measure = "ROM", data = d2, 
             m1i = kpi_treat , sd1i = kpi_treat_sd1 , n1i = replication,
             m2i = kpi_control, sd2i= kpi_contr_sd1, n2i = replication)

# convert back to data.table
d2 <- as.data.table(d2)

# remove some colums to make visualisation in console easier
d2[,c('man_treatment','man_control','year','mat','map') := NULL]

# add id based on order yi
d2[, id := frankv(yi,order = -1)]

# make first plot of individual observations
ggplot(data = d2,aes(x= id,y=yi)) + 
  geom_point() + 
  geom_line()+
  geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
  theme_bw() + xlab('study id') + ylab('lnRR, change in Crop yield') + 
  ggtitle('Observed changes in Crop yield due to Crop practices')


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
           'fertilizer_type')

# make an empty list to store output
out.dgr <- list()

# Why we do that? do the analysis in a for loop
for(i in fcols)
  
  # make a temporary variable 'evar'
  d2[,evar := get(i)]

# is the variable categorial or not
vartype = is.character(d2[,get(i)])

# do the meta-analysis for categorial variable
if(vartype==TRUE){
  
  # do the meta-analysis
  m2 <- metafor::rma.mv(yi,vi, mods = ~factor(evar)-1,data=d2,
                        random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # collect model coefficients
  m2.sum <- summary(m2)
  m2.out <- as.data.table(coefficients(m2.sum))
  m2.out[,factor := i]
  m2.out[,var := gsub('factor\\(evar\\)','',rownames(m2.sum$b))]
  
} else{
  
  # do the meta-analysis
  m2 <- metafor::rma.mv(yi,vi, mods = ~evar,data=d2,
                        random= list(~ 1|study_ID), method="REML",sparse = TRUE)}
# add kpi type
m2.out[,d2$kpi]

# add to list
out.dgr[[i]] <- copy(m2.out)

# convert list to one data.table
out.dgr <- rbindlist(out.dgr)

# make barplot (figure formatting: need to be done)
out.dgr[, pfactor := as.factor(factor)]
ggplot(data=out.dgr[var!='intercept']) +
  geom_bar(aes(x=estimate,y= reorder(var,pfactor),group = pfactor,fill=pfactor),stat="identity") + 
  geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
  ggtitle('effect of main factors on Crop yield') +
  theme(legend.position = 'bottom')

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
                           mods = ~crop_type + crop_residue + cover_crop + crop_rotation + tillage + fertilizer_type,
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