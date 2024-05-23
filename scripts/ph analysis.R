# data analysis RQ on pH
# Salim Belyazid & Gerard Ros
# 23 May 2024

#required packages
require(readxl);require(data.table);require(metafor);require(ggplot2)

# clear environment
rm(list=ls())

# read in the database (NOTE that I updated some SE and SD values
d1 <- as.data.table(read_xlsx('data/ph_database.xlsx',sheet='database'))

# update column names to simply internal references to the column names (so no spaces or brackets)
setnames(d1,old = c('CEC (meq/100g)'),new = c('cec'))
setnames(d1,tolower(colnames(d1)))

# remove text from latitude column, and make the column numeric (text is aumatically converted to NA)
d1[,lat := as.numeric(lat)]

# save the unique lon-lat for covariate extraction
# the covariates are collected in script coariate extraction
# below the preparation is in an ifelse statement, and the prepared file is read in
# subset only relevant columns, and write csv file, only done once, therefore in FALSE ifelse statement
if(FALSE){
  # subset only relevant columns
  d1.lonlat <- d1[,.(lat,lon)]
  # select only the cases where latitude is given
  d1.lonlat <- unique(d1.lonlat[!is.na(lat)])
  fwrite(d1.lonlat,'data/ph_lonlat_unique.csv')
  
} else {
  
  # after collection site properties, read these covariates in
  d1.cov <- fread('products/240523_covariates_ph.csv')
}

# combine the extracted coariates with the original dataset
d1 <- merge(d1,d1.cov,by = c('lat','lon'),all.x=TRUE)

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
cols <- c('duration','mat','map','soc','cec','ph','clay','bs')

# apply the function to expected numeric variables  
d1[, c(cols) := lapply(.SD,function(x) hfunconvert(x)),.SDcols = cols,by=id]

# regroup crop type to minimize options (to be adapted)
d1[,crop_type := tolower(crop_type)]
d1[grepl('barley|wheat|oat|cereal|rye|spelt',crop_type),ctype := 'cereal']
d1[is.na(ctype) & grepl('vegetable|carrot|onion|melon|tomato|pepper|cabbage',crop_type),ctype := 'vegetable']
d1[is.na(ctype) & grepl('bean|pea|rape|oil',crop_type),ctype := 'beans_oilcrop']
d1[is.na(ctype) & grepl('potato|sugar beet',crop_type),ctype := 'arable']
d1[is.na(ctype) & grepl('maize',crop_type),ctype := 'maize']
d1[is.na(ctype) & grepl('rice',crop_type),ctype := 'rice']
d1[is.na(ctype) & grepl('grass',crop_type),ctype := 'grass']
d1[is.na(ctype),ctype := 'unknown']

# regroup fertilizer type
d1[, fertilizer_type := tolower(fertilizer_type)]
d1[is.na(fertilizer_type), fertilizer_type := 'unknown']
d1[grepl('^mineral$|inorganic|different|npk|np|^n|amonium nitrate',fertilizer_type), fertilizer_type := 'inorganic']
d1[grepl('manure|slurry|fym|vegetal',fertilizer_type), fertilizer_type := 'organic']

# regroup nutricode (not in use)

# add unknown
d1[is.na(cover_crop), cover_crop := 'unknown']
d1[is.na(crop_rotation), crop_rotation := 'unknown']
d1[is.na(crop_residue), crop_residue := 'unknown']
d1[is.na(tillage)|grepl('conv|Plow|Plough|^Till',tillage), tillage := 'CT']
d1[grepl('no till|Mineral ',tillage), tillage := 'NT']

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
ggplot(data = d2,aes(x= id,y=yi)) + 
  geom_point() + 
  geom_line()+
  geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
  theme_bw() + xlab('study id') + ylab('lnRR, change in Crop yield') + 
  ggtitle('Observed changes in Crop yield due to Crop practices')

# given the huge changes and deline with a factor 3 is unlikely
d2 <- d2[abs(yi)<=2]

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