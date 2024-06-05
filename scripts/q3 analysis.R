# data analysis RQ3
# Gerard H. ros & Hongzheng
# 13 May 2024

# require packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)

# read in the database
d1 <- as.data.table(read_xlsx('data/RQ3_database.xlsx',sheet='RQ3_database'))

# what are the unique KPIs included
table(d1$kpi_type)

# estimate missing SD and n

  # add coefficient of variation per KPI type
  d1[,cv_treat := kpi_treat_sd / kpi_treat_mean,by= 'kpi_type']
  d1[,cv_control := kpi_control_sd / kpi_control_mean,by= 'kpi_type']
  
  # add mean CV for full dataset
  d1[,cv_treat_mean := mean(cv_treat,na.rm=TRUE),by= 'kpi_type']
  d1[,cv_control_mean := mean(cv_control,na.rm=TRUE),by= 'kpi_type']

  # replace missing SD by 1.25 * mean CV * mean
  d1[is.na(kpi_treat_sd), kpi_treat_sd := cv_treat_mean * 1.25 * kpi_treat_mean]
  d1[is.na(kpi_control_sd), kpi_control_sd := cv_control_mean * 1.25 * kpi_control_mean]
  
  # remove columns not needed any more
  d1[,c('cv_treat','cv_control','cv_treat_mean','cv_control_mean') := NULL]
  
# --- Analysis for KPI 1 -----
  
  # subset the file
  d2 <- d1[kpi_type =='DGR']

  # remove observations without a control or treatment value
  d2 <- d2[!(is.na(kpi_treat_mean) | is.na(kpi_control_mean))]
  
  # check the data.table on missing values
  summary(d2)
  
  # rename the column with a difficult name to select
  # setnames(d2,old = 'duration (days)', new = 'duration_days')
  
  # i see that duration has 3 missing values, so i replaced these values by the median value
  d2[is.na(duration_days), duration_days := median(d2$duration_days,na.rm=TRUE)]
  
  # add response ratio as effect size; you can select another by adapting "measure" argument
  d2 <- escalc(measure = "ROM", data = d2, 
               m1i = kpi_treat_mean , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control_mean, sd2i = kpi_control_sd, n2i = replication)
  
  # convert back to data.table
  d2 <- as.data.table(d2)

  # remove some colums to make visualisation in console easier
  d2[,c('location','man_treatment','man_control','control_composition','year','mat','map') := NULL]
  
  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]
  
  # make first plot of individual observations
  ggplot(data = d2,aes(x= id,y=yi)) + 
    geom_point() + 
    geom_line()+
    geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
    theme_bw() + xlab('study id') + ylab('lnRR, change in DGR') + 
    ggtitle('Observed changes in DGR due to novel protein source')

  # estimate mean response across all studies
  m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)
    
  # see the summary of the model
  summary(m1)
  
  # analyse whether the mean response differs per factor
  fcols <- c('animal_type','supplement_category','supplemental_rate','stage')
   
  # make an empty list to store output
  out.dgr <- list()
  
  # do the analysis in a for loop
  for(i in fcols){
    
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
                            random= list(~ 1|study_ID), method="REML",sparse = TRUE)
      
      # collect model coefficients
      m2.sum <- summary(m2)
      m2.out <- as.data.table(coefficients(m2.sum))
      m2.out[,factor := i]
      m2.out[,var := c('intercept',i)]
      
    }
    
    # add kpi type
    m2.out[,kpi :='DGR']
    
    # add to list
    out.dgr[[i]] <- copy(m2.out)
  }
  
  # convert list to one data.table
  out.dgr <- rbindlist(out.dgr)
  
  # make barplot (figure formatting: need to be done)
  out.dgr[, pfactor := as.factor(factor)]
  ggplot(data=out.dgr[var!='intercept']) +
    geom_bar(aes(x=estimate,y= reorder(var,pfactor),group = pfactor,fill=pfactor),stat="identity") + 
    geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
    ggtitle('effect of main factors on DGR') +
    theme(legend.position = 'bottom')+
    scale_y_discrete(limits = c("supplemental_rate","chichen","pig", "starter", "growing","finishing","whole", "rapeseed","legumes", "insect", "duckweed","microalgae"))
    
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
                             mods = ~animal_type + supplement_category + supplemental_rate + stage,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)

  # refine the model (if desired)
  # steps to do: add variables plus interactions (* for all interactions plus main factors, : for interactions only), check pvalue
  # if within a variable one subgroup is significant and others not, then adjust the groupings
  # so, you have the options:
  # only a model with additive factors => use the plus sign in the mods argument
  # if you have a model with main factors AND interactions => use the "*" sign in the mods argument
  # if you have a model with ONLY interactions => use then the ":" sign
  d2[,sup_cat := fifelse(grepl('insect',supplement_category),'insect','other')] # the forest plot showed that "insect" had a positive effect while other supplement_category showed negative effects on DGR
  d2[,early_stage := fifelse(grepl('finishing|whole',stage),'other',stage)] 
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~ supplemental_rate + animal_type +early_stage + early_stage * sup_cat,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
 # collect stats of the model
  estats(m3.full,m3.empty)
  
  # show anova whether full model is signifiantly different from empty model
  anova(m3.full,m3.empty,refit = T)
  

# --- Analysis for KPI 2 -----
  
  # subset the file
  d2 <- d1[kpi_type =='DCP']
  
  # remove observations without a control or treatment value
  d2 <- d2[!(is.na(kpi_treat_mean) | is.na(kpi_control_mean))]
  
  # add response ratio
  d2 <- escalc(measure = "ROM", data = d2, 
               m1i = kpi_treat_mean , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control_mean, sd2i = kpi_control_sd, n2i = replication)
  
  # convert back to data.table
  d2 <- as.data.table(d2)
  
  # remove some colums to make visualisation in console easier
  d2[,c('location','man_treatment','man_control','control_composition','year','mat','map') := NULL]
  
 
  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]
  
  # make first plot of individual observations
  ggplot(data = d2,aes(x= id,y=yi)) + 
    geom_point() + 
    geom_line()+
    geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
    theme_bw() + xlab('study id') + ylab('lnRR, change in DCP') + 
    ggtitle('Observed changes in DCP due to novel protein source')
  
  # estimate mean response across all studies
  m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)
  
  # see the summary of the model
  summary(m1)
  
  # analyse whether the mean response differs per factor
  # note that there should be variation, factors with a single factor gives an error
  fcols <- c('supplement_category','supplemental_rate','stage')
  
  # make an empty list to store output
  out.dcp <- list()
  
  # do the analysis in a for loop
  for(i in fcols){
    
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
                            random= list(~ 1|study_ID), method="REML",sparse = TRUE)
      
      # collect model coefficients
      m2.sum <- summary(m2)
      m2.out <- as.data.table(coefficients(m2.sum))
      m2.out[,factor := i]
      m2.out[,var := c('intercept',i)]
      
    }
    
    # add kpi type
    m2.out[,kpi :='DCP']
    
    # add to list
    out.dcp[[i]] <- copy(m2.out)
  }
  
  # convert list to one data.table
  out.dcp <- rbindlist(out.dcp)
  
  # make barplot (figure formatting: need to be done)
  out.dcp[, pfactor := as.factor(factor)]
  ggplot(data=out.dcp[var!='intercept']) +
    geom_bar(aes(x=estimate,y= reorder(var,pfactor),group = pfactor,fill=pfactor),stat="identity") + 
    geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
    ggtitle('effect of main factors on DCP') +
    theme(legend.position = 'bottom')
  
  # do one meta-regression model with multiple factors
  
  # make first an empty model
  m3.empty <- metafor::rma.mv(yi,vi, data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # make a full model with all main factors together    
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~supplement_category + supplemental_rate + stage -1,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  
  # refine the model (if desired)
  # steps to do: add variables plus interactions (* for all interactions plus main factors, : for interactions only), check pvalue
  # if within a variable one subgroup is significant and others not, then adjust the groupings
  d2[,sup_cat := fifelse(grepl('insect',supplement_category),'insect','other')]
  d2[,late_stage := fifelse(grepl('finishing',stage),'other',stage)]
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~ supplemental_rate +supplement_category+ late_stage * sup_cat,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  # analyse summary stats
  summary(m3.full)
    
  # collect stats of the model
  estats(m3.full,m3.empty)
  
  # show anova whether full model is signifiantly different from empty model
  anova(m3.full,m3.empty,refit = T)
  
  # assume you like to now the impact of protein source on animal type = 'pig', stage ='finishing',syp_cat='rapeseed',
  # and supplemental_rate = 50
  # predict: 50 * 0.0501 + 1 * -0.137 + 0 * -0.11 + 0* -0.260 + 0.1112 + 0 * 0.0727 - 0.0221 * 50 * 0
  
  
# --- Analysis for KPI 3 -----
  
  # there are quite some cases without a control. 
  ## the code failed, as dataset is not enough to develop a meta-regression model
 
  # subset the file
  d2 <- d1[kpi_type =='NUE']
  
  # remove observations without a control or treatment value
  d2 <- d2[!(is.na(kpi_treat_mean) | is.na(kpi_control_mean))]
  # check the data.table on missing values
  summary(d2)
  
  # rename the column with a difficult name to select
  setnames(d2,old = 'duration_days', new = 'duration_days')
  
  # i see that duration has 3 missing values, so i replaced these values by the median value
  d2[is.na(duration_days), duration_days := median(d2$duration_days,na.rm=TRUE)]
  
  # add response ratio as effect size; you can select another by adapting "measure" argument
  d2 <- escalc(measure = "ROM", data = d2, 
               m1i = kpi_treat_mean , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control_mean, sd2i = kpi_control_sd, n2i = replication)
  
  # convert back to data.table
  d2 <- as.data.table(d2)
  
  # remove some columns to make visualisation in console easier
  d2[,c('location','man_treatment','man_control','control_composition','year','mat','map') := NULL]
  
  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]
  
  # make first plot of individual observations
  ggplot(data = d2,aes(x= id,y=yi)) + 
    geom_point() + 
    geom_line()+
    geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
    theme_bw() + xlab('study id') + ylab('lnRR, change in NUE') + 
    ggtitle('Observed changes in NUE due to novel protein source')
  
  # estimate mean response across all studies
  m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)
  
  # see the summary of the model
  summary(m1)
  
  # analyse whether the mean response differs per factor
  fcols <- c('animal_type','supplement_category','supplemental_rate','stage')
  
  # make an empty list to store output
  out.nue <- list()
  
  # do the analysis in a for loop
  for(i in fcols){
    
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
                            random= list(~ 1|study_ID), method="REML",sparse = TRUE)
      
      # collect model coefficients
      m2.sum <- summary(m2)
      m2.out <- as.data.table(coefficients(m2.sum))
      m2.out[,factor := i]
      m2.out[,var := c('intercept',i)]
      
    }
    
    # add kpi type
    m2.out[,kpi :='NUE']
    
    # add to list
    out.nue[[i]] <- copy(m2.out)
  }
  
  sapply(lapply(d2, unique), length)
  # convert list to one data.table
  out.nue <- rbindlist(out.nue)
  
  # make barplot (figure formatting: need to be done)
  out.nue[, pfactor := as.factor(factor)]
  ggplot(data=out.nue[var!='intercept']) +
    geom_bar(aes(x=estimate,y= reorder(var,pfactor),group = pfactor,fill=pfactor),stat="identity") + 
    geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
    ggtitle('effect of main factors on NUE') +
    theme(legend.position = 'bottom')+
    scale_y_discrete(limits = c("supplemental_rate","chichen","pig", "starter", "growing","finishing","whole", "rapeseed","legumes", "insect", "duckweed","microalgae"))
  
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
                             mods = ~animal_type + supplement_category + supplemental_rate + stage,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  
  # refine the model (if desired)
  # steps to do: add variables plus interactions (* for all interactions plus main factors, : for interactions only), check pvalue
  # if within a variable one subgroup is significant and others not, then adjust the groupings
  # so, you have the options:
  # only a model with additive factors => use the plus sign in the mods argument
  # if you have a model with main factors AND interactions => use the "*" sign in the mods argument
  # if you have a model with ONLY interactions => use then the ":" sign
  d2[,sup_cat := fifelse(grepl('insect',supplement_category),'insect','other')] 
  d2[,early_stage := fifelse(grepl('finishing|whole',stage),'other',stage)] 
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~ supplemental_rate + animal_type +early_stage + early_stage * sup_cat,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  
  # collect stats of the model
  estats(m3.full,m3.empty)
  
  # show anova whether full model is signifiantly different from empty model
  anova(m3.full,m3.empty,refit = T)

  # --- Analysis for KPI 4 -----
  
  # there are quite some cases without a control.
  
  # subset the file
  d2 <- d1[kpi_type =='FCR']
  
  # remove observations without a control or treatment value
  d2 <- d2[!(is.na(kpi_treat_mean) | is.na(kpi_control_mean))]
  # check the data.table on missing values
  summary(d2)
  
  # rename the column with a difficult name to select
  setnames(d2,old = 'duration_days', new = 'duration_days')
  
  # i see that duration has 3 missing values, so i replaced these values by the median value
  d2[is.na(duration_days), duration_days := median(d2$duration_days,na.rm=TRUE)]
  
  # add response ratio as effect size; you can select another by adapting "measure" argument
  d2 <- escalc(measure = "ROM", data = d2, 
               m1i = kpi_treat_mean , sd1i = kpi_treat_sd , n1i = replication,
               m2i = kpi_control_mean, sd2i = kpi_control_sd, n2i = replication)
  
  # convert back to data.table
  d2 <- as.data.table(d2)
  
  # remove some columns to make visualisation in console easier
  d2[,c('location','man_treatment','man_control','control_composition','year','mat','map') := NULL]
  
  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]
  
  # make first plot of individual observations
  ggplot(data = d2,aes(x= id,y=yi)) + 
    geom_point() + 
    geom_line()+
    geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
    theme_bw() + xlab('study id') + ylab('lnRR, change in FCR') + 
    ggtitle('Observed changes in FCR due to novel protein source')
  
  # estimate mean response across all studies
  m1 <- metafor::rma.mv(yi,vi, 
                        data=d2,
                        random= list(~ 1|study_ID), 
                        method="REML",
                        sparse = TRUE)
  
  # see the summary of the model
  summary(m1)
  
  # analyse whether the mean response differs per factor
  fcols <- c('animal_type','supplement_category','supplemental_rate','stage')
  
  # make an empty list to store output
  out.fcr <- list()
  
  # do the analysis in a for loop
  for(i in fcols){
    
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
                            random= list(~ 1|study_ID), method="REML",sparse = TRUE)
      
      # collect model coefficients
      m2.sum <- summary(m2)
      m2.out <- as.data.table(coefficients(m2.sum))
      m2.out[,factor := i]
      m2.out[,var := c('intercept',i)]
      
    }
    
    # add kpi type
    m2.out[,kpi :='FCR']
    
    # add to list
    out.fcr[[i]] <- copy(m2.out)
  }
  
  # convert list to one data.table
  out.fcr <- rbindlist(out.fcr)
  
  # make barplot (figure formatting: need to be done)
  out.fcr[, pfactor := as.factor(factor)]
  ggplot(data=out.fcr[var!='intercept']) +
    geom_bar(aes(x=estimate,y= reorder(var,pfactor),group = pfactor,fill=pfactor),stat="identity") + 
    geom_errorbar(aes(y=var,xmin = estimate - se,xmax = estimate +se),width=0.4) + theme_bw() +
    ggtitle('effect of main factors on FCR') +
    theme(legend.position = 'bottom')+
    scale_y_discrete(limits = c("supplemental_rate","chichen","pig", "starter", "growing","finishing","whole", "rapeseed","legumes", "insect", "grassjuice","duckweed","microalgae"))
  
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
                             mods = ~animal_type + supplement_category + supplemental_rate + stage,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  
  # refine the model (if desired)
  # steps to do: add variables plus interactions (* for all interactions plus main factors, : for interactions only), check pvalue
  # if within a variable one subgroup is significant and others not, then adjust the groupings
  # so, you have the options:
  # only a model with additive factors => use the plus sign in the mods argument
  # if you have a model with main factors AND interactions => use the "*" sign in the mods argument
  # if you have a model with ONLY interactions => use then the ":" sign
  #d2[,sup_cat := fifelse(grepl('rapeseed|grassjuice|insect',supplement_category),'other',supplement_category)] # the forest plot showed that the rapeseed, grassjuice and insect have a negative effect on the FCR while other supplement_catagory showed positive effect. therefore they are sub-grouped.
  d2[,early_stage := fifelse(grepl('starter',stage),'other',stage)] 
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~ supplemental_rate + animal_type +early_stage + early_stage * supplement_category,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  
  # collect stats of the model
  estats(m3.full,m3.empty)
  
  # show anova whether full model is signifiantly different from empty model
  anova(m3.full,m3.empty,refit = T)
