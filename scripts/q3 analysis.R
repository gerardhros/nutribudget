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
  d2 <- d1[kpi_type =='ADG']

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
    theme_bw() + xlab('study id') + ylab('lnRR, change in average daily weight gain (ADG)') + 
    ggtitle('Observed changes in average daily weight gain (ADG) due to alternative protein source')

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
    ggtitle('effect of main factors on average daily weight gain (ADG)') +
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
                             mods = ~animal_type + supplement_category + supplemental_rate + stage-1,
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
  #d2[,sup_cat := fifelse(grepl('insect',supplement_category),'insect','other')] # the forest plot showed that "insect" had a positive effect while other supplement_category showed negative effects on DGR
  d2[,early_stage := fifelse(grepl('finishing|whole',stage),'other',stage)] 
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~ supplemental_rate *supplement_category + early_stage-1,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  # collect stats of the model
  estats(m3.full,m3.empty)
  
  # show anova whether full model is signifiantly different from empty model
  anova(m3.full,m3.empty,refit = T)
  
 # Extract and tidy the results
  m3_tidy <- tidy(m3.full)
  
  # Define new column values
  new_column_values <- c("supplemental_rate","supplemental_category","supplemental_category","supplemental_category","supplemental_category","supplemental_category","supplemental_category","stage","stage","interaction","interaction","interaction","interaction","interaction")
  
  # Add new column to data frame
  m3_tidy$factor <- new_column_values
  
  # Add a column for significance stars
  m3_tidy$significance <- ifelse(m3_tidy$p.value < 0.001, "***", 
                              ifelse(m3_tidy$p.value < 0.01, "**", 
                                     ifelse(m3_tidy$p.value < 0.05, "*",
                                            ifelse(m3_tidy$p.value < 0.1, "."))))
  print(m3_tidy$term)
  # Named vector with abbreviations
  abbreviations <- c("supplemental_rate"="rate","supplement_categoryduckweed"="duckweed","supplement_categorygrassjuice"="grassjuice","supplement_categoryinsect"="insect","supplement_categorylegumes"="legumes","supplement_categorymicroalgae"="microalgae","supplement_categoryrapeseed"="rapeseed","early_stageother"="late_stage","early_stagestarting"="early_stage", "supplemental_rate:supplement_categorygrassjuice"="rate+grassjuice","supplemental_rate:supplement_categoryinsect"="rate+insect","supplemental_rate:supplement_categorylegumes"="rate+legumes","supplemental_rate:supplement_categorymicroalgae"="rate+microalgae","supplemental_rate:supplement_categoryrapeseed"="rate+rapeseed")
  # Define the order of the terms
  m3_tidy$term <- factor(m3_tidy$term, levels = c("supplemental_rate","supplement_categoryduckweed","supplement_categorygrassjuice","supplement_categoryinsect","supplement_categorylegumes","supplement_categorymicroalgae","supplement_categoryrapeseed","early_stageother","early_stagestarting","supplemental_rate:supplement_categorygrassjuice","supplemental_rate:supplement_categoryinsect","supplemental_rate:supplement_categorylegumes","supplemental_rate:supplement_categorymicroalgae","supplemental_rate:supplement_categoryrapeseed"))
  
  # Create the bar plot with significance stars
  ggplot(m3_tidy, aes(x = term, y = estimate)) +
    geom_bar(stat = "identity",aes(fill = factor),alpha=0.7) +
    geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + 0.001, estimate - 0.001)), 
              vjust = ifelse(m3_tidy$estimate > 0, -0.5, 1.5), size = 5) +
    theme_minimal() +
    scale_x_discrete(labels = abbreviations)+ 
    labs(x = 'animal management practice & influencing factors',
         y = 'Parameter estimate',
         title = 'Meta-regression model estimation on the average daily weight gain (ADG)') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black frame line
          axis.ticks.length = unit(-0.5, "cm"),  # Minor tick marks outside
          panel.grid = element_blank(), # Remove grid lines
          legend.position = "none")+ # Remove legend
    ylim(-1.2, 1.5)  # Set the y-axis limits


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
  
  # rename column name
  # setnames(d2,old = c('duration..days.'),new = 'duration_days')
  
  # add id based on order yi
  d2[, id := frankv(yi,order = -1)]
  
  # make first plot of individual observations
  ggplot(data = d2,aes(x= id,y=yi)) + 
    geom_point() + 
    geom_line()+
    geom_errorbar(aes(x=id,ymin = yi-vi,ymax=yi+vi))+
    theme_bw() + xlab('study id') + ylab('lnRR, change in the digestibility of crude protein (DCP)') + 
    ggtitle('Observed changes in the digestibility of crude protein (DCP) due to alternative protein source')
  
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
    ggtitle('effect of main factors on digestability of crude protein (DCP)') +
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
  #d2[,sup_cat := fifelse(grepl('insect',supplement_category),'insect','other')]  # the forest plot showed that "insect" had a positive effect while other supplement_category showed negative effects on DGR
  #d2[,late_stage := fifelse(grepl('finishing',stage),'other',stage)]
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods = ~ supplemental_rate *supplement_category+stage-1,
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
  # Extract and tidy the results
  m3_tidy <- tidy(m3.full)
  
 # Define new column values
  new_column_values <- c("supplemental_rate","supplemental_category","supplemental_category","supplemental_category","supplemental_category","stage","stage","interaction","interaction")
  
  # Add new column to data frame
  m3_tidy$factor <- new_column_values
  
  # Add a column for significance stars
  m3_tidy$significance <- ifelse(m3_tidy$p.value < 0.001, "***", 
                                 ifelse(m3_tidy$p.value < 0.01, "**", 
                                        ifelse(m3_tidy$p.value < 0.05, "*", 
                                               ifelse(m3_tidy$p.value < 0.1, ".",""))))
  # Named vector with abbreviations
  abbreviations <- c("supplemental_rate"="rate","supplement_categoryinsect"="insect","supplement_categorylegumes"="legumes","supplement_categorymicroalgae"="microalgae","supplement_categoryrapeseed"="rapeseed","stagegrowing"="growing","stagewhole"="whole","supplemental_rate:supplement_categorylegumes"="rate+legumes","supplemental_rate:supplement_categoryrapeseed"="rate+rapeseed")
  # Define the order of the terms
  m3_tidy$term <- factor(m3_tidy$term, levels = c("supplemental_rate","supplement_categoryinsect","supplement_categorylegumes","supplement_categorymicroalgae","supplement_categoryrapeseed","supplemental_rate:supplement_categorylegumes","stagegrowing","stagewhole","supplemental_rate:supplement_categoryrapeseed"))
  
  # Create the bar plot with significance stars
  ggplot(m3_tidy, aes(x = term, y = estimate)) +
    geom_bar(stat = "identity",aes(fill = factor),alpha=0.7) +
    geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + 0.001, estimate - 0.001)), 
              vjust = ifelse(m3_tidy$estimate > 0, -0.5, 1.5), size = 5) +
    theme_minimal() +
    scale_x_discrete(labels = abbreviations)+ 
    labs(x = 'animal management practice factors',
         y = 'Parameter estimate',
         title = 'Meta-regression model estimation on the digestability of crude protein(DCP)') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black frame line
          axis.ticks.length = unit(-0.5, "cm"),  # Minor tick marks outside
          panel.grid = element_blank(), # Remove grid lines
          legend.position = "none")+ # Remove legend
    ylim(-0.15, 0.10)  # Set the y-axis limits
  
   
# --- Analysis for KPI 3 -----
  ## this KPI is excluded in D1.3, as the dataset is not enough to develop a meta-regression model



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
    theme_bw() + xlab('study id') + ylab('lnRR, change in feed conversion ratio (FCR)') + 
    ggtitle('Observed changes in feed conversion ratio (FCR) due to alternative protein source')
  
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
    ggtitle('effect of main factors on feed conversion ratio (FCR)') +
    theme(legend.position = 'bottom')+
    scale_y_discrete(limits = c("supplemental_rate","pig", "starting", "growing","finishing","whole", "rapeseed","legumes", "insect", "grassjuice","duckweed","microalgae"))
  
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
                             mods = ~animal_type + supplement_category + supplemental_rate + stage-1,
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
  #d2[,early_stage := fifelse(grepl('starter',stage),'other',stage)] 
  # Find the unique values for each variable
  sapply(lapply(d2, unique), length)
  m3.full <- metafor::rma.mv(yi,vi, 
                             mods =  ~ supplemental_rate* supplement_category + animal_type + stage * supplement_category-1,
                             data=d2,random= list(~ 1|study_ID), method="REML",sparse = TRUE)
  
  # analyse summary stats
  summary(m3.full)
  
  # collect stats of the model
  estats(m3.full,m3.empty)

  # show anova whether full model is significantly different from empty model
  anova(m3.full,m3.empty,refit = T)
 
  # Extract and tidy the results
  m3_tidy <- tidy(m3.full)
  
  # Define new column values
  new_column_values <- c("supplemental_rate","supplemental_category","supplemental_category","supplemental_category","supplemental_category","supplemental_category","supplemental_category","animal_type","animal_stage","animal_stage","animal_stage","interaction","interaction","interaction","interaction","interaction","interaction","interaction")
  
  # Add new column to data frame
  m3_tidy$factor <- new_column_values
  
  # Add a column for significance stars
  m3_tidy$significance <- ifelse(m3_tidy$p.value < 0.001, "***", 
                                 ifelse(m3_tidy$p.value < 0.01, "**", 
                                        ifelse(m3_tidy$p.value < 0.05, "*",
                                               ifelse(m3_tidy$p.value < 0.1, ".",""))))
  
  # Named vector with abbreviations
  abbreviations <- c("supplemental_rate"="rate","supplement_categoryduckweed"="duckweed","supplement_categorygrassjuice"="grassjuice","supplement_categoryinsect"="insect","supplement_categorylegumes"="legumes","supplement_categorymicroalgae"="microalage","supplement_categoryrapeseed"="rapeseed","supplemental_rate:supplement_categorygrassjuice"="rate+grassjuice","supplemental_rate:supplement_categoryinsect"="rate+insect","supplemental_rate:supplement_categorylegumes"="rate+legumes","supplemental_rate:supplement_categorymicroalgae"="rate+microalgae","supplemental_rate:supplement_categoryrapeseed"="rate+rapeseed","supplement_categoryinsect:stagegrowing"="insect at growing stage","supplement_categorylegumes:stagegrowing"="legumes at growing stage","animal_typepig"="pig","stagegrowing"="growing","stagestarting"="starting","stagewhole"="whole")
  # Define the order of the terms
  m3_tidy$term <- factor(m3_tidy$term, levels = c("supplemental_rate","supplement_categoryduckweed","supplement_categorygrassjuice","supplement_categoryinsect","supplement_categorylegumes","supplement_categorymicroalgae","supplement_categoryrapeseed","supplemental_rate:supplement_categorygrassjuice","supplemental_rate:supplement_categoryinsect","supplemental_rate:supplement_categorylegumes","supplemental_rate:supplement_categorymicroalgae","supplemental_rate:supplement_categoryrapeseed","supplement_categoryinsect:stagegrowing","supplement_categorylegumes:stagegrowing","animal_typepig","stagegrowing","stagestarting","stagewhole"))
  
  # Create the bar plot with significance stars
  ggplot(m3_tidy, aes(x = term, y = estimate)) +
    geom_bar(stat = "identity",aes(fill = factor),alpha=0.7) +
    geom_text(aes(label = significance, y = ifelse(estimate > 0, estimate + 0.005, estimate - 0.005)), 
              vjust = ifelse(m3_tidy$estimate > 0, -0.5, 1.5), size = 3) +
    theme_minimal() +
    scale_x_discrete(labels = abbreviations)+ 
    labs(x = 'animal management practice & influencing factors',
         y = 'Parameter estimate',
         title = 'Meta-regression model estimation on the feed conversion ratio (FCR)') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black frame line
          axis.ticks.length = unit(-0.5, "cm"),  # Minor tick marks outside
          panel.grid = element_blank(), # Remove grid lines
          legend.position = "none") + # Remove legend
    ylim(-2, 2) # Set the y-axis limits with a bit more space

