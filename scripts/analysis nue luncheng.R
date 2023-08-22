# Do a meta-analysis on experimental data
# example script, prepared by Gerard Ros, gerard.ros@wur.nl, 23-aug-23.

# Load libraries 
library(data.table); library(metafor); library(metagear)

# -- data preparation ----

  # read data from Excel and convert into a data.table
  d1 <- readxl::read_xlsx('data/luncheng_2023_nue.xlsx',sheet = "Tables")
  d1 <- as.data.table(d1)

  # when data on variance (SD) is missing for the NUE of the control and treated plot, estimate this from the CV of the other studies
  d1[, nuet_cv := mean(nuet_sd/nuet_mean,na.rm=T) * 1.25]
  d1[, nuec_cv := mean(nuec_sd/nuec_mean,na.rm=T) * 1.25]
  d1[is.na(nuet_sd), nuet_sd := nuet_mean * nuet_cv]
  d1[is.na(nuec_sd), nuec_sd := nuec_mean * nuec_cv]
  
  # clean up column names (remove dashes, and make everyting lower case)
  setnames(d1,gsub('\\/','_',gsub(' |\\(|\\)','',colnames(d1))))
  setnames(d1,tolower(colnames(d1)))

  # update the missing values for n_dose (replace by median N dose when missing)
  d1[is.na(n_dose), n_dose := median(d1$n_dose,na.rm=TRUE)]
  
  # # scale the variables to unit variance
  d1[,clay_scaled := scale(clay)]
  d1[,soc_scaled := scale(soc)]
  d1[,ph_scaled := scale(ph)]
  d1[,mat_scaled := scale(mat)]
  d1[,map_scaled := scale(map)]
  d1[,n_dose_scaled := scale(n_dose)]
  
  # calculate the effect size for NUE using the function escalc and the log-transformed response ratio as effect size
  # you might also use other effect sizes, see ?escalc
  es21 <- escalc(measure = "ROM", data = d1, 
                 m1i = nuet_mean, sd1i = nuet_sd, n1i = replication,
                 m2i = nuec_mean, sd2i = nuec_sd, n2i = replication )
  
  # convert to data.table
  d02 <- as.data.table(es21)

  # what are the treatments to be assessed
  d02.treat <- data.table(treatment =  c('ALL',unique(d02$management)))
  
  # add description of the labels
  d02.treat[treatment=='ALL',desc := 'All']
  d02.treat[treatment=='EE',desc := 'Enhanced Efficiency']
  d02.treat[treatment=='CF',desc := 'Combined fertilizer']
  d02.treat[treatment=='RES',desc := 'Residue retention']
  d02.treat[treatment=='RFP',desc := 'Fertilizer placement']
  d02.treat[treatment=='RFR',desc := 'Fertilizer rate']
  d02.treat[treatment=='ROT',desc := 'Crop rotation']
  d02.treat[treatment=='RFT',desc := 'Fertilizer timing']
  d02.treat[treatment=='OF',desc := 'Organic fertilizer']
  d02.treat[treatment=='RT',desc := 'Reduced tillage']
  d02.treat[treatment=='NT',desc := 'No tillage']
  d02.treat[treatment=='CC',desc := 'Crop cover']


# --- main factor analysis ----
  
  # do a main factor analysis of the factors controlling the response variable
  # for more details on the function used, see ?rma.mv
  
  # add a list to store the coefficients
  out2 = out3 = list()

  # make a for loop to do a main analysis per treatment
  for(i in d02.treat$treatment){
    
    if(i=='ALL'){
      
      # run without selection to estimate overall mean, with a random error structured by studyID
      r_nue <- metafor::rma.mv(yi,vi, 
                               data=d02,
                               random= list(~ 1|studyid), 
                               method="REML",
                               sparse = TRUE)
      
    } else {
      
      # run for selected treatment with a random error structured by studyID
      r_nue <- metafor::rma.mv(yi,vi, 
                               data=d02[management==i,],
                               random= list(~ 1|studyid), 
                               method="REML",
                               sparse = TRUE)
      
    }
    
    # save relevant output in a list: the model coefficient (corrected to relative change in %), the standard error, the pvalue and the label
    out2[[i]] <- data.table(mean = as.numeric((exp(r_nue$b)-1)*100),
                            se = as.numeric((exp(r_nue$se)-1)*100),
                            pval = round(as.numeric(r_nue$pval),4),
                            label =  paste0(d02.treat[treatment==i,desc],' (n=',r_nue$k,')')
                            )
    # print to console to see the progress
    print(paste0('treatment ',i,' has been simulated'))
  }

# convert all output from the list out2 to a data.table
out2 <- rbindlist(out2)

# --- plot a forest plot -----

  # plot for NUE (using default plot function from metafor)
  p1 <- forest(x = out2$mean, 
               sei = out2$se, 
               slab=out2$label, psize=0.9, cex=1, sortvar=out2$label, xlab="Change in NUE (%)", header="Treatment", col="#CC0000", lwd=2)

  # make a more nicely designed ggplot
  out2[,summary := fifelse(grepl('All',label),1,2)]
  p2 <- ggplot(data = out2, aes(x = label, y =mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,size=factor(summary),colour = factor(summary)))+
        geom_pointrange() + coord_flip() + geom_hline(aes(yintercept=0),lty=2,size=1)+
        scale_size_manual(values=c(1,0.5))+
        xlab("Study")+ylab("Percentage change in NUE")+
        scale_colour_manual(values=c("black","grey"))+
        scale_y_continuous(breaks=c(-20,0,20,40,60),limits = c(-40,60)) + theme_bw()+
        theme(legend.position="none")+
        ggtitle('Effect of agronomic measures on NUE')
  ggsave(plot = p2, filename = 'products/luncheng_nue_forestplot.png',width = 15, height = 10, units = 'cm')


# --- publication bias tests ----

  # begg’s test
  ranktest(out2$mean, sei=out2$se)

  # egger’s test
  regtest(out2$mean, sei = out2$se)

# --- meta-regression for main factors ----

  # do a first main factor analysis for log response ratio for NUE

  # what are the factors to be evaluated
  var.site <- c('mat_scaled','map_scaled','clay_scaled','soc_scaled','ph_scaled')
  var.crop <- c('g_crop_type','n_dose_scaled')
  var.trea <- c('fertilizer_type', 'crop_residue', 'tillage', 'cover_crop_and_crop_rotation', 'fertilizer_strategy')

  # the columns to be assessed
  var.sel <- c(var.trea,var.crop,var.site)

  # run without a main factor selection to estimate overall mean
  r_nue_0 <- rma.mv(yi,vi, data = d02,random= list(~ 1|studyid), method="REML",sparse = TRUE)

  # objects to store the effects per factor as wel summary stats of the meta-analytical models
  out1.est = out1.sum = list()

  # evaluate the impact of treatment (column tillage) on NUE given site properties
  for(i in var.sel){
    
    # check whether the column is a numeric or categorical variable
    vartype = is.character(d02[,get(i)])
    
    # run with the main factor treatment
    if(vartype == TRUE){
      
      # run a meta-regression model for main categorial variable (note that intercept is removed)
      r_nue_1 <- rma.mv(yi,vi, 
                        mods = ~factor(varsel)-1, 
                        data = d02[,.(yi,vi,studyid,varsel = get(i))],
                        random = list(~ 1|studyid), method="REML",sparse = TRUE)
      
    } else {
      
      # run a meta-regression model for main numerical variable
      r_nue_1 <- rma.mv(yi,vi, 
                        mods = ~varsel, 
                        data = d02[,.(yi,vi,studyid,varsel = get(i))],
                        random = list(~ 1|studyid), method="REML",sparse = TRUE)
    }
    
    # save output in a list: the estimated impact of the explanatory variable
    out1.est[[i]] <- data.table(var = i,
                                varname = gsub('factor\\(varsel\\)','',rownames(r_nue_1$b)),
                                mean = round(as.numeric(r_nue_1$b),3),
                                se = round(as.numeric(r_nue_1$se),3),
                                ci.lb = round(as.numeric(r_nue_1$ci.lb),3),
                                ci.ub = round(as.numeric(r_nue_1$ci.ub),3),
                                pval = round(as.numeric(r_nue_1$pval),3))
    
    # save output in a list: the summary stats collected
    out1.sum[[i]] <- data.table(var = i,
                                AIC = r_nue_1$fit.stats[4,2],
                                ll = r_nue_1$fit.stats[1,2],
                                ll_impr = round(100 * (1-r_nue_1$fit.stats[1,2]/r_nue_0$fit.stats[1,2]),2),
                                r2_impr = round(100*max(0,(sum(r_nue_0$sigma2)-sum(r_nue_1$sigma2))/sum(r_nue_0$sigma2)),2),
                                pval = round(anova(r_nue_1,r_nue_0)$pval,3)
    )
    
  }
  
  # merge output into a data.table for both summary stats as well the coefficients
  out1.sum <- rbindlist(out1.sum)
  out1.est <- rbindlist(out1.est)


# --- meta-regression for main factors with interactions

  # make a function to extract relevant model statistics
  estats <- function(model_new,model_base){
    out <- data.table(AIC = model_new$fit.stats[4,2],
                      ll = model_new$fit.stats[1,2],
                      ll_impr = round(100 * (1-model_new$fit.stats[1,2]/model_base$fit.stats[1,2]),2),
                      r2_impr = round(100*max(0,(sum(model_base$sigma2)-sum(model_new$sigma2))/sum(model_base$sigma2)),2),
                      pval = round(anova(r_nue_1,r_nue_0)$pval,3))
    return(out)
  }

  # update a few groups in the data because numbers per category are too small
  d02[tillage=='reduced', tillage := 'no-till']
  d02[,fertilizer_type := factor(fertilizer_type,levels = c('mineral','organic', 'combined','enhanced'))] 
  d02[,fertilizer_strategy := factor(fertilizer_strategy,levels = c("conventional", "placement","rate","timing"))]
  d02[,g_crop_type := factor(g_crop_type, levels = c('maize','wheat','rice'))]

  # make a new local copy of the data.table
  d4 <- copy(d02)

  # add one_dot encoding manually for a few variables
  d4[,r4pl := fifelse(fertilizer_strategy=='placement','yes','no')]
  d4[,r4ti := fifelse(fertilizer_strategy=='timing','yes','no')]
  d4[,r4do := fifelse(fertilizer_strategy=='rate','yes','no')]
  d4[,ctm := fifelse(g_crop_type=='maize','yes','no')]
  d4[,ctw := fifelse(g_crop_type=='wheat','yes','no')]
  d4[,ctr := fifelse(g_crop_type=='rice','yes','no')]
  d4[,ndose2 := scale(n_dose^2)]


  # run without a main factor selection to estimate overall mean
  r_nue_0 <- rma.mv(yi,vi, data = d02,random= list(~ 1|studyid), method="REML",sparse = TRUE)

  # this is a model, manually developed by adding variables/ testing impact on model performnance
  # the improvement for each improvement can be tested by the function estats. if improved, accept model parameter and add new one
  # the addition is done based on expected impact and actual improvement
  r_nue_4 <- rma.mv(yi,vi,
                  mods = ~fertilizer_type + r4pl + r4ti + r4do + crop_residue + tillage +
                    cover_crop_and_crop_rotation + n_dose_scaled + clay_scaled + ph_scaled + map_scaled + mat_scaled + soc_scaled+
                    soc_scaled : n_dose_scaled + ctm:r4pl + ctm + ctw + ctr + ctm:mat_scaled  + ndose2 -1,
                  data = d4,
                  random = list(~ 1|studyid), method="REML",sparse = TRUE)


# show stats and improvements
out = estats(model_new = r_nue_4,model_base = r_nue_0)
print(paste0('model improved the log likelyhood with ',round(out$ll_impr,1),'%'))
summary(r_nue_4)

k <- r_nue_4$k
wi <- 1/r_nue_4$vi
vt <- (k-1) / (sum(wi) - sum(wi^2)/sum(wi))
PR2 <- r_nue_0$sigma2 / (sum(r_nue_4$sigma2) + vt)

