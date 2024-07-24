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
  es21 <- escalc(measure = "SMD", data = d1, 
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

if(FALSE){
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
}

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
  d2 <- copy(d02)

  # add one_dot encoding manually for a few variables
  d2[,r4pl := fifelse(fertilizer_strategy=='placement','yes','no')]
  d2[,r4ti := fifelse(fertilizer_strategy=='timing','yes','no')]
  d2[,r4do := fifelse(fertilizer_strategy=='rate','yes','no')]
  d2[,ctm := fifelse(g_crop_type=='maize','yes','no')]
  d2[,ctw := fifelse(g_crop_type=='wheat','yes','no')]
  d2[,ctr := fifelse(g_crop_type=='rice','yes','no')]
  d2[,ndose2 := scale(n_dose^2)]


  # run without a main factor selection to estimate overall mean
  r_nue_0 <- rma.mv(yi,vi, data = d02,random= list(~ 1|studyid), method="REML",sparse = TRUE)

  # this is a model, manually developed by adding variables/ testing impact on model performnance
  # the improvement for each improvement can be tested by the function estats. if improved, accept model parameter and add new one
  # the addition is done based on expected impact and actual improvement
  r_nue_4 <- rma.mv(yi,vi,
                  mods = ~fertilizer_type + r4pl + r4ti + r4do + crop_residue + tillage +
                    cover_crop_and_crop_rotation + n_dose_scaled + clay_scaled + ph_scaled + map_scaled + mat_scaled + soc_scaled+
                    soc_scaled : n_dose_scaled + ctm:r4pl + ctm + ctw + ctr + ctm:mat_scaled  + ndose2 -1,
                  data = d2,
                  random = list(~ 1|studyid), method="REML",sparse = TRUE)


# show stats and improvements
out = estats(model_new = r_nue_4,model_base = r_nue_0)
print(paste0('model improved the log likelyhood with ',round(out$ll_impr,1),'%'))
summary(r_nue_4)

# --- application of the model on the INTEGRATOR dataset ----

  # save the model and use this model for application
  saveRDS(r_nue_4,file='products/mamodel_nue.rds')

  # read in the model
  m1 <- readRDS('products/mamodel_nue.rds')

  # load external datasets not stored at github given its size
  floc <- 'D:/DATA/17 nutribudget/'
  d4 <- fread(paste0(floc,'db_final_europe.csv'))
  
  # load in covariates
  d4.cov <- fread('products/240723_covariates_ncu.csv')
  
  # merge both files(note that not all NCUs are in covariates, to be checked later)
  d4 <- merge(d4,d4.cov,by.x='ncu',by.y='gncu2010_ext',all.x=TRUE)

  # predict the NUE change for each crop and site
  # note that this is done manually (from reading summary(m1) to avoid all re-naming of variables)
  
  # model coefficients
  m1.coeff <- as.data.table(broom::tidy(m1))

  # what is the mean and SD NUE change
  SDp <- d2[, mean(sqrt(((replication -1) * nuet_sd^2 + (replication - 1)*nuec_sd^2)/(2*replication - 2)))]
  SMD <- d2[, mean((nuet_mean - nuec_mean)/SDp)]
  
# rescale the variables to unit variance
  d4[, clay_scaled := (clay - mean(d1$clay,na.rm=T))/sd(d1$clay,na.rm=T)]
  d4[, soc_scaled := (soc_ref - mean(d1$soc,na.rm=T))/sd(d1$soc,na.rm=T)]
  d4[, ph_scaled := (ph - mean(d1$ph,na.rm=T))/sd(d1$ph,na.rm=T)]
  d4[, n_dose_scaled := (n_fert+n_man - mean(d1$n_dose,na.rm=T))/sd(d1$n_dose,na.rm=T)]
  d4[, mat_scaled := (tmp_mean - mean(d1$mat,na.rm=T))/sd(d1$mat,na.rm=T)]
  d4[, map_scaled := scale(pre_mean)]
  d4[, n_dose2 := ((n_fert+n_man)^2 - mean(d1$n_dose^2,na.rm=T))/sd(d1$n_dose^2,na.rm=T)]
  
  # replacing missing inputs with median value
  d4[,parea.rtct := pmax(0,parea.rtct,na.rm=T)]
  d4[,parea.ntct := pmax(0,parea.ntct,na.rm=T)]
  d4[is.na(mat_scaled), mat_scaled := median(d4$mat_scaled,na.rm=T)]
  d4[is.na(map_scaled), map_scaled := median(d4$map_scaled,na.rm=T)]

  # add one_dot encoding manually for a few variables
  d4[, nue_cur := n_covercrop * yield_ref * 0.001 /(n_fert + n_man + n_fix + n_dep)]
  d4[, nin := n_fert + n_man + n_fix + n_dep]
  
  # site properties used in the model
  d4[,ctmyes := fifelse(grepl('maiz',tolower(crop_name)),1,0)]
  d4[,ctwyes := fifelse(grepl('barl|cerea|wheat|oats|rey',tolower(crop_name)),1,0)]
  d4[,ctryes := fifelse(grepl('rice',tolower(crop_name)),1,0)]
  d4[,r4plyes := fifelse(nue_cur < 0.4 | n_fert > 40,1,0)]
  d4[,r4tiyes := fifelse(nue_cur < 0.4 | n_fert > 40,1,0)]
  d4[,r4doyes := fifelse(nue_cur < 0.4 | n_fert > 40,1,0)]
  
  # estimate the baseline response given the soil and climatic properties
  d4[, SMD := m1.coeff[term == 'clay_scaled',estimate] * clay_scaled]
  d4[, SMD := SMD + m1.coeff[term == 'ph_scaled',estimate] * ph_scaled]
  d4[, SMD := SMD + m1.coeff[term == 'map_scaled',estimate] * map_scaled]
  d4[, SMD := SMD + m1.coeff[term == 'mat_scaled',estimate] * mat_scaled]
  d4[, SMD := SMD + m1.coeff[term == 'soc_scaled',estimate] * soc_scaled]
  d4[, SMD := SMD + m1.coeff[term == 'n_dose_scaled',estimate] * n_dose_scaled]
  d4[ctmyes == 1, SMD := SMD + m1.coeff[term == 'ctmyes',estimate]]
  d4[ctwyes == 1, SMD := SMD + m1.coeff[term == 'ctwyes',estimate]]
  d4[, SMD := SMD + m1.coeff[term == 'ndose2',estimate] * n_dose2]
  d4[, SMD := SMD + m1.coeff[term == 'n_dose_scaled:soc_scaled',estimate] * n_dose_scaled*soc_scaled]
  d4[r4plyes == 1 & ctmyes ==1, SMD := SMD + m1.coeff[term == 'r4plyes:ctmyes',estimate]]
  d4[ctwyes == 1, SMD := SMD + m1.coeff[term == 'mat_scaled:ctmyes',estimate]*mat_scaled]
  d4[r4plyes == 1, SMD := SMD + m1.coeff[term == 'r4plyes',estimate]]
  d4[r4tiyes == 1, SMD := SMD + m1.coeff[term == 'r4tiyes',estimate]]
  d4[r4doyes == 1, SMD := SMD + m1.coeff[term == 'r4doyes',estimate]]
  
  # effect of lowering N dose (38% reduction)
  d4[, SMD2 := 0]
  d4[, SMD2 := SMD2 + m1.coeff[term == 'n_dose_scaled',estimate] * (n_dose_scaled -1)]
  d4[, SMD2 := SMD2 + m1.coeff[term == 'n_dose_scaled:soc_scaled',estimate] * (n_dose_scaled-1)*soc_scaled]
  d4[, SMD2 := SMD2 + m1.coeff[term == 'ndose2',estimate] * (n_dose2 -1)]
  
  # effect of management
  d4[, SMD_FD := (SMD + (parea.mifhof + parea.hifhof+parea.mifnof + parea.hifnof) * SMD2/ area_ncu_ha)*SDp]
  d4[, SMD_NT := (SMD + parea.ntct  *  m1.coeff[grepl('tillageno',term),estimate] / area_ncu_ha)*SDp ]
  d4[, SMD_CR := (SMD + (pmax(0,parea.cr,na.rm=T) * m1.coeff[grepl('crop_residueyes',term),estimate])/area_ncu_ha)*SDp]
  d4[, SMD_CC := (SMD + (pmax(0,parea.cc,na.rm=T) * m1.coeff[grepl('cover_crop_and_crop_rotationyes',term),estimate])/area_ncu_ha)*SDp]
  d4[, SMD_FT := (SMD + (parea.mifhof + parea.hifhof+parea.mifnof + parea.hifnof) * m1.coeff[grepl('fertilizer_typemineral',term),estimate]/ area_ncu_ha)*SDp]
  d4[, SMD_EF := (SMD + (parea.mifhof + parea.hifhof+parea.mifnof + parea.hifnof) * m1.coeff[grepl('fertilizer_typeenhanced',term),estimate]/ area_ncu_ha)*SDp]
  d4[, SMD_OF := (SMD + (parea.nifnof + parea.mifnof + parea.hifnof) * m1.coeff[grepl('fertilizer_typeorganic',term),estimate]/ area_ncu_ha)*SDp]
  d4[, SMD_CF := (SMD + (parea.nifnof + parea.mifhof + parea.hifhof+parea.mifnof + parea.hifnof) * m1.coeff[grepl('fertilizer_typecombined',term),estimate]/ area_ncu_ha)*SDp]
  d4[, SMD_COMBI := SMD_FD + SMD_NT + SMD_CR + SMD_CC + SMD_FT + SMD_EF + SMD_OF + SMD_CF - 7*SMD*SDp]
  
  # take weighted mean per ncu
  cols <- colnames(d4)[grepl('^SMD_|nue_|^n_sp|^nin',colnames(d4))]
  
  d5 <- d4[,lapply(.SD,function(x) weighted.mean(x,w=area_ncu)),.SDcols = cols,by='ncu']

  # estimate the change in Nsp
  #cols <- colnames(d4)[grepl('^SMD_',colnames(d4))]
  #d5[,c(cols) := lapply(.SD,function(x) n_sp_ref - x * 0.01 *nin),.SDcols = cols]
  
# --- plot impact of measures on crop yield ----

  # load in NUTS shapefile
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

# plot impact of measures on soil N surplus
  #pbreak <- c(-100,20,20,60,80,100,1000)
  #plabel <- c('< 20','20-40','40-60','60-80','80-100','>100')
  pbreak <- c(-1000,2,4,6,8,10,15,25,1000)
  plabel <- c('< 2','2-4','4-6','6-8','8-10','10-15','15-25','>25')
  p1 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_NT, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n(abs-%)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of no-till on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p1, filename = paste0(floc,'up_nue_meas_nt.jpg'),width = 12,height=12,units='cm')  

p2 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_CR, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n(abs-%)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of crop rotation on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p2, filename = paste0(floc,'up_nue_meas_cr.jpg'),width = 12,height=12,units='cm') 

p3 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_CC, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n (abs-%) (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of cover crops on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p3, filename = paste0(floc,'up_nue_meas_cc.jpg'),width = 12,height=12,units='cm') 

p4 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_EF, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n (abs-%) (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of efficiency fertilizers on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p4, filename = paste0(floc,'up_nue_meas_ef.jpg'),width = 12,height=12,units='cm') 

p5 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_CF, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n (abs-%) (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of combined fertilizers on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p5, filename = paste0(floc,'up_nue_meas_cf.jpg'),width = 12,height=12,units='cm') 

p6 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_FT, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n (abs-%) (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of timing fertilizers on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p6, filename = paste0(floc,'up_nue_meas_ft.jpg'),width = 12,height=12,units='cm') 

p6 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_FD, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n (abs-%) (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of fertilizer dose on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p6, filename = paste0(floc,'up_nue_meas_fd.jpg'),width = 12,height=12,units='cm') 

p7 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= cut(SMD_COMBI, pbreak,labels = plabel))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Effect on\nNUE\n (abs-%) (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Effect of combined measures on NUE") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p7, filename = paste0(floc,'up_nue_meas_combi.jpg'),width = 12,height=12,units='cm') 

r.ncu[,ta_bau := fifelse(pmin(n_sp_sw_crit, n_sp_gw_crit) <= n_sp_ref,1,0)]
r.ncu[,ta_best := fifelse(pmin(n_sp_sw_crit, n_sp_gw_crit) <= n_sp_ref - (1 - SMD_COMBI) * 0.01 * nin,1,0)]
p8 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= fifelse(pmin(n_sp_sw_crit, n_sp_gw_crit) <= n_sp_ref ,'realised','not realised'))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Target N-surplus\n achieved')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Target N surplus achieved, baseline") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p8, filename = paste0(floc,'up_nsp_target_bau.jpg'),width = 12,height=12,units='cm') 

p9 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= fifelse(pmin(n_sp_sw_crit, n_sp_gw_crit) <= n_sp_ref - (1 - SMD_COMBI) * 0.01 *nin,'realised','not realised'))) +
      scale_fill_viridis_d(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'Target N-surplus\nachieved')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Target N surplus achieved, combi of measures") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p9, filename = paste0(floc,'up_nsp_target_combi.jpg'),width = 12,height=12,units='cm') 

# give quantile
data.table(quantile = paste0('Q',c(0.05,0.25,0.50,0.75,0.95)),
           FD = round(quantile(d5$SMD_FD,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           CR = round(quantile(d5$SMD_CR,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           CC = round(quantile(d5$SMD_CC,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           FT = round(quantile(d5$SMD_FT,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           EF = round(quantile(d5$SMD_EF,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           OF = round(quantile(d5$SMD_OF,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           CF = round(quantile(d5$SMD_CF,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           NT = round(quantile(d5$SMD_NT,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2),
           COMBI = round(quantile(d5$SMD_COMBI,c(0.05,0.25,0.5,0.75,0.95),na.rm=T),2))

table(pmin(d5$n_sp_sw_crit, d5$n_sp_gw_crit) <= d5$n_sp_ref)*100/nrow(d5)
table(pmin(d5$n_sp_sw_crit, d5$n_sp_gw_crit) <= d5$n_sp_ref - (1 - d5$SMD_COMBI) * 0.01 *d5$nin )*100/nrow(d5)

# save predicted output
d6 <- copy(d5)
setnames(d6,tolower(gsub('SMD_','nue_',colnames(d6))))
saveRDS(d6,paste0(floc,'ncu_nue_meas.rds'))
