# script to load data from previous INTEGRATOR database
# the databases have been prepared as being part of the PhD of Maddy Young
# script prepared by Gerard Ros, gerard.ros@wur.nl

# settings
rm(list=ls())

# require packages
require(data.table); require(readxl)

# --- load existing databases INTEGRATOR ----
  
  # load external datasets not stored at github given its size
  floc <- 'D:/DATA/17 nutribudget/'
  
  # load the database which contains NUTS levels per NCU (file prepared in make_ncu_nuts.r)
  ncu_nuts_id <- readRDS(paste0(floc,'ncu_nuts_3035.rds'))
  
  # select unique combo of ncu-nuts0-nuts1-nuts2 to join NCUs to NUTS2 level (and not NUTS3)
  ncu_nuts_id <- unique(ncu_nuts_id[,.(ncu,NUTS0,NUTS1,NUTS2)])
  
  # read excel with NCU data
  d1 <- as.data.table(readxl::read_xlsx(paste0(floc,'ncu_data_crop_area.xlsx'),skip=1))

  # setnames with ****yield target and reference switched back correctly****
  setnames(d1,c('country','crop','texture','ncu','area_ncu',
                'yield_target','yield_ref','soc_ref','soc_target','n_sp_ref','n_sp_sw_crit','n_sp_gw_crit',
                'density','cn','clay','ph','n_fert','n_man','n_fix','n_dep','n_covercrop',
                'nh3_man','nh3_fert','n_leach_frac'))

  # rename code 40 (extensive grassland) to code 33 (grassland)
  d1[crop==40, crop := 33]

  # crop naming
  d1.crop <- as.data.table(readxl::read_xls(paste0(floc,'CAPRI_CROP_DEF.xls')))
  setnames(d1.crop, c('crop_id','crop_code','crop_name'))

  # merge with d1 so that the ncu data also contains crop name (and code)
  d1 <- merge(d1, d1.crop, by.x = 'crop', by.y = 'crop_id',all.x = TRUE)

  # merge with covariates
  d2 <- fread(paste0(floc,'ncu_covariates.csv'))
  d1 <- merge(d1,d2,by = c('ncu','crop_code','country'),all.x=TRUE)

  
# --- load existing databases EUROSTAT ----
  
  # load the excel table with country codes for the EU database (file prepared by Maddy)
  # select first column with names (code and full name are together)
  nuts_name_code <- as.data.table(read_xls(paste0(floc,'eurostat_codes.xls'),range = 'A1:A320'))
  # rename column
  setnames(nuts_name_code,'country_description')
  # split/extract NUTS code from the country description (full name)
  nuts_name_code[,nuts_code := tstrsplit(country_description,' - ',keep = 1)]
  # split/extract NUTS name from country description
  nuts_name_code[,nuts_name := sub(".*? - ", "", country_description)]
  # two small countries are same size as NUTS2 level, resulting in 2 identical rows => remove NUTS0 code
  nuts_name_code <- nuts_name_code[!nuts_code %in% c('LU','MT')]
  
# --- get areas per land use, crop rotation
  
  # get areas for arable land, permanent grassland and permanent crops per NUTS
  d2.arable <- as.data.table(read_xlsx(paste0(floc,'eurostat_till_cov_rot.xlsx'),sheet='Arable land area',skip=4, range='A5:G335'))
  # change column names
  setnames(d2.arable,c('nuts_name','ar_farm','ar_ag_total','ar_arable','ar_perm_grass','ar_perm_crops','ar_gardens'))
  # remove ':' and other special characters and replace with NA
  d2.arable <- d2.arable[,lapply(.SD,function(x) gsub('\\:|Special value|not available',NA,x)),.SDcols = colnames(d2.arable)]
  cols <- colnames(d2.arable)[grepl('ar_',colnames(d2.arable))]
  d2.arable <- d2.arable[,c(cols) := lapply(.SD,function(x) pmax(0,as.numeric(x),na.rm=TRUE)),.SDcols = cols]
  # remove rows that have no nuts name
  d2.arable <- d2.arable[!is.na(nuts_name)]
  # remove full duplicates
  d2.arable <- unique(d2.arable)
  # update the EU databases with NUTS code
  d2.arable <- merge(d2.arable,nuts_name_code,by ='nuts_name',all.x = TRUE)
  # remove the rows that have totals on NUTS0 level (countries)
  d2.arable <- d2.arable[nchar(nuts_code)>2]
  # remove data where the total area is missing
  d2.arable <- d2.arable[!(ar_ag_total==0)]
  #add garden area to permanent crops and remove gardens col
  d2.arable <- d2.arable[,ar_perm_crops := ar_perm_crops + ar_gardens]
  d2.arable[,c('ar_gardens','ar_farm') := NULL]

  # get areas for crop rotation
    
  # read crop rotation (TOTAL = "ARABLE LAND AREA" PER NUTS2)
  d2.rot <- as.data.table(read_xlsx(paste0(floc,'eurostat_till_cov_rot.xlsx'),sheet='Crop rotation',skip=4))
  # change column names
  setnames(d2.rot,c('nuts_name','area_arable_rot','area_rot_cont','area_rot_25','area_rot_50','area_rot_75','area_rot_100','area_rot_na'))
  # remove ':' and other special characters and replace with NA
  d2.rot <- d2.rot[,lapply(.SD,function(x) gsub('\\:|Special value|not available',NA,x)),.SDcols = colnames(d2.rot)]
  # remove rows that have no nuts name
  d2.rot <- d2.rot[!is.na(nuts_name)]
  # remove full duplicates
  d2.rot <- unique(d2.rot)
  # update the EU databases with NUTS code
  d2.rot <- merge(d2.rot,nuts_name_code,by ='nuts_name',all.x = TRUE)
  # remove the rows that have totals on NUTS0 level (countries)
  d2.rot <- d2.rot[nchar(nuts_code)>2]
  # remove data where the total area is missing
  d2.rot <- d2.rot[!is.na(area_arable_rot)]
  # replace NA with change between total and sum
  d2.rot.melt <- melt(d2.rot,id.vars = c('nuts_name','nuts_code','area_arable_rot','country_description'),
                      variable.name = 'crop_rot')
  #(one row for each rotation category)
  d2.rot.melt[,value := as.numeric(value)]
  d2.rot.melt[,area_arable_rot := as.numeric(area_arable_rot)]
  #(format numeric)
  d2.rot.melt[,value2 := pmax(0,(area_arable_rot - sum(value,na.rm=T)) / sum(is.na(value))),by = nuts_code]
  #(grouping by nuts2 code, divide total area remaining from filled categories by the number of empty categories)
  d2.rot.melt[is.na(value), value := value2]
  #(assign this to each empty category)
  d2.rot <- dcast(d2.rot.melt,nuts_code + nuts_name + area_arable_rot ~ crop_rot, value.var = 'value')
  
  # read soil cover and residues
  d2.cover <- as.data.table(read_xlsx(paste0(floc,'eurostat_till_cov_rot.xlsx'),sheet='Soil cover and residues',skip=4))
  # change column names
  setnames(d2.cover,c('nuts_name','area_arable_cov','area_cov_excl','area_cov_winter','area_cov_cover','area_cov_per',
                      'area_cov_res', 'area_cov_bare'))
  # remove ':' and other special characters and replace with NA
  d2.cover <- d2.cover[,lapply(.SD,function(x) gsub('\\:|Special value|not available',NA,x)),.SDcols = colnames(d2.cover)]
  # convert some columns to numeric, ans set to 0 when NA
  cols <- colnames(d2.cover)[grepl('area_',colnames(d2.cover))]
  d2.cover <- d2.cover[,c(cols) := lapply(.SD,function(x) pmax(0,as.numeric(x),na.rm=TRUE)),.SDcols = cols]
  # remove rows that have no nuts name
  d2.cover <- d2.cover[!is.na(nuts_name)]
  # remove full duplicates
  d2.cover <- unique(d2.cover)
  # update the EU databases with NUTS code
  d2.cover <- merge(d2.cover,nuts_name_code,by ='nuts_name',all.x = TRUE)
  # remove the rows that have totals on NUTS0 level
  d2.cover <- d2.cover[nchar(nuts_code)>2]
  # add totals when input is available
  d2.cover[area_arable_cov == 0, area_arable_cov := area_cov_excl + area_cov_winter + area_cov_cover + area_cov_bare + area_cov_per + area_cov_res]
  # remove data where the total area is missing
  d2.cover <- d2.cover[area_arable_cov > 0]

  # area per tillage
  
  # read tillage
  d2.tillage <- as.data.table(read_xlsx(paste0(floc,'eurostat_till_cov_rot.xlsx'),sheet='Tillage Practices',skip=4))
  # change column names
  setnames(d2.tillage,c('nuts_name','area_arable_till','area_till_excl','area_till_conv','area_till_conserv','area_till_no'))
  # remove ':' and other special characters and replace with NA
  d2.tillage <- d2.tillage[,lapply(.SD,function(x) gsub('\\:|Special value|not available',NA,x)),.SDcols = colnames(d2.tillage)]
  # convert some columns to numeric, and set to 0 when NA
  cols <- colnames(d2.tillage)[grepl('area_',colnames(d2.tillage))]
  d2.tillage <- d2.tillage[,c(cols) := lapply(.SD,function(x) pmax(0,as.numeric(x),na.rm=TRUE)),.SDcols = cols]
  # remove rows that have no nuts name
  d2.tillage <- d2.tillage[!is.na(nuts_name)]
  # remove full duplicates
  d2.tillage <- unique(d2.tillage)
  # update the EU databases with NUTS code
  d2.tillage <- merge(d2.tillage,nuts_name_code,by ='nuts_name',all.x = TRUE)
  # remove the rows that have totals on NUTS0 level
  d2.tillage <- d2.tillage[nchar(nuts_code)>2]
  # add totals when total is missing but individual categories are available
  d2.tillage[area_arable_till == 0,area_arable_till := area_till_excl + area_till_conv + area_till_conserv + area_till_no]
  # remove data where the total area is missing
  d2.tillage <- d2.tillage[area_arable_till > 0]

# make likelyhood for measures per crop rotation
# to be used in the connection from INTEGRATOR with Eurostat
  
  # prepare likelihood table for all cropping practices (later tillage depends on cropping)
  # grepl function searches for the character strings (using lower case for all) and adds a likelihood col
  lkh.crop <- data.table(crop_name = unique(d1$crop_name))
  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.rot := 0.4]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.rot := 0.8]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.rot := 0.8]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.rot := 0.8]
  #temporary industrial crops
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.rot := 0.5]
  #permanent industrial crops
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.rot := 0.01]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$|vegeta',tolower(crop_name)), lh.rot := 0.01]

  # prepare likelihood table for soil cover

  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.win := 0.99]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.win := 0.7]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.win := 0.7]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.win := 0.4]
  #temporary industrial
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.win := 0.8]
  #permanent industrial
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.win := 0.7]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$|vegeta',tolower(crop_name)), lh.win := 0.01]

  # prepare likelihood table for cover or intermediate crop

  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.cov := 0.99]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.cov := 0.8]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.cov := 0.9]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.cov := 0.5]
  #temporary industrial
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.cov := 0.6]
  #permanent industrial
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.cov := 0.01]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$|vegeta',tolower(crop_name)), lh.cov := 0.01]

  # prepare likelihood table for residues

  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.res := 0.99]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.res := 0.9]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.res := 0.7]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.res := 0.7]
  #temporary industrial
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.res := 0.7]
  #permanent industrial
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.res := 0.6]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$|vegeta',tolower(crop_name)), lh.res := 0.01]


  # prepare likelihood table for perennials

  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.per := 0.99]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.per := 0.01]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.per := 0.01]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.per := 0.01]
  #temporary industrial
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.per := 0.01]
  #permanent industrial
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.per := 0.01]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$|vegeta',tolower(crop_name)), lh.per := 0.99]

  # prepare likelihood table for bare soil

  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.bar := 0.99]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.bar := 0.7]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.bar := 0.7]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.bar := 0.7]
  #temporary industrial
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.bar := 0.7]
  #permanent industrial
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.bar := 0.7]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$|vegeta',tolower(crop_name)), lh.bar := 0.7]

  # prepare likelihood table for under glass or protective cover

  #small grain cereals
  lkh.crop[grepl('wheat|bar|rey|cere|oat|rice',tolower(crop_name)), lh.exc := 0.99]
  #maize
  lkh.crop[grepl('maize',tolower(crop_name)), lh.exc := 0.01]
  #legumes
  lkh.crop[grepl('soya|puls',tolower(crop_name)), lh.exc := 0.01]
  #tubers
  lkh.crop[grepl('potat|sugar',tolower(crop_name)), lh.exc := 0.01]
  #temporary industrial
  lkh.crop[grepl('rape|cotton|oil|tobac|industrial|sunf',tolower(crop_name)), lh.exc := 0.01]
  #permanent industrial
  lkh.crop[grepl('industrial',tolower(crop_name)), lh.exc := 0.01]
  #temporary grassland, permanent grassland, permanent crops, nursery
  lkh.crop[grepl('oliv|orang|wine|fruit|grap|tempor|gras$',tolower(crop_name)), lh.exc := 0.01]
  lkh.crop[grepl('vegeta',tolower(crop_name)), lh.exc := 0.5]

  # prepare table likelihood for conventional tillage
  lkh.till <- data.table(crop_name = unique(d1$crop_name))
  lkh.till[grepl('wheat|bar|rey|cere',tolower(crop_name)), lh := 0.8]
  lkh.till[grepl('potat|maize|soya|sugar|vegeta|oat|industrial',tolower(crop_name)), lh := 1.0]
  lkh.till[grepl('puls|rice|sunf',tolower(crop_name)), lh := 0.4]
  lkh.till[grepl('oliv|orang|wine|fruit|oil|tobac|grap',tolower(crop_name)), lh := 0]
  lkh.till[grepl('grassland|gras$',tolower(crop_name)), lh := 0]
  lkh.till[grepl('rape|cotton',tolower(crop_name)), lh := 0.2]

# ---- integration INTEGRATOR with Eurostat data -----

  # merge NCU data with given NUTS based on country and ncu
  d3 <- merge.data.table(d1,
                         ncu_nuts_id,
                         by.x = c('ncu','country'),
                         by.y = c('ncu','NUTS0'),
                         all.x = TRUE ,
                         all.y = FALSE)
  
  # add a country dependent correction factor for area per NCU
  d3[,ncu_areacountry := sum(area_ncu),by=.(ncu,NUTS2)] #total of each ncu type matched at NUTS2 level
  d3[,ncu_areacountrycf := ncu_areacountry / sum(area_ncu),by=.(ncu,country)] #divide by total of ncu type at country level
  
  # adapt the area of the ncu
  d3[, area_ncu := ncu_areacountrycf * area_ncu]
  
  # join with EU database
  d3 <- merge(d3,d2.arable,by.x = 'NUTS2',by.y='nuts_code',all.x=TRUE)
  d3 <- merge(d3,d2.rot[,mget(colnames(d2.rot)[!grepl('nuts_name|country_',colnames(d2.rot))])],by.x = 'NUTS2',by.y='nuts_code',all.x=TRUE)
  d3 <- merge(d3,d2.tillage[,mget(colnames(d2.tillage)[!grepl('nuts_name|country_',colnames(d2.tillage))])],by.x = 'NUTS2',by.y='nuts_code',all.x=TRUE)
  d3 <- merge(d3,d2.cover[,mget(colnames(d2.cover)[!grepl('nuts_name|country_',colnames(d2.cover))])],by.x = 'NUTS2',by.y='nuts_code',all.x=TRUE)
  
  # about 67 ha (< 0.01%) has no NUTS code, so delete these ncu's (~400 rows)  #309,961
  d3 <- d3[!is.na(NUTS2)]
  
  # add total area of integrator per NUTS2 area
  d3[,area_int_nuts := sum(area_ncu*100,na.rm = TRUE), by =.(NUTS2)]
  
  # join likelihood continuous cropping to the ncu database
  d3 <- merge(d3,lkh.crop,by='crop_name',all.x=TRUE)
  
  # helper function to connect datasets given the likelihood
  dsRA_new <- function(area_ncu,area_nuts2,lkh, nuts2,area_nuts2_all = NULL){
    
    # make local copy of the inputs
    dt <- data.table(area_ncu, area_nuts2,lkh, nuts2,area_nuts2_all)
    
    # assume that total area for the measure equals the measure to be reallocated
    if(is.null(area_nuts2_all)){
      dt[,area_nuts2_all := area_nuts2]
    }
    
    # convert unit for integrator NCU (km2) to match eurostat (ha)
    dt[,area_ncu := area_ncu * 100]
    
    # add unique id
    dt[, id := .I]
    
    # estimate the total area in ncu on nutslevel
    dt[,a_0 := sum(area_ncu,na.rm = T),by=nuts2]
    
    # estimate total likely area on NUTS2 level given NCU crop area and likelihood for measure X
    dt[,a_1 := sum(lkh * area_ncu, na.rm = TRUE), by = nuts2]
    
    # correct this proportial to the area that is likely used for management measure X
    # the area of some NCUs already maximized/filled up, so cannot go over the total crop area
    # if likely area is less than NCU area, it saves this (see below)
    # gives an NCU likely area which is proportional to total area at the NUTS2 level
    dt[,a_2 := pmin(area_ncu,lkh * area_ncu * area_nuts2_all / a_1,na.rm = TRUE), by = nuts2]
    
    # how much NUTS2 area leftover should still be reallocated for management measure X
    dt[,a_3 := area_nuts2_all - sum(a_2), by = nuts2]
    
    # how much NCU area is still left within crops that have a likelihood for measure X
    # calculates area where more may be applied, if lkh is 0 then assigns 0 value
    dt[,a_4 := fifelse(lkh > 0,(area_ncu - a_2),0), by = nuts2]
    
    # redistribute the area that need to be reallocated to the area that is left
    # calculates exact proportion of NUTS2 area to go to each; if lkh is 0 then a_5 becomes 0
    dt[,a_5 := a_3 * a_4/sum(a_4),by = nuts2]
    
    # estimate total area of measure X to add to each NCU; if no more room then 0 added
    dt[,a_6 := a_2 + pmax(0,a_5,na.rm=T)]
    
    # estimate the fraction of measure x
    dt[,afin := a_6 * area_nuts2 / area_nuts2_all]
    
    # ensure correct order
    setorder(dt,id)
    
    # extract value
    value <- dt[,afin]
    
    # return value
    return(value)
    
  }
  
  
  ############################## calculate the area continuous cropping ##############################
  
  # estimate the NCU crop area used for continuous cropping (opposite of rotation, so likelihood is 1 minus likelihood)
  d3[, rot_area_cont := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_rot_cont,lkh = 1-lh.rot, nuts2 = NUTS2)]
  
  # estimate the NCU crop area used for 0-25% crop rotation (likelihood is likelihood for rotation)
  d3[, rot_area_cr25 := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_rot_25, lkh = lh.rot, nuts2 = NUTS2)]
  
  # estimate the NCU crop area used for 25-50% crop rotation (likelihood is likelihood for rotation)
  d3[, rot_area_cr50 := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_rot_50, lkh = lh.rot, nuts2 = NUTS2)]
  
  # estimate the NCU crop area used for 50-75% crop rotation (likelihood is likelihood for rotation)
  d3[, rot_area_cr75 := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_rot_75, lkh = lh.rot, nuts2 = NUTS2)]
  
  # estimate the NCU crop area used for >75% crop rotation (likelihood is likelihood for rotation)
  d3[, rot_area_cr100 := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_rot_100, lkh = lh.rot, nuts2 = NUTS2)]
  
  # what is total NCU crop area already under crop rotation (in ha), summing over different categories
  d3[, rot_area_tot := rot_area_cont * 0 + rot_area_cr25 * 0.125 + rot_area_cr50 * 0.375 +
       rot_area_cr75 * 0.625 + rot_area_cr100 * 0.875]
  
  # remove columns no longer needed
  cols <- c('area_arable_rot','area_rot_cont', 'area_rot_25', 'area_rot_50', 'area_rot_75', 'area_rot_100', 'area_rot_na',
            'rot_area_cont','rot_area_cr25', 'rot_area_cr50', 'rot_area_cr75', 'rot_area_cr100')
  d3[,c(cols) := NULL]
  
  # area of arable land within each crop type where crop rotation may potentially be applied (in ha)
  d3[, area_mon := area_ncu * 100 - rot_area_tot]
  
  # if there is zero applicability of rotation, the full area is assigned monoculture but it is
  # instead fully permanent crops with no potential for conversion
  d3[lh.rot == 0, area_mon := 0]
  
  #cols to check if total area matches up
  d3[,sum_rot := rot_area_tot + area_mon]
  d3[,area_ncu_ha := area_ncu*100]
  d3[,test_rot :=  area_ncu_ha - sum_rot]
  
  # what is the potential area for more cereals in crop rotation / diversification of crop rotation plan (cr)
  d3[,parea.cr := area_mon]
  
  # remove variables not needed any more
  d3[,c('sum_rot','test_rot','lh.rot','area_mon','rot_area_tot') := NULL]
  
  
  ############################## calculate the soil cover area ######################################
  
  # corrects soil cover columns so that the totals do not go over total cover given by EUROSTAT?
  d3[,area_arable_cov_sum := area_cov_bare + area_cov_cover + area_cov_per + area_cov_res + area_cov_winter +area_cov_excl]
  cols <- colnames(d3)[grepl('^area_cov', colnames(d3))]
  d3[,c(cols) := lapply(.SD,function(x) x * area_arable_cov_sum / area_arable_cov), .SDcols = cols]
  
  # normal harvested winter cover
  d3[, area_win := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_cov_winter,
                            lkh = lh.win, nuts2 = NUTS2,
                            area_nuts2_all = area_arable_cov_sum )]
  
  # cover/intermediate crop
  d3[, area_cov := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_cov_cover,
                            lkh = lh.cov, nuts2 = NUTS2,
                            area_nuts2_all = area_arable_cov_sum)]
  
  # plant residue cover
  d3[, area_res := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_cov_res,
                            lkh = lh.res, nuts2 = NUTS2,
                            area_nuts2_all = area_arable_cov_sum)]
  
  # multi-annuals
  d3[, area_per := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_cov_per,
                            lkh = lh.per, nuts2 = NUTS2,
                            area_nuts2_all = area_arable_cov_sum)]
  
  # bare soil
  d3[, area_bar := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_cov_bare,
                            lkh = lh.bar, nuts2 = NUTS2,
                            area_nuts2_all = area_arable_cov_sum)]
  
  # excluded/non applicable area (permanent cover)
  d3[, area_exc := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_cov_excl,
                            lkh = lh.exc, nuts2 = NUTS2,
                            area_nuts2_all = area_arable_cov_sum)]
  
  # area non cropland
  d3[, area_ncl := area_ncu_ha - (area_win + area_cov + area_res + area_per + area_bar + area_exc)]
  d3[, sum_cov := area_ncl + area_win + area_cov + area_res + area_per + area_bar + area_exc]
  
  # what is the potential area for cover crops application (cc) and area for crop residue application
  d3[,parea.cc := fifelse(lh.cov < 0.02,0, area_ncl + area_exc + area_bar)]
  d3[,parea.cres := fifelse(lh.res < 0.02,0, area_ncl + area_exc + area_bar + area_per + area_win)]
  
  # remove variables not needed any more
  d3[,c('lh.win','lh.cov','lh.res','lh.per','lh.bar','lh.exc','area_arable_cov_sum',
        'area_arable_cov', 'area_cov_excl', 'area_cov_winter', 'area_cov_cover', 'area_cov_per', 'area_cov_res', 'area_cov_bare') := NULL]
  
  # join likelyhood continous cropping to the ncu database
  d4 <- merge(d3,lkh.till[,.(crop_name,lh.ctill = lh)],by='crop_name',all.x=TRUE)
  
  # reduced and no till are more likely where conservation practices occur (cover and residues)
  d4[, lh.ntill := pmin(fifelse((area_win + area_per + area_bar) / sum_cov > 0.5,lh.ctill + 0.3,lh.ctill),1)]
  d4[, lh.rtill := pmin(fifelse((area_win + area_per) / sum_cov > 0.5,lh.ctill + 0.3,lh.ctill),1)]
  
  # area conventional tilled soil
  d4[, area_conv_till := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_till_conv,
                                  lkh = lh.ctill, nuts2 = NUTS2, area_nuts2_all = area_arable_till)]
  
  # area conservative tilled soil
  d4[, area_conserv_till := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_till_conserv,
                                     lkh = lh.ctill, nuts2 = NUTS2, area_nuts2_all = area_arable_till)]
  
  # area no-till soil
  d4[, area_no_till := dsRA_new(area_ncu = area_ncu, area_nuts2 = area_till_no,
                                lkh = lh.ntill, nuts2 = NUTS2, area_nuts2_all = area_arable_till)]
  
  # what is the potential area that can be used for tillage practices
  d4[,parea.rtct := fifelse(lh.ctill == 0, 0, area_conv_till)]
  d4[,parea.ntct := fifelse(lh.ctill == 0 | grepl('potato|sugar',tolower(crop_name)),0,area_conv_till+area_conserv_till)]
  
  # remove variables not needed anymore
  d4[,c('lh.ctill','lh.ntill','lh.rtill','area_arable_till', 'area_till_excl', 'area_till_conv' ,'area_till_conserv', 'area_till_no',
        'area_conv_till','area_conserv_till',
        'area_win', 'area_cov' ,  'area_res',  'area_per', 'area_bar','area_exc', 'area_ncl') := NULL]
  
  # how much organic manure is available on grassland within a NUTS region
  # n_man = N inputs from manure
  d4[, n_man_nuts := fifelse(grepl('grassland|Gras',crop_name),n_man,0)] #select manure N only from grassland
  d4[, n_man_nuts := sum(n_man_nuts,na.rm = TRUE),by = NUTS2] #sum over NUTS2 region
  
  # this amount might be adapted to the volume produced in stables
  # so correct for grazing (see e.g. Eurostat); a relative fraction of total applied
  
  # what is the total effective N input on arable land (available manure and mineral on all NON-grass crops)
  # assume only 60% of manure N available may be applied
  d4[, neff := fifelse(!grepl('grassland|Gras',crop_name),n_man * 0.6 + n_fert,0)]
  
  # add the extra manure from grassland to all cropland relative to the total N effective dose
  # so that manure at NUTS2 level is proportionally redistributed to each NCU, and none is added to grassland
  d4[, n_man_ncu := n_man_nuts * neff / sum(neff,na.rm = TRUE), by = NUTS2]
  
  ############################## calculate max change in SOC from organic inputs ######################################
  
  # assume that 50% of manure is slurry and 50% is solid manure
  # slurry has 50 kg EOS (effect organic material) per ton and solid manure has 109 kg EOS per ton
  # see also humification rations and CN ratios (12 and 17)
  # https://www.handboekbodemenbemesting.nl/nl/handboekbodemenbemesting/Handeling/Organische-stofbeheer/Organische-stofbalans/Kengetallen-organische-stof.htm
  # this is the max change in SOC allowed due to manure applicacation (kg C/ha)
  d4[, c_man_ncu := 0.5 * n_man_ncu * 12 * 0.7 * (2.1/4.0) + 0.5 * n_man_ncu * 17 * 0.7 * (6.6/7.7)]  #solid manure
  
  # area for combination of fertilizer vs inorganic fertilizers
  
  # situation NIFNOF: all cases with less than 40 kg N / ha
  d4[,parea.nifnof := fifelse(!grepl('gras',tolower(crop_name)) & n_man < 40 & n_fert < 40, area_ncu * 100, 0)]
  
  # medium fertilized, low animal N-input (MIFNOF): all cases up to 100 kg N / ha and less than 40 kg N/ha from animal manure
  d4[,parea.mifnof := fifelse(!grepl('gras',tolower(crop_name)) & n_man < 40 & n_fert >= 40 & n_fert < 100, area_ncu * 100, 0)]
  
  # highly fertilized, low animal N-input (HIFNOF): all cases with more than 100 kg N /ha and less than 40 kg N/ha from animal manure
  d4[,parea.hifnof := fifelse(!grepl('gras',tolower(crop_name)) & n_man < 40 & n_fert >= 100, area_ncu * 100, 0)]
  
  # medium fertilized, high animal N-input (MIFHOF): all cases up to 100 kg N / ha and more than 40 kg N/ha from animal manure
  d4[,parea.mifhof := fifelse(!grepl('gras',tolower(crop_name)) & n_man >= 40 & n_fert >= 40 & n_fert < 100, area_ncu * 100, 0)]
  
  # highly fertilized, high animal N-input (HIFHOF): all cases with more than 100 kg N /ha and more than 40 kg N/ha from animal manure
  d4[,parea.hifhof := fifelse(!grepl('gras',tolower(crop_name)) & n_man >= 40 & n_fert >= 100, area_ncu * 100, 0)]
  
  # add type MA model
  d4[parea.nifnof > 0,fert_type := 'nifnof']
  d4[parea.mifnof > 0,fert_type := 'mifnof']
  d4[parea.hifnof > 0,fert_type := 'hifnof']
  d4[parea.mifhof > 0,fert_type := 'mifhof']
  d4[parea.hifhof > 0,fert_type := 'hifhof']
  d4[is.na(fert_type),fert_type := 'none']
  
  # remove
  d4[,c('country_description','nuts_name') := NULL]
  
  #add letters to covariate names so that low/medium/high categories are distinguished differently for each variable
  d4[,cov_fert := paste0('f',cov_fert)]
  d4[,cov_soc := paste0('c',cov_soc)]
  
  # save file as csv
  fwrite(d4,paste0(floc,'db_final_europe.csv'))
  
  
  
  