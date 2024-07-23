# covariate extraction for meta-analysis
# may 2022

# clear environment
rm(list=ls())

# require packages
require(readxl);require(terra);require(data.table);require(sf)
require(stringr); require(foreign)

# load external datasets not stored at github given its size
floc <- 'D:/DATA/17 nutribudget/'

# get the raster to plot
r.ncu <- terra::rast(paste0(floc,'gncu2010_ext.asc'))
terra::crs(r.ncu) <- 'epsg:3035'
r.ncu <- terra::project(r.ncu,'epsg:4326',method='near')

# what rasters are available
# downloaded via QGIS for ISRIC, 0.5 degrees resolution, https://maps.isric.org/
# downloaded via CRU, https://catalogue.ceda.ac.uk/uuid/89e1e34ec3554dc98594a5732622bce9
# downloaded via https://datashare.ed.ac.uk/handle/10283/3089

# read in the rasters via hard drive
r1 <- list.files('D:/DATA/01 soil',pattern = 'tif|nc',full.names = T)
r1 <- r1[!grepl('stack',r1)]
r2 <- list.files('D:/DATA/02 climate',pattern = 'tif|nc',full.names = T)
r2 <- r2[grepl('pet|pre|tmp',r2)]
r3 <- list.files('D:/DATA/03 metzger',pattern = 'tif|nc',full.names = T)

# read in the raster files and convert to spatrasters
isric <- sds(r1)
metzger <- rast(r3)
isric <- rast(isric)
pre <- sds(r2[grepl('pre',r2)])
tmp <- sds(r2[grepl('tmp',r2)])
pet <- sds(r2[grepl('pet',r2)])
pre <- rast(pre)
tmp <- rast(tmp)
pet <- rast(pet)

pre.mean <- app(pre,mean)
pre.sd <- app(pre,sd)
tmp.mean <- app(tmp,mean)
tmp.sd <- app(tmp,sd)
pet.mean <- app(pet,mean)
pet.sd <- app(pet,sd)
climate <- c(pre.mean,pre.sd,tmp.mean,tmp.sd,pet.mean,pet.sd)
names(climate) <- c('pre_mean','pre_sd','tmp_mean','tmp_sd','pet_mean','pet_sd')

# update names of isric raster to avoid duplication in names
names(isric) <- str_split_fixed(names(isric),"_isric_",2)[,2]

# --- extract isric data ----

  # reproject isiric data to r.ncu
  r.soil <- terra::resample(isric,r.ncu,method='bilinear')

  # reproject metzger to r.ncu
  r.metz <- terra::resample(metzger,r.ncu,method='near')

  # reproject climate to r.ncu
  r.clim <- terra::resample(climate,r.ncu,method='bilinear')

  # combine in one raster
  r.ncu.cov1 <- c(r.ncu,r.soil,r.metz,r.clim)

  # convert to data.table
  d1.ncu <- as.data.frame(r.ncu.cov1,xy=TRUE)
  d1.ncu <- as.data.table(d1.ncu)
  d1.ncu <- d1.ncu[!is.na(gncu2010_ext)&!is.na(pre_mean)]
  
  # calculate mean per ncu
  d2.ncu <- d1.ncu[,lapply(.SD,mean,na.rm=T),by='gncu2010_ext']
  
# save the file
  fwrite(d2.ncu,'products/240723_covariates_ncu.csv')
  
 
