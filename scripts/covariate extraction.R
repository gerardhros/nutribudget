# covariate extraction for meta-analysis
# may 2022

# clear environment
rm(list=ls())

# require packages
require(readxl);require(terra);require(data.table);require(sf)
require(stringr); require(foreign)

# variable selection: options yield, ph, soc
var = 'soc' 

# load in the data with location
if(var=='yield'){
  d1 <- fread('data/yield_lonlat_unique.csv')
  setnames(d1,c('lon_deg','lat_deg'),c('lon','lat'))
} else if(var=='ph'){
  d1 <- fread('data/ph_lonlat_unique.csv')
} else if(var=='soc'){
  d1 <- fread('data/soc_lonlat_unique.csv')
}

# convert to a spatial object
s1 <- st_as_sf(d1,coords = c('lon','lat'),crs = 4326)
s1 <- vect(s1)

# read in dbf file for metzger climatic regions
s2 <- foreign::read.dbf('D:/DATA/03 metzger/GenS_v3.dbf')

# what rasters are available
# downloaded via QGIS for ISRIC, 0.5 degrees resolution, https://maps.isric.org/
# downloaded via CRU, https://catalogue.ceda.ac.uk/uuid/89e1e34ec3554dc98594a5732622bce9
# downloaded via https://datashare.ed.ac.uk/handle/10283/3089

# read in the rasters via hard drive
r1 <- list.files('D:/DATA/01 soil',pattern = 'tif|nc',full.names = T)
r1 <- r1[!grepl('stack',r1)]
r2 <- list.files('D:/DATA/02 climate',pattern = 'tif|nc',full.names = T)
r3 <- list.files('D:/DATA/03 metzger',pattern = 'tif|nc',full.names = T)

# read in the raster files and convert to spatrasters
isric <- sds(r1)
climate <- sds(r2)
metzger <- rast(r3)
isric <- rast(isric)
climate <- rast(climate)

# --- extract isric data ----

  # update names of isric raster to avoid duplication in names
  names(isric) <- str_split_fixed(names(isric),"_isric_",2)[,2]

  # extract data for the spatial objects
  d1.isric <- terra::extract(x = isric, y = s1)
  d2.isric <- terra::extract(x = isric, y = buffer(s1,width = 10000), fun = mean, na.rm=T)
  d3.isric <- terra::extract(x = isric, y = buffer(s1,width = 20000), fun = mean, na.rm=T)
  d4.isric <- terra::extract(x = isric, y = buffer(s1,width = 30000), fun = mean, na.rm=T)
  d5.isric <- terra::extract(x = isric, y = buffer(s1,width = 50000), fun = mean, na.rm=T)

  # convert to data.table to facilitatie re-arranging
  setDT(d1.isric);setDT(d2.isric);setDT(d3.isric);setDT(d4.isric);setDT(d5.isric)

    # function to adapt colnames
    acn <- function(x,var='or'){c('ID',paste0(var,'_',gsub('isric_|_mean_|','',x[-1])))}

    # adapt colnames
    setnames(d1.isric,acn(colnames(d1.isric)))
    setnames(d2.isric,acn(colnames(d2.isric),'e1'))
    setnames(d3.isric,acn(colnames(d3.isric),'e2'))
    setnames(d4.isric,acn(colnames(d4.isric),'e3'))
    setnames(d5.isric,acn(colnames(d5.isric),'e4'))

  c1.isric <- merge(d1.isric,d2.isric,by = "ID")
  c1.isric <- merge(c1.isric,d3.isric,by = "ID")
  c1.isric <- merge(c1.isric,d4.isric,by = "ID")
  c1.isric <- merge(c1.isric,d5.isric,by = "ID")

  c1.isric <- melt(c1.isric,id = 'ID',
                   measure=patterns("or_", "e1_","e2_","e3_","e4_"),
                   variable.factor = FALSE,
                   value.name = c("or", "e1","e2","e3","e4"))
  c1.isric[,variable := sort(names(isric))[as.integer(variable)]]
  c1.isric[,value := as.numeric(or)]
  c1.isric[or < 1 & e1 > 1,value := e1]
  c1.isric[or < 1 & e2 > 1,value := e2]
  c1.isric[or < 1 & e3 > 1,value := e3]
  c1.isric[or < 1 & e4 > 1,value := e4]

  c2.isric <- dcast(c1.isric,ID~variable, value.var = 'value')


# --- extract climate data ----


  # extract climate data nc files
  d1.climate <- terra::extract(x = climate, y = s1)

  # convert to data.table
  d1.climate <- as.data.table(d1.climate)

  # rearrange data
  d2.climate <- melt(d1.climate,id.vars = 'ID', variable.name = 'variable')
  d2.climate <- d2.climate[!grepl('_stn_',variable)]
  d2.climate[, cvar :=  stringr::str_extract_all(variable,"(?<=[0-9]{4}\\.[0-9]{4}\\.).+(?=\\.dat_)",simplify = T)]
  d2.climate[, years :=  stringr::str_extract_all(variable,"[0-9]{4}\\.[0-9]{4}",simplify = T)]
  d2.climate[, month :=  stringr::str_extract_all(variable,"(?<=[a-z]{3}_)\\d+",simplify = T)]

  # estimate mean global climate properties over 1991-2019
  # temperature in degrees (mean = tmp, max = tmx, min = tmn)
  # potential evaporation in mm/day
  # precipitation in mm/month
  d3.climate <- dcast(d2.climate,ID+years+month~cvar,value.var = 'value')

  # derive the mean and SD per gridcel over period 1991-2019
  c1.climate <- d3.climate[,list(pre_mean = mean(pre),
                                 pre_sd = sd(pre),
                                 tmp_mean = mean(tmp),
                                 tmp_sd = sd(tmp),
                                 pet_mean = mean(pet),
                                 pet_sd = sd(pet)
                                 ),by='ID']
  c2.climate <- copy(c1.climate)

# --- extract metzger data

  # extract for measurement points
  d1.metzger <- terra::extract(x = metzger, y = s1)

  # read metzger decription
  s2.dt <- as.data.table(s2)

  # merge description
  d1.metzger <- as.data.table(d1.metzger)
  c2.metzger <- merge(d1.metzger,s2.dt,by.x = 'gens_v3', by.y = 'GEnS_seq')

  # subset
  c2.metzger <- c2.metzger[,.(ID,GEnZname,GEnZ,GEnS)]

# merge the data files
  d1[,ID := 1:.N]
  dt <- d1[,.(ID,lon,lat)]
  dt <- merge(dt,c2.isric, by = 'ID')
  dt <- merge(dt,c2.climate, by = 'ID')
  dt <- merge(dt,c2.metzger,by='ID')

# save the file
  if(var=='yield'){
    fwrite(dt,'products/240523_covariates_yield.csv')
  } else if(var=='ph'){
    fwrite(dt,'products/240617_covariates_ph.csv')
  } else if(var=='soc'){
    fwrite(dt,'products/240617_covariates_soc.csv')
  }
  
 
