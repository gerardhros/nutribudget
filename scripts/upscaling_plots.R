# plots needed for upscaling, deliverable D1.5
# plot the current situation for SOC, pH, yield, and NUE

# require packages

require(ggplot2); require(sf); require(rnaturalearth)
require(rnaturalearthdata); require(terra); require(data.table)

# set theme
theme_set(theme_bw())

# load external datasets not stored at github given its size
floc <- 'D:/DATA/17 nutribudget/'

# load in the database with values for 14 crops per NCU
d1 <- fread(paste0(floc,'db_final_europe.csv'))

# laod in NUTS shapefile
s.nuts <- st_read(paste0(floc,'eu_nuts.gpkg'),layer='eu_nuts')

# calculate the mean properties for a series of columns
cols <- c('area_ncu','yield_target', 'yield_ref', 'soc_ref','soc_target', 'n_sp_ref','n_sp_sw_crit', 'n_sp_gw_crit',
          'ph', 'n_fert','n_man','n_dep')
d2 <- d1[,mget(c(cols,'ncu'))]
d2 <- d2[,lapply(.SD,function(x) weighted.mean(x,w=area_ncu)),by='ncu']

# get the raster to plot
r1 <- terra::rast(paste0(floc,'gncu2010_ext.asc'))
terra::crs(r1) <- 'epsg:3035'
r1 <- terra::project(r1,'epsg:4326',method='near')

# convert to data.frame
r1.p <- as.data.frame(r1,xy=TRUE)
r1.p <- as.data.table(r1.p)

# join/merge d1 with r1.p
r.ncu <- merge(r1.p, d2, by.x = 'gncu2010_ext', by.y = 'ncu')

# get summary stats
sstat <- function(x,r){round(quantile(x,c(0.05,0.5,0.95),na.rm=T),r)}

# ---- plot data on EU level, current -----
p1 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= yield_ref/1000)) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'yield (t/ha)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Crop yield") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p1, filename = paste0(floc,'up_yieldref.jpg'),width = 12,height=12,units='cm')

p2 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= pmin(20,soc_ref))) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'SOC (g/kg)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Soil organic carbon") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p2, filename = paste0(floc,'up_socref.jpg'),width = 12,height=12,units='cm')

p3 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= pmax(0,n_sp_ref))) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'N surplus\n(kg/ha)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Nitrogen surplus agriculture") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p3, filename = paste0(floc,'up_nspref.jpg'),width = 12,height=12,units='cm')

p4 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= ph)) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'soil pH')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Soil pH") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p4, filename = paste0(floc,'up_soilph.jpg'),width = 12,height=12,units='cm')

list(yield = sstat(r.ncu$yield_ref,0),
     soc = sstat(r.ncu$soc_ref,1),
     nsp = sstat(r.ncu$n_sp_ref,0),
     ph = sstat(r.ncu$ph,1))

# ---- plot data on EU level, distance to target -----
p1 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= pmin(1,yield_ref/yield_target))) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'dtt-yield (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Crop yield, distance-to-target") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p1, filename = paste0(floc,'up_dtt_yield.jpg'),width = 12,height=12,units='cm')

p2 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= pmin(4,soc_ref/soc_target))) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'dtt-SOC (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("SOC, distance-to-target") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p2, filename = paste0(floc,'up_dtt_soc.jpg'),width = 12,height=12,units='cm')

p3 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= pmin(10,pmax(0,n_sp_ref/n_sp_sw_crit)))) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'dtt-N surplus (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Nitrogen surplus, distance-to-target") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p3, filename = paste0(floc,'up_dtt_nsp.jpg'),width = 12,height=12,units='cm')

p4 <- ggplot() +
      geom_sf(data=s.nuts,color='black',fill=NA,show.legend = FALSE)+
      geom_tile(data = r.ncu,aes(x=x,y=y,fill= ph/5.5)) +
      scale_fill_viridis_c(direction=-1)+
      theme(legend.position.inside = c(0.1,0.8))+
      labs(fill = 'dtt-soil pH (-)')+
      xlab("Longitude") + ylab("Latitude") +
      ggtitle("Soil pH, distance-to-target") +
      coord_sf(crs = 4326) + theme_bw()
ggsave(plot = p4, filename = paste0(floc,'up_dtt_ph.jpg'),width = 12,height=12,units='cm')

list(yield = sstat(r.ncu$yield_ref/r.ncu$yield_target,1),
     soc = sstat(r.ncu$soc_ref/r.ncu$soc_target,1),
     nsp = sstat(r.ncu$n_sp_ref/r.ncu$n_sp_sw_crit,1),
     ph = sstat(r.ncu$ph/5.5,1))
