
#############################################################################
# Extract PLAND from EarthData

setwd('path/to/AphidDecline')

library(raster)
library(rgdal)
library(rgeos)
library(gdalUtils)
proj1 = '+proj=longlat +datum=WGS84 +no_defs'


#####
# Idaho

cropland = raster('./rasters/EarthData/Idaho/GFSAD30NACE_2010_N40W120_001_2017274094200.tif')

aphid_sites = spTransform(readOGR(dsn='./shapefiles/PNW_aphid_sites_50km.shp',layer='PNW_aphid_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1))
aphid_grid = spTransform(readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
aphid_grid2 = aphid_grid[match(unique(aphid_sites@data$grid_id),aphid_grid@data$grid_id),]

pland = data.frame('grid_id'=aphid_grid2@data$grid_id,'Cropland'=-999)
for (i in 1:length(aphid_grid2)){
	cat(i,'\n')
	grid1 = aphid_grid2[i,]
	cropland1 = crop(cropland,grid1)
	cropland2 = as.matrix(cropland1)
	pland1 = length(which(cropland2==2)) / length(cropland2)
	pland[i,1] = grid1@data$grid_id
	pland[i,2] = pland1
}
write.table(pland,'./data/croplandED_50km_IDAHO.txt',sep='\t',row.names=F,quote=F)


#####
# MIDWEST

files = list.files('./rasters/EarthData/Midwest',full.names=T); files = files[grep('.tif',files)]
mosaic_rasters(files,'./rasters/Midwest_cropland_EarthData.tif')
cropland = raster('./rasters/Midwest_cropland_EarthData.tif')

aphid_sites = spTransform(readOGR(dsn='./shapefiles/aphid_sites_50km.shp',layer='aphid_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1))
aphid_grid = spTransform(readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
aphid_grid2 = aphid_grid[match(unique(aphid_sites@data$grid_id),aphid_grid@data$grid_id),]

pland = data.frame('grid_id'=aphid_grid2@data$grid_id,'Cropland'=-999)
for (i in 1:length(aphid_grid2)){
	cat(i,'\n')
	grid1 = aphid_grid2[i,]
	cropland1 = crop(cropland,grid1)
	cropland2 = as.matrix(cropland1)
	pland1 = length(which(cropland2==2)) / length(cropland2)
	pland[i,1] = grid1@data$grid_id
	pland[i,2] = pland1
}
write.table(pland,'./data/croplandED_50km_MIDWEST.txt',sep='\t',row.names=F,quote=F)


#####
# UK

files = list.files('./rasters/EarthData/UK',full.names=T); files = files[grep('.tif',files)]
mosaic_rasters(files,'./rasters/UK_cropland_EarthData.tif')
cropland = raster('./rasters/UK_cropland_EarthData.tif')

aphid_sites = spTransform(readOGR(dsn='./shapefiles/aphid_sites_50km_UK.shp',layer='aphid_sites_50km_UK',verbose=F,stringsAsFactors=F),CRS(proj1))
aphid_grid = spTransform(readOGR(dsn='./shapefiles/UK_grid_50km.shp',layer='UK_grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
aphid_grid2 = aphid_grid[match(unique(aphid_sites@data$grid_id),aphid_grid@data$grid_id),]

pland = data.frame('grid_id'=aphid_grid2@data$grid_id,'Cropland'=-999)
for (i in 1:length(aphid_grid2)){
	cat(i,'\n')
	grid1 = aphid_grid2[i,]
	cropland1 = crop(cropland,grid1)
	cropland2 = as.matrix(cropland1)
	pland1 = length(which(cropland2==2)) / length(cropland2)
	pland[i,1] = grid1@data$grid_id
	pland[i,2] = pland1
}
write.table(pland,'./data/croplandED_50km_UK.txt',sep='\t',row.names=F,quote=F)

