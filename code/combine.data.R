
# Combine dataframes for aphid abundance trends, climate, land cover, and aphid traits

setwd('path/to/AphidDecline')

library(rgdal)
library(rgeos)
proj1 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'


#####
# MIDWEST

# Abundance trends
trends.dir = list.files('./inla_output/MW',full.names=T)
trends = c()
for (i in 1:length(trends.dir)){
	tr = read.table(trends.dir[i],sep='\t',as.is=T,check.names=F,header=T)
	species1 = strsplit(strsplit(trends.dir[i],'/')[[1]][11],'_')[[1]][1]
	add1 = data.frame('Species'=species1,'grid_id'=tr$grid_id,'Abundance.trend'=tr$tau)
	trends = data.frame(rbind(trends,add1),stringsAsFactors=F)
}

#Environment
climate = read.table('./data/climate_trends_CRUTS4.03_50km_MIDWEST.txt',sep='\t',as.is=T,check.names=F,header=T)
landscape = read.table('./data/cropland_change_50km_MIDWEST.txt',sep='\t',as.is=T,check.names=F,header=T)
pland = read.table('./data/croplandED_50km_MIDWEST.txt',sep='\t',as.is=T,check.names=F,header=T)
landscape = merge(landscape,pland,by='grid_id'); colnames(landscape) = c('grid_id','Cropland','Cropland.trend','PLAND')

# Traits
traits = read.table('./data/traits_MIDWEST.txt',sep='\t',as.is=T,check.names=F,header=T)
traits$summer_host = gsub(', ','_',traits$Summer_host)
traits$hetero_mono = traits$Host_alternating; traits$hetero_mono[which(traits$hetero_mono=='Heteroecious')] = 'H'; traits$hetero_mono[which(traits$hetero_mono=='Monoecious')] = 'M'; traits$hetero_mono[which(is.na(traits$hetero_mono))] = 'V'
traits$pest = traits$Pest_status; traits$pest[which(!is.na(traits$pest))] = 'Y'; traits$pest[which(is.na(traits$pest))] = 'N'
traits$primary_host = traits$Primary_host; traits$primary_host[which(is.na(traits$primary_host))] = 'None'; traits$primary_host = gsub(', ','_',traits$primary_host)
traits$holocyclic = traits$Parthenogenesis; traits$holocyclic[which(traits$holocyclic=='Holocyclic')] = 'Y'; traits$holocyclic[which(traits$holocyclic=='Anholocyclic')] = 'N'
traits$host_breadth = traits$Host_breadth
traits$host_breadth2 = traits$host_breadth; traits$host_breadth2[which(traits$host_breadth2=='O' | traits$host_breadth2=='P')] = 'OP'
traits$native = traits$Invasive; traits$native[which(traits$native=='Invasive')] = 'N'; traits$native[which(traits$native=='Native')] = 'Y'
traits$origin = traits$native; traits$origin[which(traits$origin=='Y')] = 'Nearctic'; traits$origin[which(traits$origin=='N')] = 'Palearctic'
traits$body_length = apply(cbind(traits$Minimum_length,traits$Maximum_length),1,median)
traits = traits[,c(1,18:27)]

# Spatial grid
grid50 = spTransform(readOGR('./shapefiles/grid_50km.shp','grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1))
centroids = gCentroid(grid50,byid=T)
grid50@data$Longitude = centroids@coords[,1]
grid50@data$Latitude = centroids@coords[,2]
grid50@data = grid50@data[,-1]

merge1m = merge(climate,landscape,by='grid_id')
merge2 = merge(merge1m,grid50@data,by='grid_id')
merge3 = merge(trends,traits,by='Species')
merge_MW = merge(merge3,merge2,by='grid_id') #merged data for the Midwest - use to rbind into final dataframe below
merge_MW$Dataset = 'Midwest'


#####
# IDAHO

# Abundance trends
trends.dir = list.files('./inla_output/ID',full.names=T)
trends = c()
for (i in 1:length(trends.dir)){
	tr = read.table(trends.dir[i],sep='\t',as.is=T,check.names=F,header=T)
	species1 = strsplit(strsplit(trends.dir[i],'/')[[1]][11],'_')[[1]][1]
	add1 = data.frame('Species'=species1,'grid_id'=tr$grid_id,'Abundance.trend'=tr$tau)
	trends = data.frame(rbind(trends,add1),stringsAsFactors=F)
}

# Environment
climate = read.table('./data/climate_trends_CRUTS4.03_50km_IDAHO.txt',sep='\t',as.is=T,check.names=F,header=T)
landscape = read.table('./data/cropland_change_50km_IDAHO.txt',sep='\t',as.is=T,check.names=F,header=T)
pland = read.table('./data/croplandED_50km_IDAHO.txt',sep='\t',as.is=T,check.names=F,header=T)
landscape = merge(landscape,pland,by='grid_id'); colnames(landscape) = c('grid_id','Cropland','Cropland.trend','PLAND')

# Traits
traits2 = read.table('./curated_data/traits_Idaho.txt',sep='\t',as.is=T,check.names=F,header=T)
traits2$body_length = apply(cbind(traits2$Minimum_length,traits2$Maximum_length),1,median)
traits2$host_breadth2 = traits2$host_breadth; traits2$host_breadth2[which(traits2$host_breadth2=='O' | traits2$host_breadth2=='P')] = 'OP'
traits2 = traits2[,c(1:9,12:13)]

merge1i = merge(climate,landscape,by='grid_id')
merge2 = merge(merge1i,grid50@data,by='grid_id',all=F)
merge3 = merge(trends,traits2,by='Species')
merge_ID = merge(merge3,merge2,by='grid_id') #merged data for the Midwest - use to rbind into final dataframe below
merge_ID$Dataset = 'Idaho'

# Update species synonyms
merge_ID$Species[which(merge_ID$Species=='Aphis helianthi')] = 'Aphis asclepiadis'
merge_ID$Species[which(merge_ID$Species=='Aphis armoraciae')] = 'Protaphis middletonii'
merge_ID$Species[which(merge_ID$Species=='Amphorophora crataegi')] = 'Utamphorophora crataegi'
merge_ID$Species[which(merge_ID$Species=='Aphis avicularis')] = 'Aphis polygonata'
merge_ID$Species[which(merge_ID$Species=='Pleotrichophorus rusticatus')] = 'Pleotrichophorus pullus'
merge_ID$Species[which(merge_ID$Species=='Rhopalosiphum insertum')] = 'Rhopalosiphum oxyacanthae'
merge_ID$Species[which(merge_ID$Species=='Hyadaphis tartaricae')] = 'Hyadaphis tataricae'

US_merge = rbind(merge_MW,merge_ID)
merge1mi = rbind(merge1m,merge1i)
shp.mi = merge(grid50,merge1mi,by='grid_id',all.x=F)


#####
# UK

# Abundance trends
trends.dir = list.files('./inla_output/UK',full.names=T)
trends = c()
for (i in 1:length(trends.dir)){
	tr = read.table(trends.dir[i],sep='\t',as.is=T,check.names=F,header=T)
	species1 = strsplit(strsplit(trends.dir[i],'/')[[1]][10],'_')[[1]][1]
	add1 = data.frame('Species'=species1,'grid_id'=tr$grid_id,'Abundance.trend'=tr$tau)
	trends = data.frame(rbind(trends,add1),stringsAsFactors=F)
}

# Environment
climate = read.table('./data/climate_trends_CRUTS4.03_50km_UK.txt',sep='\t',as.is=T,check.names=F,header=T)
pland = read.table('./data/croplandED_50km_UK.txt',sep='\t',as.is=T,check.names=F,header=T)
landscape = data.frame('grid_id'=pland[,1],'Cropland'=NA,'Cropland.trend'=NA,'PLAND'=pland[,2])

# Traits
traits3 = read.table('./curated_data/traits_UK.txt',sep='\t',as.is=T,check.names=F,header=T)
traits3$body_length = apply(cbind(traits3$Minimum_length,traits3$Maximum_length),1,median)
traits3$summer_host = gsub(', ','_',traits3$summer_host)
traits3$pest[which(traits3$pest!='N')] = 'Y'
traits3$primary_host = gsub(', ','_',traits3$primary_host); traits3$primary_host[which(is.na(traits3$primary_host))] = 'None'
traits3$native = NA #invasive status annotations are doubtful
traits3$host_breadth2 = traits3$host_breadth; traits3$host_breadth2[which(traits3$host_breadth2=='O' | traits3$host_breadth2=='P')] = 'OP'
traits3 = traits3[,c(1:2,4,3,5:9,12:13)]

# Spatial grid
grid50 = spTransform(readOGR('./shapefiles/UK_grid_50km.shp','UK_grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1))
centroids = gCentroid(grid50,byid=T)
grid50@data$Longitude = centroids@coords[,1]
grid50@data$Latitude = centroids@coords[,2]
grid50@data = grid50@data[,-c(1:5)]

merge1u = merge(climate,landscape,by='grid_id')
merge2 = merge(merge1u,grid50@data,by='grid_id',all=F)
merge3 = merge(trends,traits3,by='Species')
merge_UK = merge(merge3,merge2,by='grid_id') #merged data for the Midwest - use to rbind into final dataframe below
merge_UK$Dataset = 'UK'; str(merge_UK)

shp.u = merge(grid50,merge1u,by='grid_id',all.x=F)


# Merge UK, MIDWEST, and IDAHO datasets
all_merge = rbind(merge_MW,merge_ID,merge_UK)
write.table(all_merge,'./data/aphids_trends_envars_traits_50km_m5_MW_ID_UK.txt',sep='\t',quote=F,row.names=F)

shp = rbind(shp.mi,spTransform(shp.u,CRS(proj4string(shp.mi))))
writeOGR(shp,'./shapefiles/envars_midwest_idaho_uk.shp','envars_midwest_idaho_uk',driver='ESRI Shapefile',overwrite=T)

