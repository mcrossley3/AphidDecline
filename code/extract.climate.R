
###########################################################################
# Annual means from CRU TS 4.03, per 50x50km grid cell, for aphid apocalypse ms

# http://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.03/data/ #source of data - pre and tmp
# https://doi.org/10.1002/joc.3711

setwd('path/to/AphidDecline')

library(rgdal)
library(raster)
library(rgeos)
source('./code/AR_reml.R')

proj1 = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' #proj4string for CRU_TS raster


################
# MIDWEST

aphid_sites = spTransform(readOGR(dsn='./shapefiles/aphid_sites_50km.shp',layer='aphid_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1))
aphid_grid = spTransform(readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
aphid_grid2 = aphid_grid[match(unique(aphid_sites@data$grid_id),aphid_grid@data$grid_id),]

# Temperature
brick.tmp = brick('path/to/cru_ts4.03.1901.2018.tmp.dat.nc') #netcdf info http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm
extract.tmp = extract(brick.tmp,aphid_grid2,fun=mean,df=T,na.rm=T)
tmp.years = apply(array(colnames(extract.tmp)),1,function(x){sub('X','',strsplit(x,'[.]')[[1]][1])})
u.tmp.years = unique(tmp.years[which(tmp.years!='ID')])
tmp.yearly.means = data.frame('grid_id'=NA,'Year'=-999,'Mean.tmp'=-999)
tx = 1
for (g in 1:length(aphid_grid2@data$grid_id)){
	noquote(print(paste0(g,' of ',length(aphid_grid2@data$grid_id))))
	for (y in 1:length(u.tmp.years)){
		tmp.yearly.means[tx,1] = aphid_grid2@data$grid_id[g]
		tmp.yearly.means[tx,2] = u.tmp.years[y]
		tmp.yearly.means[tx,3] = mean(unlist(extract.tmp[g,which(tmp.years==u.tmp.years[y])]))
		tx = tx + 1
	}
}

# Precipitation
brick.pre = brick('path/to/cru_ts4.03.1901.2018.pre.dat.nc')
extract.pre = extract(brick.pre,aphid_grid2,fun=mean,df=T,na.rm=T)
pre.years = apply(array(colnames(extract.pre)),1,function(x){sub('X','',strsplit(x,'[.]')[[1]][1])})
u.pre.years = unique(pre.years[which(pre.years!='ID')])
pre.yearly.means = data.frame('grid_id'=NA,'Year'=-999,'Mean.pre'=-999)
tx = 1
for (g in 1:length(aphid_grid2@data$grid_id)){
	noquote(print(paste0(g,' of ',length(aphid_grid2@data$grid_id))))
	for (y in 1:length(u.pre.years)){
		pre.yearly.means[tx,1] = aphid_grid2@data$grid_id[g]
		pre.yearly.means[tx,2] = u.pre.years[y]
		pre.yearly.means[tx,3] = mean(unlist(extract.pre[g,which(pre.years==u.pre.years[y])]))
		tx = tx + 1
	}
}

# Check for gids with missing data
na.count = apply(array(aphid_grid2@data$grid_id),1,function(x){length(which(is.na(tmp.yearly.means$Mean.tmp[which(tmp.yearly.means$grid_id==x)])))})
na.gids = aphid_grid2@data$grid_id[which(na.count>0)]; na.gids
na.count = apply(array(aphid_grid2@data$grid_id),1,function(x){length(which(is.na(pre.yearly.means$Mean.pre[which(pre.yearly.means$grid_id==x)])))})
na.gids = aphid_grid2@data$grid_id[which(na.count>0)]; na.gids

# Dataframe
annualmeans = data.frame('grid_id'=aphid_grid2@data$grid_id,'tmp.modern'=-999,'tmp.trend'=-999,'pre.modern'=-999,'pre.trend'=-999)
for (i in 1:nrow(annualmeans)){
	print(i)
	# Estimate trend in precip
	pre1 = pre.yearly.means[which(pre.yearly.means$grid_id==annualmeans$grid_id[i]),]
	mean.pre.present = mean(pre1$Mean.pre[match(2005:2018,pre1$Year)])
	t.scale = 1901:2018 - 1901
	X1 = pre1$Mean.pre
	Z1 = (X1 – mean(X1, na.rm=T))/sd(X1, na.rm=T) # z-transform
	arr.Z = AR_reml(Z1 ~ t.scale) #Z-transformed time trends
	pre.trend = arr.Z$coef[2] #slope of trend line
	png(paste0('./plots/climate_change/MIDWEST_pre_',annualmeans$grid_id[i],'.png'))
	plot(1901:2018,X1,type='l',lwd=3,xlab='Year',ylab='Avg. precip. (mm)',cex.lab=1.5,cex.axis=1.2)
	dev.off()
	# Estimate trend in temp
	tmp1 = tmp.yearly.means[which(tmp.yearly.means$grid_id==annualmeans$grid_id[i]),]
	mean.tmp.present = mean(tmp1$Mean.tmp[match(2005:2018,tmp1$Year)])
	X2 = tmp1$Mean.tmp
	Z2 = (X2 – mean(X2, na.rm=T))/sd(X2, na.rm=T) # z-transform
	arr.Z = AR_reml(Z2 ~ t.scale) #Z-transformed time trends
	tmp.trend = arr.Z$coef[2] #slope of trend line
	png(paste0('./plots/climate_change/MIDWEST_tmp_',annualmeans$grid_id[i],'.png'))
	plot(1901:2018,X2,type='l',lwd=3,xlab='Year',ylab='Avg. temp. (C)',cex.lab=1.5,cex.axis=1.2)
	dev.off()
	# dataframe
	annualmeans[i,2] = mean.tmp.present
	annualmeans[i,3] = tmp.trend
	annualmeans[i,4] = mean.pre.present
	annualmeans[i,5] = pre.trend
}
write.table(annualmeans,'./data/climate_trends_CRUTS4.03_50km_MIDWEST.txt',sep='\t',quote=F,row.names=F)


################
# IDAHO

aphid_sites = spTransform(readOGR(dsn='./shapefiles/PNW_aphid_sites_50km.shp',layer='PNW_aphid_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1))
aphid_grid = spTransform(readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
aphid_grid2 = aphid_grid[match(unique(aphid_sites@data$grid_id),aphid_grid@data$grid_id),]

# Temperature
brick.tmp = brick('path/to/cru_ts4.03.1901.2018.tmp.dat.nc') #netcdf info http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm
extract.tmp = extract(brick.tmp,aphid_grid2,fun=mean,df=T,na.rm=T)
tmp.years = apply(array(colnames(extract.tmp)),1,function(x){sub('X','',strsplit(x,'[.]')[[1]][1])})
u.tmp.years = unique(tmp.years[which(tmp.years!='ID')])
tmp.yearly.means = data.frame('grid_id'=NA,'Year'=-999,'Mean.tmp'=-999)
tx = 1
for (g in 1:length(aphid_grid2@data$grid_id)){
	noquote(print(paste0(g,' of ',length(aphid_grid2@data$grid_id))))
	for (y in 1:length(u.tmp.years)){
		tmp.yearly.means[tx,1] = aphid_grid2@data$grid_id[g]
		tmp.yearly.means[tx,2] = u.tmp.years[y]
		tmp.yearly.means[tx,3] = mean(unlist(extract.tmp[g,which(tmp.years==u.tmp.years[y])]))
		tx = tx + 1
	}
}

# Precipitation
brick.pre = brick('path/to/cru_ts4.03.1901.2018.pre.dat.nc')
extract.pre = extract(brick.pre,aphid_grid2,fun=mean,df=T,na.rm=T)
pre.years = apply(array(colnames(extract.pre)),1,function(x){sub('X','',strsplit(x,'[.]')[[1]][1])})
u.pre.years = unique(pre.years[which(pre.years!='ID')])
pre.yearly.means = data.frame('grid_id'=NA,'Year'=-999,'Mean.pre'=-999)
tx = 1
for (g in 1:length(aphid_grid2@data$grid_id)){
	noquote(print(paste0(g,' of ',length(aphid_grid2@data$grid_id))))
	for (y in 1:length(u.pre.years)){
		pre.yearly.means[tx,1] = aphid_grid2@data$grid_id[g]
		pre.yearly.means[tx,2] = u.pre.years[y]
		pre.yearly.means[tx,3] = mean(unlist(extract.pre[g,which(pre.years==u.pre.years[y])]))
		tx = tx + 1
	}
}

# Check for gids with missing data
na.count = apply(array(aphid_grid2@data$grid_id),1,function(x){length(which(is.na(tmp.yearly.means$Mean.tmp[which(tmp.yearly.means$grid_id==x)])))})
na.gids = aphid_grid2@data$grid_id[which(na.count>0)]; na.gids
na.count = apply(array(aphid_grid2@data$grid_id),1,function(x){length(which(is.na(pre.yearly.means$Mean.pre[which(pre.yearly.means$grid_id==x)])))})
na.gids = aphid_grid2@data$grid_id[which(na.count>0)]; na.gids

# Dataframe
annualmeans = data.frame('grid_id'=aphid_grid2@data$grid_id,'tmp.modern'=-999,'tmp.trend'=-999,'pre.modern'=-999,'pre.trend'=-999)
for (i in 1:nrow(annualmeans)){
	print(i)
	# Estimate trend in precip
	pre1 = pre.yearly.means[which(pre.yearly.means$grid_id==annualmeans$grid_id[i]),]
	mean.pre.present = mean(pre1$Mean.pre[match(1985:2001,pre1$Year)])
	t.scale = 1901:2001 - 1901
	X1 = pre1$Mean.pre[match(1901:2001,u.pre.years)]
	Z1 = (X1 – mean(X1, na.rm=T))/sd(X1, na.rm=T) # z-transform
	arr.Z = AR_reml(Z1 ~ t.scale) #Z-transformed time trends
	pre.trend = arr.Z$coef[2] #slope of trend line
	png(paste0('./plots/climate_change/IDAHO_pre_',annualmeans$grid_id[i],'.png'))
	plot(1901:2001,X1,type='l',lwd=3,xlab='Year',ylab='Avg. precip. (mm)',cex.lab=1.5,cex.axis=1.2)
	dev.off()
	# Estimate trend in temp
	tmp1 = tmp.yearly.means[which(tmp.yearly.means$grid_id==annualmeans$grid_id[i]),]
	mean.tmp.present = mean(tmp1$Mean.tmp[match(1985:2001,tmp1$Year)])
	X2 = tmp1$Mean.tmp[match(1901:2001,u.tmp.years)]
	Z2 = (X2 – mean(X2, na.rm=T))/sd(X2, na.rm=T) # z-transform
	arr.Z = AR_reml(Z2 ~ t.scale) #Z-transformed time trends
	tmp.trend = arr.Z$coef[2] #slope of trend line
	png(paste0('./plots/climate_change/IDAHO_tmp_',annualmeans$grid_id[i],'.png'))
	plot(1901:2001,X2,type='l',lwd=3,xlab='Year',ylab='Avg. temp. (C)',cex.lab=1.5,cex.axis=1.2)
	dev.off()
	# dataframe
	annualmeans[i,2] = mean.tmp.present
	annualmeans[i,3] = tmp.trend
	annualmeans[i,4] = mean.pre.present
	annualmeans[i,5] = pre.trend
}
write.table(annualmeans,'./data/climate_trends_CRUTS4.03_50km_IDAHO.txt',sep='\t',quote=F,row.names=F)


################
# UK

aphid_sites = spTransform(readOGR(dsn='./shapefiles/aphid_sites_50km_UK.shp',layer='aphid_sites_50km_UK',verbose=F,stringsAsFactors=F),CRS(proj1))
aphid_grid = spTransform(readOGR(dsn='./shapefiles/UK_grid_50km.shp',layer='UK_grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
aphid_grid2 = aphid_grid[match(unique(aphid_sites@data$grid_id),aphid_grid@data$grid_id),]

# Temperature
brick.tmp = brick('path/to/cru_ts4.03.1901.2018.tmp.dat.nc') #netcdf info http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm
extract.tmp = extract(brick.tmp,aphid_grid2,fun=mean,df=T,na.rm=T)
tmp.years = apply(array(colnames(extract.tmp)),1,function(x){sub('X','',strsplit(x,'[.]')[[1]][1])})
u.tmp.years = unique(tmp.years[which(tmp.years!='ID')])
tmp.yearly.means = data.frame('grid_id'=NA,'Year'=-999,'Mean.tmp'=-999)
tx = 1
for (g in 1:length(aphid_grid2@data$grid_id)){
	noquote(print(paste0(g,' of ',length(aphid_grid2@data$grid_id))))
	for (y in 1:length(u.tmp.years)){
		tmp.yearly.means[tx,1] = aphid_grid2@data$grid_id[g]
		tmp.yearly.means[tx,2] = u.tmp.years[y]
		tmp.yearly.means[tx,3] = mean(unlist(extract.tmp[g,which(tmp.years==u.tmp.years[y])]))
		tx = tx + 1
	}
}

# Precipitation
brick.pre = brick('path/to/cru_ts4.03.1901.2018.pre.dat.nc')
extract.pre = extract(brick.pre,aphid_grid2,fun=mean,df=T,na.rm=T)
pre.years = apply(array(colnames(extract.pre)),1,function(x){sub('X','',strsplit(x,'[.]')[[1]][1])})
u.pre.years = unique(pre.years[which(pre.years!='ID')])
pre.yearly.means = data.frame('grid_id'=NA,'Year'=-999,'Mean.pre'=-999)
tx = 1
for (g in 1:length(aphid_grid2@data$grid_id)){
	noquote(print(paste0(g,' of ',length(aphid_grid2@data$grid_id))))
	for (y in 1:length(u.pre.years)){
		pre.yearly.means[tx,1] = aphid_grid2@data$grid_id[g]
		pre.yearly.means[tx,2] = u.pre.years[y]
		pre.yearly.means[tx,3] = mean(unlist(extract.pre[g,which(pre.years==u.pre.years[y])]))
		tx = tx + 1
	}
}

# Check for gids with missing data
na.count = apply(array(aphid_grid2@data$grid_id),1,function(x){length(which(is.na(tmp.yearly.means$Mean.tmp[which(tmp.yearly.means$grid_id==x)])))})
na.gids = aphid_grid2@data$grid_id[which(na.count>0)]; na.gids
na.count = apply(array(aphid_grid2@data$grid_id),1,function(x){length(which(is.na(pre.yearly.means$Mean.pre[which(pre.yearly.means$grid_id==x)])))})
na.gids = aphid_grid2@data$grid_id[which(na.count>0)]; na.gids

# Dataframe
annualmeans = data.frame('grid_id'=aphid_grid2@data$grid_id,'tmp.modern'=-999,'tmp.trend'=-999,'pre.modern'=-999,'pre.trend'=-999)
for (i in 1:nrow(annualmeans)){
	print(i)
	# Estimate trend in precip
	pre1 = pre.yearly.means[which(pre.yearly.means$grid_id==annualmeans$grid_id[i]),]
	mean.pre.present = mean(pre1$Mean.pre[match(1965:2018,pre1$Year)])
	t.scale = 1901:2018 - 1901
	X1 = pre1$Mean.pre
	Z1 = (X1 – mean(X1, na.rm=T))/sd(X1, na.rm=T) # z-transform
	arr.Z = AR_reml(Z1 ~ t.scale) #Z-transformed time trends
	pre.trend = arr.Z$coef[2] #slope of trend line
	png(paste0('./plots/climate_change/UK_pre_',annualmeans$grid_id[i],'.png'))
	plot(1901:2018,X1,type='l',lwd=3,xlab='Year',ylab='Avg. precip. (mm)',cex.lab=1.5,cex.axis=1.2)
	dev.off()
	# Estimate trend in temp
	tmp1 = tmp.yearly.means[which(tmp.yearly.means$grid_id==annualmeans$grid_id[i]),]
	mean.tmp.present = mean(tmp1$Mean.tmp[match(1965:2018,tmp1$Year)])
	X2 = tmp1$Mean.tmp
	Z2 = (X2 – mean(X2, na.rm=T))/sd(X2, na.rm=T) # z-transform
	arr.Z = AR_reml(Z2 ~ t.scale) #Z-transformed time trends
	tmp.trend = arr.Z$coef[2] #slope of trend line
	png(paste0('./plots/climate_change/UK_tmp_',annualmeans$grid_id[i],'.png'))
	plot(1901:2018,X2,type='l',lwd=3,xlab='Year',ylab='Avg. temp. (C)',cex.lab=1.5,cex.axis=1.2)
	dev.off()
	# dataframe
	annualmeans[i,2] = mean.tmp.present
	annualmeans[i,3] = tmp.trend
	annualmeans[i,4] = mean.pre.present
	annualmeans[i,5] = pre.trend
}
write.table(annualmeans,'./data/climate_trends_CRUTS4.03_50km_UK.txt',sep='\t',quote=F,row.names=F)

