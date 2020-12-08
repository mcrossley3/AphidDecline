
setwd('path/to/AphidDecline')

library(lubridate)
library(rgdal)
library(spdep)

convert.date = function(mdy){ #function that converts dates into a format readable by as.Date function
	mdy.split = strsplit(mdy,'/')[[1]]
	YYYY = mdy.split[3]
	MM = apply(array(mdy.split[1]),1,function(x){if(nchar(x)<2){paste0(0,x)}else{x}})
	DD = apply(array(mdy.split[2]),1,function(x){if(nchar(x)<2){paste0(0,x)}else{x}})
	return(paste(YYYY,MM,DD,sep='-'))
}

##################################################
# MIDWEST

# Import aphid data
STN = read.table('./data/STN_counts_curated_time-consistent2.txt',sep='\t',as.is=T,check.names=F,header=T) #data available upon reasonable request
STN = STN[which(STN$Site!='Eureka'),]
STN$Date = as.Date(apply(array(STN$Date),1,function(z){convert.date(z)}))

# Aggregate counts by species*site*year, standardizing by effort (trap-days)
STN.annual = data.frame('Site'=NA,'Longitude'=-999,'Latitude'=-999,'Year'=-999,'Species'=NA,'Trap.days'=-999,'N.aphids'=-999,'Aphids.day'=-999,'Julian.start'=NA,'Julian.end'=NA)
sites = sort(unique(STN$Site))
species = sort(unique(STN$Species))
ix = 1
for (s in 1:length(sites)){
	noquote(print(sites[s]))
	dat1 = STN[which(STN$Site==sites[s]),]
	years = sort(unique(dat1$Year))
	lat = unique(dat1$Latitude)
	lon = unique(dat1$Longitude)
	for (y in 1:length(years)){
		dat2 = dat1[which(dat1$Year==years[y]),]
		# How many weeks were sampled at this site*year?
		dates = unique(dat2$Date)
		wk.count = apply(array(dates),1,function(x){length(which(is.na(dat2$Female.count[which(dat2$Date==x)])))})
		n.wk = length(which(wk.count==0))
		n.day = n.wk * 7
		for (i in 1:length(species)){
			pos = which(dat2$Species==species[i])
			if (length(pos)==0){ #species absent = true zero
				n.aphids = n.aphid.day = 0
			} else {
				dat3 = dat2[pos,]
				n.aphids = sum(dat3$Female.count,na.rm=T)
				n.aphid.day = n.aphids / n.day
			}
			STN.annual[ix,1] = sites[s] #Site
			STN.annual[ix,2] = lon #Longitude
			STN.annual[ix,3] = lat #Latitude
			STN.annual[ix,4] = years[y] #Year
			STN.annual[ix,5] = species[i] #Species
			STN.annual[ix,6] = n.day #Trap.days
			STN.annual[ix,7] = n.aphids #N.aphids
			STN.annual[ix,8] = n.aphid.day #Aphids.day
			STN.annual[ix,9] = yday(sort(dates)[1]) #First.trap.date
			STN.annual[ix,10] = yday(sort(dates,decreasing=T)[1]) #Last.trap.date
			ix = ix + 1
		}
	}
}
write.table(STN.annual,'./curated_data/aphid_counts_annual_MIDWEST.txt',sep='\t',quote=F,row.names=F)

check = STN.annual[which(STN.annual$Species=='Aphis glycines' & STN.annual$Site=='Arlington'),]

#####
# Filter out species*site that had < 5 non-zero records

STN.annual = read.table('./data/aphid_counts_annual_MIDWEST.txt',sep='\t',as.is=T,check.names=F,header=T)
STN.annual$SpeciesSite = paste(STN.annual$Species,STN.annual$Site,sep='_')
ss = unique(STN.annual$SpeciesSite)
ss.nonzero = apply(array(ss),1,function(x){length(which(STN.annual$N.aphids[which(STN.annual$SpeciesSite==x)]>0))})
ss.keep = ss[which(ss.nonzero>4)]
keep.pos = unlist(apply(array(ss.keep),1,function(x){which(STN.annual$SpeciesSite==x)}))
STN.annual2 = STN.annual[keep.pos,]; dim(STN.annual); dim(STN.annual2)
write.table(STN.annual2,'./curated_data/aphid_counts_annual_MIDWEST_m5.txt',sep='\t',quote=F,row.names=F)


###########################################################
# Idaho

convert.date2 = function(mdy){ #function that converts dates into a format readable by as.Date function
	YYYY = mdy$YEAR
	MM = mdy$MONTH
	DD = mdy$DAY
	return(paste(YYYY,MM,DD,sep='-'))
}

count.weeks = function(dds){
	jumps = c()
	for (d1 in 2:length(dds)){
		jumps = c(jumps,dds[d1] - dds[d1-1])
	}
}

# Import aphid data
STN = read.csv('./data/PNW_suctiontrap_counts.csv',as.is=T,check.names=F,header=T) # data available upon reasonable request
STN = STN[which(!is.na(STN$Species)),]
STN = STN[which(!is.na(STN$N)),]
species = sort(unique(STN$Species))
sites = sort(unique(STN$Site))

# Match coordinates to sites
coords = read.csv('./data/PNW_coordinates.csv',as.is=T,check.names=F,header=T) # data available upon reasonable request
STN$Lon = STN$Lat = NA
for (i in 1:length(coords$Site)){
	STN$Lon[which(STN$Site==coords$Site[i])] = coords$Longitude[i]
	STN$Lat[which(STN$Site==coords$Site[i])] = coords$Latitude[i]
}
STN$Year = STN$YEAR
STN$Aphid.count = STN$N
sum(STN$Aphid.count,na.rm=T) #658,048 aphids identified
STN$Date = as.Date(convert.date2(STN[,3:5]))
STN$Julian = yday(STN$Date)

# Aggregate counts by species*site*year, standardizing by effort (trap-days)
STN.annual = data.frame('Site'=NA,'Longitude'=-999,'Latitude'=-999,'Year'=-999,'Species'=NA,'Trap.days'=-999,'N.aphids'=-999,'Aphids.day'=-999,'Julian.start'=NA,'Julian.end'=NA)
sites = sort(unique(STN$Site))
species = sort(unique(STN$Species))
ix = 1
for (s in 1:length(sites)){
	noquote(print(sites[s]))
	dat1 = STN[which(STN$Site==sites[s]),]
	years = sort(unique(dat1$Year))
	lat = unique(dat1$Lat)
	lon = unique(dat1$Lon)
	for (y in 1:length(years)){
		dat2 = dat1[which(dat1$Year==years[y]),]
		# How many weeks were sampled at this site*year?
		dates = unique(dat2$Julian)		
		wk.count = apply(array(dates),1,function(x){length(which(is.na(dat2$Aphid.count[which(dat2$Julian==x)])))})
		n.wk = length(which(wk.count==0))
		n.day = n.wk * 7
		for (i in 1:length(species)){
			pos = which(dat2$Species==species[i])
			if (length(pos)==0){ #species absent = true zero
				n.aphids = n.aphid.day = 0
			} else {
				dat3 = dat2[pos,]
				n.aphids = sum(dat3$Aphid.count,na.rm=T)
				n.aphid.day = n.aphids / n.day
			}
			STN.annual[ix,1] = sites[s] #Site
			STN.annual[ix,2] = lon #Longitude
			STN.annual[ix,3] = lat #Latitude
			STN.annual[ix,4] = years[y] #Year
			STN.annual[ix,5] = species[i] #Species
			STN.annual[ix,6] = n.day #Trap.days
			STN.annual[ix,7] = n.aphids #N.aphids
			STN.annual[ix,8] = n.aphid.day #Aphids.day
			STN.annual[ix,9] = sort(dates)[1] #First.trap.date
			STN.annual[ix,10] = sort(dates,decreasing=T)[1] #Last.trap.date
			ix = ix + 1
		}
	}
}
write.table(STN.annual,'./curated_data/aphid_counts_annual_NORTHWEST.txt',sep='\t',quote=F,row.names=F)

check = STN.annual[which(STN.annual$Species=='Diuraphis noxia' & STN.annual$Site=='Moscow'),]

#####
# Filter out species*site that had < 5 non-zero records

STN.annual = read.table('./curated_data/aphid_counts_annual_NORTHWEST.txt',sep='\t',as.is=T,check.names=F,header=T)
STN.annual$SpeciesSite = paste(STN.annual$Species,STN.annual$Site,sep='_')
ss = unique(STN.annual$SpeciesSite)
ss.nonzero = apply(array(ss),1,function(x){length(which(STN.annual$N.aphids[which(STN.annual$SpeciesSite==x)]>0))})
ss.keep = ss[which(ss.nonzero>4)]
keep.pos = unlist(apply(array(ss.keep),1,function(x){which(STN.annual$SpeciesSite==x)}))
STN.annual2 = STN.annual[keep.pos,]; dim(STN.annual); dim(STN.annual2)
write.table(STN.annual2,'./curated_data/aphid_counts_annual_NORTHWEST_m5.txt',sep='\t',quote=F,row.names=F)


##################################################
# UK

# Import aphid data
rfiles = list.files('./data/Rothamsted/data_files',full.names=T) # data available upon reasonable request
rfiles = rfiles[-grep('Inoperative',rfiles)]
UK = read.csv(rfiles[1],as.is=T,check.names=F,header=T)
for (i in 2:length(rfiles)){
	add = read.csv(rfiles[i],as.is=T,check.names=F,header=T)
	UK = data.frame(rbind(UK,add),stringsAsFactors=F,check.names=F)
}
UK$Site = UK$TrapName
UK$Species = UK$binomial

# Rename sub-species
key = read.table('./data/Rothamsted_aphid_species_key.txt',sep='\t',as.is=T,check.names=F,header=F)
species = sort(unique(UK$Species))
for (i in 1:length(species)){
	pos = which(key[,1]==species[i])
	if (length(pos)==0){
		print(species[i])
	} else {
		nn = key[pos,2]
		UK$Species[which(UK$Species==species[i])] = nn
	}
}

#Remove genus-level and non-aphid records
nrow(UK)
adelgids = grep('ADELGIDAE',UK$Species); length(adelgids)
spp = grep(' spp',UK$Species); length(spp)
unknown = which(UK$Species=='Unknown' | UK$Species=='Others'); length(unknown)
UK = UK[-c(adelgids,spp,unknown),]; nrow(UK)
sort(unique(UK$Species)) #364 species

# Import dates when traps were inoperative
rfiles = list.files('./data/Rothamsted/data_files',full.names=T)
rfiles = rfiles[grep('Inoperative',rfiles)]
UK.inop = read.csv(rfiles[1],as.is=T,check.names=F,header=T)
for (i in 2:length(rfiles)){
	add = read.csv(rfiles[i],as.is=T,check.names=F,header=T)
	UK.inop = data.frame(rbind(UK.inop,add),stringsAsFactors=F,check.names=F)
}

traps = c('Rothamsted','Wye','Brooms Barn','Newcastle','Silwood','Elgin','Hereford','Preston','Ayr','Writtle','Long Ashton II','Elgin II','Dundee','East Craigs')
lats = c(51.806997,51.185507,52.260681,55.213254,51.409410,57.645285,52.124201,53.854383,55.476826,51.733599,51.424815,57.650850,56.457147,55.949386)
lons = c(-0.360091,0.944941,0.568430,-1.685083,-0.643357,-3.365517,-2.638156,-2.766990,-4.567922,0.429233,-2.667435,-3.425272,-3.073650,-3.311634)
coords = data.frame('Site'=traps,'Lon'=lons,'Lat'=lats)

# Align trap names with coords key
UK$Site[which(UK$Site=="Brooms' Barn S")] = 'Brooms Barn'
UK$Site[which(UK$Site=='Elgin One')] = 'Elgin'
e2 = UK$Site[grep('Elgin',UK$Site)]
UK$Site[which(UK$Site== e2[19375])] = 'Elgin II'
UK$Site[which(UK$Site=='Hereford (Rosemaund)')] = 'Hereford'
a2 = UK$Site[grep('Ashton',UK$Site)]
UK$Site[which(UK$Site==a2[1])] = 'Long Ashton II'
UK$Site[which(UK$Site=='Writtle S')] = 'Writtle'
sort(unique(UK$Site))

# Aggregate counts by species*site*year, standardizing by effort (trap-days)
UK.annual = data.frame('Site'=NA,'Longitude'=-999,'Latitude'=-999,'Year'=-999,'Species'=NA,'Trap.days'=-999,'N.aphids'=-999,'Aphids.day'=-999,'Julian.start'=NA,'Julian.end'=NA)
sites = sort(unique(UK$Site))
species = sort(unique(UK$Species))
leap.years = c(1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016)
days.in.year = cbind(1965:2018,rep(365,length(1965:2018)))
days.in.year[match(leap.years,days.in.year[,1]),2] = 366
ix = 1
for (s in 1:length(sites)){
	noquote(print(sites[s]))
	dat1 = UK[which(UK$Site==sites[s]),]
	years = sort(unique(dat1$CalYear))
	lat = coords$Lat[which(coords$Site==sites[s])]
	lon = coords$Lon[which(coords$Site==sites[s])]
	for (y in 1:length(years)){
		noquote(print(paste0('   ',years[y])))
		dat2 = dat1[which(dat1$CalYear==years[y]),]
		# How many weeks were sampled at this site*year?
		inop1 = UK.inop[which(UK.inop$Trap==sites[s] & UK.inop$Year==years[y]),] #record of inoperative days for site*year
		dy = days.in.year[which(days.in.year[,1]==years[y]),2] #how many days in year y
		opm = match(1:dy,inop1$CalDay)
		op.days = which(is.na(opm))
		n.day = length(op.days)
		for (i in 1:length(species)){
			pos = which(dat2$Species==species[i])
			if (length(pos)==0){ #species absent = true zero
				n.aphids = n.aphid.day = 0
			} else {
				dat3 = dat2[pos,]
				n.aphids = sum(dat3$WeeklyCount,na.rm=T)
				n.aphid.day = n.aphids / n.day
				UK.annual[ix,1] = sites[s] #Site
				UK.annual[ix,2] = lon #Longitude
				UK.annual[ix,3] = lat #Latitude
				UK.annual[ix,4] = years[y] #Year
				UK.annual[ix,5] = species[i] #Species
				UK.annual[ix,6] = n.day #Trap.days
				UK.annual[ix,7] = n.aphids #N.aphids
				UK.annual[ix,8] = n.aphid.day #Aphids.day
				UK.annual[ix,9] = min(op.days) #First.trap.date
				UK.annual[ix,10] = max(op.days) #Last.trap.date
				ix = ix + 1
			}
		}
	}
}
write.table(UK.annual,'./data/aphid_counts_annual_UK.txt',sep='\t',quote=F,row.names=F)

check = UK.annual[which(UK.annual$Species=='Rhopalosiphum padi' & UK.annual$Site=='Wye'),]

#####
# Filter out species*site that had < 5 or 10 non-zero records

STN.annual = read.table('./data/aphid_counts_annual_UK.txt',sep='\t',as.is=T,check.names=F,header=T)
STN.annual$SpeciesSite = paste(STN.annual$Species,STN.annual$Site,sep='_')
ss = unique(STN.annual$SpeciesSite)
ss.nonzero = apply(array(ss),1,function(x){length(which(STN.annual$N.aphids[which(STN.annual$SpeciesSite==x)]>0))})
ss.keep = ss[which(ss.nonzero>4)]
keep.pos = unlist(apply(array(ss.keep),1,function(x){which(STN.annual$SpeciesSite==x)}))
STN.annual2 = STN.annual[keep.pos,]; dim(STN.annual); dim(STN.annual2)
write.table(STN.annual2,'./curated_data/aphid_counts_annual_UK_m5.txt',sep='\t',quote=F,row.names=F)

