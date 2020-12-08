
###############################################################
# MIDWEST: Prepare data for INLA

setwd('path/to/AphidDecline')

library(rgdal)
library(spdep)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs84 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

aphids = read.table('./data/aphid_counts_annual_MIDWEST_m5.txt',sep='\t',as.is=T,check.names=F,header=T) # data available upon reasonable request

sites = unique(aphids$Site)
out = data.frame('Site'=sites,'Lon'=aphids$Longitude[match(sites,aphids$Site)],'Lat'=aphids$Latitude[match(sites,aphids$Site)],stringsAsFactors=F)
shp = spTransform(SpatialPointsDataFrame(coords=cbind(as.numeric(out$Lon),as.numeric(out$Lat)),data=out,proj4string=CRS(wgs84)),CRS(proj1))

# Create 50km gridded dataset
grid50 = readOGR(dsn='./shapefiles/50km_grid_clip.shp',layer='50km_grid_clip',verbose=F,stringsAsFactors=F)
grid50@data$grid_id = seq(1,length(grid50),1)
writeOGR(grid50,dsn='./shapefiles/grid_50km.shp',layer='grid_50km',driver='ESRI Shapefile',overwrite=T)
bover = over(shp,grid50)
shp@data$grid_id = bover$grid_id
writeOGR(shp,dsn='./shapefiles/aphid_sites_m5_50km.shp',layer='aphid_sites_m5_50km',driver='ESRI Shapefile',overwrite=T)
gids = sort(unique(shp@data$grid_id))
aphids$grid_id = NA
for (g in 1:length(gids)){
	site1 = shp@data$Site[which(shp@data$grid_id==gids[g])]
	for (i in 1:length(site1)){
		row.pos = which(aphids$Site==site1[i])
		aphids$grid_id[row.pos] = gids[g]	
	}
}
write.table(aphids,'./data/aphid_counts_annual_MIDWEST_m5_50km.txt',sep='\t',row.names=F)
nb50 <- poly2nb(grid50, row.names=grid50$grid_id); nb50
is.symmetric.nb(nb50, verbose = FALSE, force = TRUE)
nb2INLA("./shapefiles/nb50.graph", nb50)


#####
# MIDWEST: Estimate abundance trends with INLA

library(rgdal)
library(sp)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(rgeos)
library(maptools)
library(spdep)
library(inlabru)
library(INLA)
library(brinla)

options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

# plot theme
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# map theme
theme_map <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.background = element_blank(),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text.align=1,
          legend.text = element_text(size=9),
          legend.title=element_text(hjust=0, size=11),
          legend.justification=c(0, 0.5),
          plot.title = element_text(size=14, hjust = 0.7))
}


# Import data
aphid_counts = read.table('./data/aphid_counts_annual_MIDWEST_m5_50km.txt',sep='\t',as.is=T,check.names=F,header=T)
aphid_sites = spTransform(readOGR(dsn='./shapefiles/aphid_sites_m5_50km.shp',layer='aphid_sites_m5_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
aphid_grid = readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F) #shapefile with grid
grid50 <- as(aphid_grid, "sf")
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors
states = spTransform(readOGR(dsn='./shapefiles/cb_2018_us_state_500k.shp',layer='cb_2018_us_state_500k',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
states2 = states[match(c('IL','IN','IA','KY','KS','LA','MO','MI','MN','ND','NE','OH','SD','WI'),states@data$STUSPS),]
eco_sf = as(states2, 'sf')

species = sort(unique(aphid_counts$Species))

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.aphids'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
for (i in 1:length(species)){

	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]

	# Subset data for species i
	sp_counts = aphid_counts[which(aphid_counts$Species==species[i]),]
	
	# Remove sites with all zeroes
	sp_sites = unique(sp_counts$Site)
	site.counts = apply(array(sp_sites),1,function(x){length(which(sp_counts$N.aphids[which(sp_counts$Site==x)]>0))})
	keep.sites = sp_sites[which(site.counts>0)]
	keep.rows = unlist(apply(array(keep.sites),1,function(x){which(sp_counts$Site==x)}))
	sp_counts = sp_counts[keep.rows,]
	sp.summary$Total.N.aphids[i] = sum(sp_counts$N.aphids)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]

	# Create spatial dataframe
	sp_counts.shp = spTransform(SpatialPointsDataFrame(coords=cbind(sp_counts$Longitude,sp_counts$Latitude),data=sp_counts,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')),CRS(proj1))
	sp_sites = unique(sp_counts$Site)
	sp_circles = aphid_sites[which(!is.na(match(aphid_sites$Site,sp_sites))),]
	
	# Summarize circles per cell (used for mapping later)
	gids = sort(unique(sp_counts.shp@data$grid_id))
	gids_count = apply(array(gids),1,function(x){gd = sp_counts.shp[which(sp_counts.shp$grid_id==x),];length(unique(gd$Site))})
	circles_per_cell = data.frame('grid_id'=gids,'number_circles'=gids_count)

	# Transform effort & standardize years
	sp_counts$ln.Trap.days = log(as.numeric(sp_counts$Trap.days))
	sp_counts$std_Year = sp_counts$Year - max(sp_counts$Year)

	# Index and sort
	sp_counts$eps_i <- sp_counts$grid_id
	sp_counts$alpha_i = sp_counts$grid_id
	sp_counts$tau_i = sp_counts$grid_id
	sp_counts$kappa_k = as.integer(factor(sp_counts$Site))
	sp_counts = arrange(sp_counts, grid_id, std_Year)
	n_circs = max(sp_counts$kappa_k, na.rm=T)
	n_cells = max(sp_counts$alpha_i, na.rm=T)

	# Make negative binomial model with raw counts
	form1 <- N.aphids ~ -1 + # remove grand mean
	  # cell ICAR random intercepts
	  f(alpha_i, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # cell ICAR random effort slopes
	  f(eps_i, ln.Trap.days, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # cell ICAR random year slopes
	  f(tau_i, std_Year, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # random circle intercepts
	  f(kappa_k, model="iid", constr=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

	# Run model
	out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3)

	# Overdispersion parameter
	sm1 = summary(out1)[[3]]
	sp.summary$Overdispersion[i] = sm1[1,4] # > 1 implies overdispersion relative to a poisson distribution. Justifies use of the negative binomial (formula: exp(-log(1/sm1[1,4]))

	# Vector of grid IDs of cells with counts
	cells_with_counts = unique(sp_counts$grid_id)

	# get alpha summaries
	alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
	alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
	alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
	alph_iw <- alph_ul - alph_ll
#	par(mfrow=c(1,3))
#	hist(alph); summary(alph)
#	hist(alph_ll); summary(alph_ll)
#	hist(alph_ul); summary(alph_ul)

	# get epsilon summaries
	eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
	eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
	eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
	eps_iw <- eps_ul - eps_ll
#	par(mfrow=c(1,3))
#	hist(eps); summary(eps); round(sum(eps<1)/length(eps), 2)
#	hist(eps_ll); summary(eps_ll)
#	hist(eps_ul); summary(eps_ul)

	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c1 = cor.test(eps, alph, method="spearman")
		sp.summary$Epsilon_Alpha_rho[i] = c1[[4]]
		sp.summary$Epsilon_Alpha_rho_p[i] = c1[[3]]
	}
	
	# get tau summaries
	tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
	tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
	tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
	tau_iw <- tau_ul - tau_ll
#	par(mfrow=c(2,2))
#	hist(tau); summary(tau); round(sum(tau>=0)/length(tau), 2)
#	hist(tau_ll); summary(tau_ll)
#	hist(tau_ul); summary(tau_ul)
#	hist(tau_iw); summary(tau_iw)
	
	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c2 = cor.test(alph, tau, method="spearman")
		sp.summary$Alpha_Tau_rho[i] = c2[[4]]
		sp.summary$Alpha_Tau_rho_p[i] = c2[[3]]
	}

	# Visualize goodness of fit with Probability Integral Transform
#	sum(out1$cpo$failure, na.rm=T) -2 * sum(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
	pit1 <- data.frame(PIT=out1$cpo$pit) %>%
	  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
	pit2 <- ggplot(data=pit1, aes(x=PIT)) +
	  geom_histogram(col="white") +
	  xlab("Probability integral transform (PIT)") +
	  ylab("Count"); pit2; summary(pit1$PIT)
	ggsave(paste0(species[i],'_pit_50km.png'),plot=pit2,device='png',path='./inla_model_output/PIT/MW/m5/',width=8,height=8,units="in",dpi=300)

	# collect posterior summaries into one dataframe
	post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,eps, eps_ll, eps_ul, eps_iw, eps_sig=NA,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
	post_sum$eps_sig <- ifelse((post_sum$eps_ll < 0 & post_sum$eps_ul > 0),post_sum$eps_sig <- NA,post_sum$eps_sig <- post_sum$eps)
#	post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
	post_sum$tau_sig <- ifelse((post_sum$tau_ll < 0 & post_sum$tau_ul > 0),post_sum$tau_sig <- NA,post_sum$tau_sig <- post_sum$tau)
	# Add medians of slope estimates to data frame
	sm2 = summary(post_sum)
	sm2[3,] = trimws(gsub('Median :','',sm2[3,]),which='both')
	sp.summary$mdn_Alph[i] = as.numeric(sm2[3,2])
	sp.summary$mdn_Alph_ll[i] = as.numeric(sm2[3,3])
	sp.summary$mdn_Alph_ul[i] = as.numeric(sm2[3,4])
	sp.summary$mdn_Alph_iw[i] = as.numeric(sm2[3,5])
	sp.summary$mdn_Eps[i] = as.numeric(sm2[3,6])
	sp.summary$mdn_Eps_ll[i] = as.numeric(sm2[3,7])
	sp.summary$mdn_Eps_ul[i] = as.numeric(sm2[3,8])
	sp.summary$mdn_Eps_iw[i] = as.numeric(sm2[3,9])
	sp.summary$mdn_Tau[i] = as.numeric(sm2[3,11])
	sp.summary$mdn_Tau_ll[i] = as.numeric(sm2[3,12])
	sp.summary$mdn_Tau_ul[i] = as.numeric(sm2[3,13])
	sp.summary$mdn_Tau_iw[i] = as.numeric(sm2[3,14])

	# Make cell level maps
	results_cells <- merge(aphid_grid, post_sum, by="grid_id", all=F)
	write.table(results_cells@data,paste0('./inla_output/MW/',species[i],'_trends_50km.txt'),sep='\t',quote=F,row.names=F) # Write model output to file
	res_sf <- as(results_cells, "sf")

	# map tau
	tau_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray20",size=1.2) +
	  geom_sf(data=res_sf, aes(fill=tau), col='gray40', size=0.1) +
	  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
						   mid = "white",
						   high = ("royalblue4"), midpoint = 0, space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# map epsilon
	eps_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
	  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
	  scale_fill_gradient2("Sampling\neffort", low = muted("purple4"), mid = "white",
						   high = muted("green4"), midpoint = median(res_sf$eps), space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# map alpha
	alph_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray20",size=1.2) +
	  geom_sf(data=res_sf, aes(fill=alph), col='gray40', size=0.1) +
	  scale_fill_gradient2("Relative\nabund.", low = "tan4", mid = "white",
						   high = "green4", midpoint = (max(res_sf$alph) + min(res_sf$alph)) / 2, space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# print cell maps
	ggsave(paste0(species[i],'_effort_50km.png'),plot =eps_p1,device='png',path='./inla_model_output/maps/effort/MW/m5',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_change_per_year_50km.png'),plot =tau_p1,device='png',path='./inla_model_output/maps/change_per_year/MW/m5',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_relative_abundance_50km.png'),plot =alph_p1,device='png',path='./inla_model_output/maps/relative_abundance/MW/m5',width=8,height=8,units="in",dpi=300)
	
}
write.table(sp.summary,'./data/aphid_trends_summary_MIDWEST_m5_50km.txt',sep='\t',quote=F,row.names=F)


###############################################################
# Idaho: Prepare data for INLA

library(lubridate)
library(rgdal)
library(spdep)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs84 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# Import aphid data
STN = read.table('./data/aphid_counts_annual_NORTHWEST_m5.txt',sep='\t',as.is=T,check.names=F,header=T)

sites = unique(STN$Site)
out = data.frame('Site'=sites,'Lon'=STN$Lon[match(sites,STN$Site)],'Lat'=STN$Lat[match(sites,STN$Site)],stringsAsFactors=F)
out = out[which(!is.na(out$Lon)),]
shp = spTransform(SpatialPointsDataFrame(coords=cbind(as.numeric(out$Lon),as.numeric(out$Lat)),data=out,proj4string=CRS(wgs84)),CRS(proj1))

# Create 50km gridded dataset
grid50 = readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F)
bover = over(shp,grid50)
shp@data$grid_id = bover$grid_id
writeOGR(shp,dsn='./shapefiles/PNW_aphid_sites_50km.shp',layer='PNW_aphid_sites_50km',driver='ESRI Shapefile',overwrite=T)
gids = sort(unique(shp@data$grid_id))
STN$grid_id = NA
for (g in 1:length(gids)){
	site1 = shp@data$Site[which(shp@data$grid_id==gids[g])]
	for (i in 1:length(site1)){
		row.pos = which(STN$Site==site1[i])
		STN$grid_id[row.pos] = gids[g]	
	}
}
write.table(STN,'./curated_data/aphid_counts_annual_NORTHWEST_m5_50km.txt',sep='\t',row.names=F)
#nb50 <- poly2nb(grid50, row.names=grid50$grid_id); nb50
#is.symmetric.nb(nb50, verbose = FALSE, force = TRUE)
#nb2INLA("./shapefiles/nb50.graph", nb50)


#####
# Idaho: Estimate abundance trends with INLA

library(rgdal)
library(sp)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(rgeos)
library(maptools)
library(spdep)
library(inlabru)
library(INLA)
library(brinla)

options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

# plot theme
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# map theme
theme_map <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.background = element_blank(),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text.align=1,
          legend.text = element_text(size=9),
          legend.title=element_text(hjust=0, size=11),
          legend.justification=c(0, 0.5),
          plot.title = element_text(size=14, hjust = 0.7))
}


# Import data
aphid_counts = read.table('./data/PNW_aphid_data_gridded_m5_50km.txt',sep='\t',as.is=T,check.names=F,header=T)
aphid_counts = aphid_counts[which(!is.na(aphid_counts$Latitude)),]
aphid_sites = spTransform(readOGR(dsn='./shapefiles/PNW_aphid_sites_50km.shp',layer='PNW_aphid_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
aphid_grid = readOGR(dsn='./shapefiles/grid_50km.shp',layer='grid_50km',verbose=F,stringsAsFactors=F) #shapefile with grid
grid50 <- as(aphid_grid, "sf")
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors
states = spTransform(readOGR(dsn='./shapefiles/cb_2018_us_state_500k.shp',layer='cb_2018_us_state_500k',verbose=F,stringsAsFactors=F),CRS(proj1))
states2 = states[match(c('ID'),states@data$STUSPS),]
eco_sf = as(states2, 'sf')

species = sort(unique(aphid_counts$Species)) #88 species

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.aphids'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
for (i in 1:length(species)){

	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]

	# Subset data for species i
	sp_counts = aphid_counts[which(aphid_counts$Species==species[i]),]
	
	# Remove sites with all zeroes
	sp_sites = unique(sp_counts$Site)
	site.counts = apply(array(sp_sites),1,function(x){length(which(sp_counts$N.aphids[which(sp_counts$Site==x)]>0))})
	keep.sites = sp_sites[which(site.counts>0)]
	keep.rows = unlist(apply(array(keep.sites),1,function(x){which(sp_counts$Site==x)}))
	sp_counts = sp_counts[keep.rows,]
	sp.summary$Total.N.aphids[i] = sum(sp_counts$N.aphids)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]

	# Create spatial dataframe
	sp_counts.shp = spTransform(SpatialPointsDataFrame(coords=cbind(sp_counts$Longitude,sp_counts$Latitude),data=sp_counts,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')),CRS(proj1))
	sp_sites = unique(sp_counts$Site)
	sp_circles = aphid_sites[which(!is.na(match(aphid_sites$Site,sp_sites))),]
	
	# Summarize circles per cell (used for mapping later)
	gids = sort(unique(sp_counts.shp@data$grid_id))
	gids_count = apply(array(gids),1,function(x){gd = sp_counts.shp[which(sp_counts.shp$grid_id==x),];length(unique(gd$Site))})
	circles_per_cell = data.frame('grid_id'=gids,'number_circles'=gids_count)

	# Transform effort & standardize years
	sp_counts$ln.Trap.days = log(as.numeric(sp_counts$Trap.days))
	sp_counts$std_Year = sp_counts$Year - max(sp_counts$Year)

	# Index and sort
	sp_counts$eps_i <- sp_counts$grid_id
	sp_counts$alpha_i = sp_counts$grid_id
	sp_counts$tau_i = sp_counts$grid_id
	sp_counts$kappa_k = as.integer(factor(sp_counts$Site))
	sp_counts = arrange(sp_counts, grid_id, std_Year)
	n_circs = max(sp_counts$kappa_k, na.rm=T)
	n_cells = max(sp_counts$alpha_i, na.rm=T)

	# Make negative binomial model with raw counts
	form1 <- N.aphids ~ -1 + # remove grand mean
	  # cell ICAR random intercepts
	  f(alpha_i, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # cell ICAR random effort slopes
	  f(eps_i, ln.Trap.days, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # cell ICAR random year slopes
	  f(tau_i, std_Year, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # random circle intercepts
	  f(kappa_k, model="iid", constr=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

	# Run model
	out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3)

	# Overdispersion parameter
	sm1 = summary(out1)[[3]]
	sp.summary$Overdispersion[i] = sm1[1,4] # > 1 implies overdispersion relative to a poisson distribution. Justifies use of the negative binomial (formula: exp(-log(1/sm1[1,4]))

	# Vector of grid IDs of cells with counts
	cells_with_counts = unique(sp_counts$grid_id)

	# get alpha summaries
	alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
	alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
	alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
	alph_iw <- alph_ul - alph_ll
#	par(mfrow=c(1,3))
#	hist(alph); summary(alph)
#	hist(alph_ll); summary(alph_ll)
#	hist(alph_ul); summary(alph_ul)

	# get epsilon summaries
	eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
	eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
	eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
	eps_iw <- eps_ul - eps_ll
#	par(mfrow=c(1,3))
#	hist(eps); summary(eps); round(sum(eps<1)/length(eps), 2)
#	hist(eps_ll); summary(eps_ll)
#	hist(eps_ul); summary(eps_ul)

	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c1 = cor.test(eps, alph, method="spearman")
		sp.summary$Epsilon_Alpha_rho[i] = c1[[4]]
		sp.summary$Epsilon_Alpha_rho_p[i] = c1[[3]]
	}
	
	# get tau summaries
	tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
	tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
	tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
	tau_iw <- tau_ul - tau_ll
#	par(mfrow=c(2,2))
#	hist(tau); summary(tau); round(sum(tau>=0)/length(tau), 2)
#	hist(tau_ll); summary(tau_ll)
#	hist(tau_ul); summary(tau_ul)
#	hist(tau_iw); summary(tau_iw)
	
	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c2 = cor.test(alph, tau, method="spearman")
		sp.summary$Alpha_Tau_rho[i] = c2[[4]]
		sp.summary$Alpha_Tau_rho_p[i] = c2[[3]]
	}

	# Visualize goodness of fit with Probability Integral Transform
#	sum(out1$cpo$failure, na.rm=T) -2 * sum(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
	pit1 <- data.frame(PIT=out1$cpo$pit) %>%
	  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
	pit2 <- ggplot(data=pit1, aes(x=PIT)) +
	  geom_histogram(col="white") +
	  xlab("Probability integral transform (PIT)") +
	  ylab("Count"); pit2; summary(pit1$PIT)
	ggsave(paste0(species[i],'_pit_50km.png'),plot=pit2,device='png',path='./inla_model_output/PIT/PNW/m5',width=8,height=8,units="in",dpi=300)

	# collect posterior summaries into one dataframe
	post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,eps, eps_ll, eps_ul, eps_iw, eps_sig=NA,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
	post_sum$eps_sig <- ifelse((post_sum$eps_ll < 0 & post_sum$eps_ul > 0),post_sum$eps_sig <- NA,post_sum$eps_sig <- post_sum$eps)
#	post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
	post_sum$tau_sig <- ifelse((post_sum$tau_ll < 0 & post_sum$tau_ul > 0),post_sum$tau_sig <- NA,post_sum$tau_sig <- post_sum$tau)
	# Add medians of slope estimates to data frame
	sm2 = summary(post_sum)
	sm2[3,] = trimws(gsub('Median :','',sm2[3,]),which='both')
	sp.summary$mdn_Alph[i] = as.numeric(sm2[3,2])
	sp.summary$mdn_Alph_ll[i] = as.numeric(sm2[3,3])
	sp.summary$mdn_Alph_ul[i] = as.numeric(sm2[3,4])
	sp.summary$mdn_Alph_iw[i] = as.numeric(sm2[3,5])
	sp.summary$mdn_Eps[i] = as.numeric(sm2[3,6])
	sp.summary$mdn_Eps_ll[i] = as.numeric(sm2[3,7])
	sp.summary$mdn_Eps_ul[i] = as.numeric(sm2[3,8])
	sp.summary$mdn_Eps_iw[i] = as.numeric(sm2[3,9])
	sp.summary$mdn_Tau[i] = as.numeric(sm2[3,11])
	sp.summary$mdn_Tau_ll[i] = as.numeric(sm2[3,12])
	sp.summary$mdn_Tau_ul[i] = as.numeric(sm2[3,13])
	sp.summary$mdn_Tau_iw[i] = as.numeric(sm2[3,14])

	# Make cell level maps
	results_cells <- merge(aphid_grid, post_sum, by="grid_id", all=F)
	write.table(results_cells@data,paste0('./inla_model_output/txt/PNW/m5/',species[i],'_trends_50km.txt'),sep='\t',quote=F,row.names=F) # Write model output to file
	res_sf <- as(results_cells, "sf")

	# map tau
	tau_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray20",size=1.2) +
	  geom_sf(data=res_sf, aes(fill=tau), col='gray40', size=0.1) +
	  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
						   mid = "white",
						   high = ("royalblue4"), midpoint = 0, space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# map epsilon
	eps_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
	  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
	  scale_fill_gradient2("Sampling\neffort", low = muted("purple4"), mid = "white",
						   high = muted("green4"), midpoint = median(res_sf$eps), space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# map alpha
	alph_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray20",size=1.2) +
	  geom_sf(data=res_sf, aes(fill=alph), col='gray40', size=0.1) +
	  scale_fill_gradient2("Relative\nabund.", low = "tan4", mid = "white",
						   high = "green4", midpoint = (max(res_sf$alph) + min(res_sf$alph)) / 2, space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# print cell maps
	ggsave(paste0(species[i],'_effort_50km.png'),plot =eps_p1,device='png',path='./inla_model_output/maps/effort/PNW/m5',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_change_per_year_50km.png'),plot =tau_p1,device='png',path='./inla_model_output/maps/change_per_year/PNW/m5',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_relative_abundance_50km.png'),plot =alph_p1,device='png',path='./inla_model_output/maps/relative_abundance/PNW/m5',width=8,height=8,units="in",dpi=300)
	
}
write.table(sp.summary,'./data/aphid_trends_summary_NORTHWEST_m5_50km.txt',sep='\t',quote=F,row.names=F)


###############################################################
# UK: Prepare data for INLA

library(rgdal)
library(spdep)

proj1 = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs' #EPSG:27700, OSGB 1936 / British National Grid -- United Kingdom Ordnance Survey
wgs84 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

aphids = read.table('./data/aphid_counts_annual_UK_m5.txt',sep='\t',as.is=T,check.names=F,header=T)

sites = unique(aphids$Site)
out = data.frame('Site'=sites,'Lon'=aphids$Longitude[match(sites,aphids$Site)],'Lat'=aphids$Latitude[match(sites,aphids$Site)],stringsAsFactors=F)
shp = spTransform(SpatialPointsDataFrame(coords=cbind(as.numeric(out$Lon),as.numeric(out$Lat)),data=out,proj4string=CRS(wgs84)),CRS(proj1))

# Create 50km gridded dataset
grid50 = spTransform(readOGR(dsn='./shapefiles/UK_grid_50km.shp',layer='UK_grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #created in QGIS with "create grid" and GBR_adm0_EPSG27700.shp as input
grid50@data$grid_id = grid50@data$id
writeOGR(grid50,dsn='./shapefiles/UK_grid_50km.shp',layer='UK_grid_50km',driver='ESRI Shapefile',overwrite=T)
nb50 <- poly2nb(grid50, row.names=grid50$grid_id); nb50
is.symmetric.nb(nb50, verbose = FALSE, force = TRUE)
nb2INLA("./shapefiles/UK_nb50.graph", nb50)

bover = over(shp,grid50)
shp@data$grid_id = bover$grid_id
writeOGR(shp,dsn='./shapefiles/aphid_sites_50km_UK.shp',layer='aphid_sites_50km_UK',driver='ESRI Shapefile',overwrite=T)
gids = sort(unique(shp@data$grid_id))
aphids$grid_id = NA
for (g in 1:length(gids)){
	site1 = shp@data$Site[which(shp@data$grid_id==gids[g])]
	for (i in 1:length(site1)){
		row.pos = which(aphids$Site==site1[i])
		aphids$grid_id[row.pos] = gids[g]	
	}
}
write.table(aphids,'./data/aphid_counts_annual_UK_m5_50km.txt',sep='\t',row.names=F)


#####
# UK: Estimate abundance trends with INLA

library(rgdal)
library(sp)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(rgeos)
library(maptools)
library(spdep)
library(inlabru)
library(INLA)
library(brinla)

options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

proj1 = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs' #EPSG:27700, OSGB 1936 / British National Grid -- United Kingdom Ordnance Survey

# plot theme
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# map theme
theme_map <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.background = element_blank(),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text.align=1,
          legend.text = element_text(size=9),
          legend.title=element_text(hjust=0, size=11),
          legend.justification=c(0, 0.5),
          plot.title = element_text(size=14, hjust = 0.7))
}


# Import data
aphid_counts = read.table('./data/aphid_counts_annual_UK_m5_50km.txt',sep='\t',as.is=T,check.names=F,header=T)
aphid_counts$N.aphids = ceiling(aphid_counts$N.aphids)
aphid_sites = spTransform(readOGR(dsn='./shapefiles/aphid_sites_50km_UK.shp',layer='aphid_sites_50km_UK',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
aphid_grid = spTransform(readOGR(dsn='./shapefiles/UK_grid_50km.shp',layer='UK_grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
grid50 <- as(aphid_grid, "sf")
g50 <- inla.read.graph('./shapefiles/UK_nb50.graph') #graph network of grid neighbors
uk.outline = spTransform(readOGR(dsn='./shapefiles/GBR_adm0_EPSG27700.shp',layer='GBR_adm0_EPSG27700',verbose=F,stringsAsFactors=F),CRS(proj1))
eco_sf = as(uk.outline, 'sf')

species = sort(unique(aphid_counts$Species))
# Check
scount = apply(array(species),1,function(x){length(which(aphid_counts$Species==x))})
singletons = species[which(scount<5)]; singletons #should be none

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.aphids'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
notrends = data.frame('Site'=NA,'Species'=NA,'tau'=0)
nx = 1
for (i in 1:length(species)){

	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]

	# Subset data for species i
	sp_counts = aphid_counts[which(aphid_counts$Species==species[i]),]
	
	# Remove sites with all zeroes
	sp_sites = unique(sp_counts$Site)
	site.counts = apply(array(sp_sites),1,function(x){length(which(sp_counts$N.aphids[which(sp_counts$Site==x)]>0))})
	keep.sites = sp_sites[which(site.counts>0)]
	keep.rows = unlist(apply(array(keep.sites),1,function(x){which(sp_counts$Site==x)}))
	sp_counts = sp_counts[keep.rows,]
	sp.summary$Total.N.aphids[i] = sum(sp_counts$N.aphids)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]
	# Deal with species*sites with same count across years
	ucounts = apply(array(sp_sites),1,function(x){length(unique(sp_counts$N.aphids[which(sp_counts$Site==x)]))})
	rm.sites = sp_sites[which(ucounts==1)]
	if (length(rm.sites)>0){
		rm.pos = as.numeric(unlist(apply(array(rm.sites),1,function(x){which(sp_counts$Site==x)})))
		sp_counts = sp_counts[-rm.pos,] #remove these sites because INLA cannot handle them
		for (nt in 1:length(rm.sites)){
			notrends[nx,1] = rm.sites[nt]
			notrends[nx,2] = species[i]
			notrends[nx,3] = 0
			nx = nx + 1
		}
	} else {
	}

	# Create spatial dataframe
	sp_counts.shp = spTransform(SpatialPointsDataFrame(coords=cbind(sp_counts$Longitude,sp_counts$Latitude),data=sp_counts,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')),CRS(proj1))
	sp_sites = unique(sp_counts$Site)
	sp_circles = aphid_sites[which(!is.na(match(aphid_sites$Site,sp_sites))),]
	
	# Summarize circles per cell (used for mapping later)
	gids = sort(unique(sp_counts.shp@data$grid_id))
	gids_count = apply(array(gids),1,function(x){gd = sp_counts.shp[which(sp_counts.shp$grid_id==x),];length(unique(gd$Site))})
	circles_per_cell = data.frame('grid_id'=gids,'number_circles'=gids_count)

	# Transform effort & standardize years
	sp_counts$ln.Trap.days = log(as.numeric(sp_counts$Trap.days))
	sp_counts$std_Year = sp_counts$Year - max(sp_counts$Year)

	# Index and sort
	sp_counts$eps_i <- sp_counts$grid_id
	sp_counts$alpha_i = sp_counts$grid_id
	sp_counts$tau_i = sp_counts$grid_id
	sp_counts$kappa_k = as.integer(factor(sp_counts$Site))
	sp_counts = arrange(sp_counts, grid_id, std_Year)
	n_circs = max(sp_counts$kappa_k, na.rm=T)
	n_cells = max(sp_counts$alpha_i, na.rm=T)

	# Make negative binomial model with raw counts
	form1 <- N.aphids ~ -1 + # remove grand mean
	  # cell ICAR random intercepts
	  f(alpha_i, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # cell ICAR random effort slopes
	  f(eps_i, ln.Trap.days, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # cell ICAR random year slopes
	  f(tau_i, std_Year, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
	  # random circle intercepts
	  f(kappa_k, model="iid", constr=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

	# Run model
	out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3, verbose=TRUE)

	# Overdispersion parameter
	sm1 = summary(out1)[[3]]
	sp.summary$Overdispersion[i] = sm1[1,4] # > 1 implies overdispersion relative to a poisson distribution. Justifies use of the negative binomial (formula: exp(-log(1/sm1[1,4]))

	# Vector of grid IDs of cells with counts
	cells_with_counts = unique(sp_counts$grid_id)

	# get alpha summaries
	alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
	alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
	alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
	alph_iw <- alph_ul - alph_ll
#	par(mfrow=c(1,3))
#	hist(alph); summary(alph)
#	hist(alph_ll); summary(alph_ll)
#	hist(alph_ul); summary(alph_ul)

	# get epsilon summaries
	eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
	eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
	eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
	eps_iw <- eps_ul - eps_ll
#	par(mfrow=c(1,3))
#	hist(eps); summary(eps); round(sum(eps<1)/length(eps), 2)
#	hist(eps_ll); summary(eps_ll)
#	hist(eps_ul); summary(eps_ul)

	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c1 = cor.test(eps, alph, method="spearman")
		sp.summary$Epsilon_Alpha_rho[i] = c1[[4]]
		sp.summary$Epsilon_Alpha_rho_p[i] = c1[[3]]
	}
	
	# get tau summaries
	tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
	tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
	tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
	tau_iw <- tau_ul - tau_ll
#	par(mfrow=c(2,2))
#	hist(tau); summary(tau); round(sum(tau>=0)/length(tau), 2)
#	hist(tau_ll); summary(tau_ll)
#	hist(tau_ul); summary(tau_ul)
#	hist(tau_iw); summary(tau_iw)
	
	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c2 = cor.test(alph, tau, method="spearman")
		sp.summary$Alpha_Tau_rho[i] = c2[[4]]
		sp.summary$Alpha_Tau_rho_p[i] = c2[[3]]
	}

	# Visualize goodness of fit with Probability Integral Transform
#	sum(out1$cpo$failure, na.rm=T) -2 * sum(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
	pit1 <- data.frame(PIT=out1$cpo$pit) %>%
	  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
	pit2 <- ggplot(data=pit1, aes(x=PIT)) +
	  geom_histogram(col="white") +
	  xlab("Probability integral transform (PIT)") +
	  ylab("Count"); pit2; summary(pit1$PIT)
	ggsave(paste0(species[i],'_pit_50km.png'),plot=pit2,device='png',path='./inla_model_output/PIT/UK',width=8,height=8,units="in",dpi=300)

	# collect posterior summaries into one dataframe
	post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,eps, eps_ll, eps_ul, eps_iw, eps_sig=NA,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
	post_sum$eps_sig <- ifelse((post_sum$eps_ll < 0 & post_sum$eps_ul > 0),post_sum$eps_sig <- NA,post_sum$eps_sig <- post_sum$eps)
	post_sum$tau_sig <- ifelse((post_sum$tau_ll < 0 & post_sum$tau_ul > 0),post_sum$tau_sig <- NA,post_sum$tau_sig <- post_sum$tau)
	# Add medians of slope estimates to data frame
	sm2 = summary(post_sum)
	sm2[3,] = trimws(gsub('Median :','',sm2[3,]),which='both')
	sp.summary$mdn_Alph[i] = as.numeric(sm2[3,2])
	sp.summary$mdn_Alph_ll[i] = as.numeric(sm2[3,3])
	sp.summary$mdn_Alph_ul[i] = as.numeric(sm2[3,4])
	sp.summary$mdn_Alph_iw[i] = as.numeric(sm2[3,5])
	sp.summary$mdn_Eps[i] = as.numeric(sm2[3,6])
	sp.summary$mdn_Eps_ll[i] = as.numeric(sm2[3,7])
	sp.summary$mdn_Eps_ul[i] = as.numeric(sm2[3,8])
	sp.summary$mdn_Eps_iw[i] = as.numeric(sm2[3,9])
	sp.summary$mdn_Tau[i] = as.numeric(sm2[3,11])
	sp.summary$mdn_Tau_ll[i] = as.numeric(sm2[3,12])
	sp.summary$mdn_Tau_ul[i] = as.numeric(sm2[3,13])
	sp.summary$mdn_Tau_iw[i] = as.numeric(sm2[3,14])

	# Make cell level maps
	results_cells <- merge(aphid_grid, post_sum, by="grid_id", all=F)
	write.table(results_cells@data,paste0('./inla_model_output/txt/UK/',species[i],'_trends_50km.txt'),sep='\t',quote=F,row.names=F) # Write model output to file
	res_sf <- as(results_cells, "sf")

	# map tau
	tau_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray20",size=1.2) +
	  geom_sf(data=res_sf, aes(fill=tau), col='gray40', size=0.1) +
	  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
						   mid = "white",
						   high = ("royalblue4"), midpoint = 0, space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# map epsilon
	eps_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
	  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
	  scale_fill_gradient2("Sampling\neffort", low = muted("purple4"), mid = "white",
						   high = muted("green4"), midpoint = median(res_sf$eps), space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# map alpha
	alph_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray20",size=1.2) +
	  geom_sf(data=res_sf, aes(fill=alph), col='gray40', size=0.1) +
	  scale_fill_gradient2("Relative\nabund.", low = "tan4", mid = "white",
						   high = "green4", midpoint = (max(res_sf$alph) + min(res_sf$alph)) / 2, space = "Lab",
						   na.value = "grey40", guide = "colourbar") +
	  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

	# print cell maps
	ggsave(paste0(species[i],'_effort_50km.png'),plot =eps_p1,device='png',path='./inla_model_output/maps/effort/UK/',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_change_per_year_50km.png'),plot =tau_p1,device='png',path='./inla_model_output/maps/change_per_year/UK/',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_relative_abundance_50km.png'),plot =alph_p1,device='png',path='./inla_model_output/maps/relative_abundance/UK/',width=8,height=8,units="in",dpi=300)
	
}
write.table(sp.summary,'./data/aphid_trends_summary_UK_50km.txt',sep='\t',quote=F,row.names=F)

