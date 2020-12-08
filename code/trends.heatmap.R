
setwd('path/to/AphidDecline')

library(ggplot2)

aphids = read.table('./data/aphids_trends_envars_traits_50km_m5_MW_ID_UK.txt',sep='\t',as.is=T,check.names=F,header=T)
aphids2 = aphids[which(!is.na(aphids$holocyclic) & aphids$hetero_mono!='V' & aphids$holocyclic!='V' & !is.na(aphids$body_length)),]

# align traits that differ among datasets, for the sake of sorting species in heatmap
aphids2$pest[which(aphids2$Species=='Hyalopterus pruni')] = 'Y'
aphids2$hetero_mono[which(aphids2$Species=='Macrosiphum euphorbiae')] = 'H'
aphids2$hetero_mono[which(aphids2$Species=='Rhopalosiphum padi')] = 'M'
aphids2$holocyclic[which(aphids2$Species=='Macrosiphum euphorbiae')] = 'Y'
aphids2$holocyclic[which(aphids2$Species=='Rhopalosiphum padi')] = 'N'

aphids2$Lon = as.character(aphids2$Longitude)
lons = unique(aphids2$Lon[order(aphids2$Longitude,decreasing=T)])
species = sort(unique(aphids2$Species[which(aphids2$Dataset=='Midwest' | aphids2$Dataset=='Idaho')]))
dat = data.frame('Species'=NA,'Site'=NA,'Trend'=-999,'Hetero'=NA,'Pest'=NA,'Holo'=NA)
dx = 1
for (i in 1:length(species)){
	for (j in 1:length(lons)){
		check = which(aphids2$Species==species[i] & aphids2$Lon==lons[j])
		if (length(check)==0){
			trend = hetero = holo = pest = NA
		} else {
			trend = aphids2$Abundance.trend[check]
			hetero = aphids2$hetero_mono[check]
			pest = aphids2$pest[check]
			holo = aphids2$holocyclic[check]
		}
		dat[dx,1] = species[i]
		dat[dx,2] = lons[j]
		dat[dx,3] = trend
		dat[dx,4] = hetero
		dat[dx,5] = pest
		dat[dx,6] = holo
		dx = dx + 1
	}
}
dat$Hetero_Pest_Holo_Species = paste(dat$Hetero,dat$Pest,dat$Holo,dat$Species,sep='_')
dat$Hetero_Pest_Holo_Species[grep('NA_NA',dat$Hetero_Pest_Holo_Species)] = NA
fac = data.frame('SHP'=sort(unique(dat$Hetero_Pest_Holo_Species)),'Pos'=1:length(which(!is.na(unique(dat$Hetero_Pest_Holo_Species)))),'Species'=apply(array(sort(unique(dat$Hetero_Pest_Holo_Species))),1,function(x){strsplit(x,'_')[[1]][4]}))
sps = unique(fac$Species)
sps.count = apply(array(sps),1,function(x){length(which(fac$Species==x))}); sps[which(sps.count>1)]
dat$Species2 = factor(dat$Species,levels=fac$Species)
dat$Site2 = as.factor(as.numeric(dat$Site))


# Bin abundance trends
qs = c(-50,-30,-10,0,10,30,50,70)
qsl = c('-50 - -31','-30 - -11','-10 - -1','0 - 9','10 - 29','30 - 49','50 - 70')
dat$Trend.bin = NA
for (i in 1:(length(qs)-1)){
	dat$Trend.bin[which(dat$Trend >= qs[i] & dat$Trend < qs[i+1])] = qsl[i]
}
dat$Trend.bin = factor(dat$Trend.bin,levels=qsl)


# Heatmap (tall)
rwb <- colorRampPalette(colors = c("orangered2", "lightyellow1", "blue"))(length(qsl))
textcol <- "black"
p <- ggplot(dat,aes(x=Site2,y=Species2,fill=Trend.bin))+
  geom_tile(colour="black",size=0.2)+
  guides(fill=guide_legend(title="Abundance\ntrend %/yr"))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),breaks=species)+
  scale_fill_manual(values=rwb,na.value = "grey20")+
  geom_vline(xintercept=16.5,colour='white') +
  geom_vline(xintercept=58.5,colour='white') +
  theme_grey(base_size=10)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol,size=30,face='bold'),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=30),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.8,"cm"),
        axis.text.x=element_text(size=15,colour=textcol,angle=-45,hjust=0,vjust=0.5),
        axis.text.y=element_text(size=25,vjust=0.2,colour=textcol),
		axis.title=element_text(size=30,face="bold"),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.5,0.4,0.4,"cm"))
ggsave(p,path='./plots/',filename="aphid trends heatmap_tall.png",height=40,width=20,units="in",dpi=300)

