

library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(performance)
library(ggplot2)
library(scales) #for alpha()
library(pdp) #for grid.arrange()
library(performance)
library(sjPlot) #for plot_model()
library(visreg) #for visreg

# Model abundance trends as a function of aphid traits and environmental variables

# Import data
aphids = read.table('./curated_data/aphids_trends_envars_traits_50km_m5_MW_ID_UK.txt',sep='\t',as.is=T,check.names=F,header=T)
aphids2 = aphids[which(!is.na(aphids$holocyclic) & aphids$hetero_mono!='V' & aphids$holocyclic!='V' & !is.na(aphids$body_length)),]

length(unique(aphids2$Species))
length(unique(aphids2$grid_id))
nrow(aphids2)
mw = aphids2[which(aphids2$Dataset=='Midwest'),] 
id = aphids2[which(aphids2$Dataset=='Idaho'),]
uk = aphids2[which(aphids2$Dataset=='UK'),]
nrow(mw); length(which(mw$Abundance.trend>0)); length(which(mw$Abundance.trend<0)); 212/1423
nrow(id); length(which(id$Abundance.trend>0)); length(which(id$Abundance.trend<0)); 275/474
nrow(uk); length(which(uk$Abundance.trend>0)); length(which(uk$Abundance.trend<0)); 1002/1075

# Z-transform continuous variables and ensure that categorical variables are factors
aphids2$pre.modern = (aphids2$pre.modern - mean(aphids2$pre.modern)) / sd(aphids2$pre.modern) #average cumulative precipitation during sampling period
aphids2$tmp.modern = (aphids2$tmp.modern - mean(aphids2$tmp.modern)) / sd(aphids2$tmp.modern) #average temperature during sampling period
aphids2$pre.trend = (aphids2$pre.trend - mean(aphids2$pre.trend)) / sd(aphids2$pre.trend) #trend in cumulative precipitation 1901-sampling period
aphids2$tmp.trend = (aphids2$tmp.trend - mean(aphids2$tmp.trend)) / sd(aphids2$tmp.trend) #trend in average precipitation 1901-sampling period
aphids2$PLAND = (aphids2$PLAND - mean(aphids2$PLAND)) / sd(aphids2$PLAND) #proportion cropland in grid cell - note different datasource than for USA
aphids2$Species = as.factor(aphids2$Species)
aphids2$Genus = as.factor(apply(array(aphids2$Genus),1,function(x){strsplit(x,' ')[[1]][1]}))
aphids2$Subfamily = as.factor(apply(array(aphids2$Subfamily),1,function(x){strsplit(x,' ')[[1]][1]}))
aphids2$grid_id = as.factor(aphids2$grid_id)
aphids2$Dataset = as.factor(aphids2$Dataset)
aphids2$hetero_mono = factor(aphids2$hetero_mono,levels=c('M','H')) #heteroecious or monoecious?
aphids2$host_breadth = as.factor(aphids2$host_breadth2) #monophagous vs not monophagous
aphids2$holocyclic = as.factor(aphids2$holocyclic) #holociclic or anholocyclic?
aphids2$pest = as.factor(aphids2$pest) #pest or non-pest?
aphids2$body_length = (aphids2$body_length - mean(aphids2$body_length)) / sd(aphids2$body_length) #median body length of apterae (in rare cases, of alatae, where apterae measurements are missing)
aphids2$pos <- numFactor(scale(aphids2$Longitude), scale(aphids2$Latitude)) # create a numeric factor recording the coordinates of the sampled locations
aphids2$ID <- factor(rep(1, nrow(aphids2))) # create a dummy group factor to be used as a random term

# Test for spatial autocorrelation in residuals
#mod1 = glmmTMB(Abundance.trend ~ 
#	holocyclic + hetero_mono + pest + host_breadth + body_length +
#	PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
#	(1 | Genus/Species) + (1 | grid_id) + (1 | Dataset), 
#	data=aphids2, REML=F)
#res = simulateResiduals(mod1)
#res2 = recalculateResiduals(res, group=aphids2$grid_id)
#testSpatialAutocorrelation(res2, x=aggregate(aphids2$Longitude,list(aphids2$grid_id),mean)$x, y=aggregate(aphids2$Latitude,list(aphids2$grid_id),mean)$x)
# significant evidence of spatial autocorrelation.

# Find optimal global model taxonomy random effect structure
mod1 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod2 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Genus) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod3 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Subfamily) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod4 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod5 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Subfamily/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod6 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Subfamily/Genus) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod7 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Subfamily/Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
a1 = AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7); a1[order(a1$AIC,decreasing=F),]
# global model with Genus/Species is AIC-best

mod1 = glmmTMB(Abundance.trend ~ holocyclic + hetero_mono + pest + host_breadth + body_length + PLAND + pre.modern + pre.trend + tmp.modern + tmp.trend + 
	(1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)

# Find AIC-best model using dredge
dd1 <- dredge(mod1)
subset(dd1, delta < 2)
write.csv(dd1,"./AICcdredge_abundance_glmmTMB_all3_genus.csv", row.names = FALSE)

# Write AIC-best models to text file
all3 = read.csv('./AICcdredge_abundance_glmmTMB_all3_genus.csv',as.is=T,check.names=F,header=T)
all3 = all3[which(all3$delta<2),3:12]; colnames(all3) = gsub('\\)','',gsub('cond\\(','',colnames(all3)))
con1 = file('AIC_best_models_all3_genus.txt','w')
for (i in 1:nrow(all3)){
	begin1 = paste0('mod',i,' = glmmTMB(Abundance.trend ~ ')
	end1 = ' + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)'
	vars =  colnames(all3)[which(!is.na(all3[i,]))]
	middle1 = paste0(vars,collapse=' + ')
	writeLines(paste0(begin1,middle1,end1),con1)
}
writeLines(paste0("modav = model.avg(list(",paste0(paste("mod",seq(1:nrow(all3)),sep=''),collapse=','),"))"),con1)
close(con1)

mod1 = glmmTMB(Abundance.trend ~ hetero_mono + holocyclic + pest + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod2 = glmmTMB(Abundance.trend ~ hetero_mono + holocyclic + pest + pre.modern + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod3 = glmmTMB(Abundance.trend ~ hetero_mono + holocyclic + host_breadth + pest + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod4 = glmmTMB(Abundance.trend ~ body_length + hetero_mono + holocyclic + pest + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod5 = glmmTMB(Abundance.trend ~ hetero_mono + holocyclic + host_breadth + pest + pre.modern + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod6 = glmmTMB(Abundance.trend ~ hetero_mono + pest + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod7 = glmmTMB(Abundance.trend ~ body_length + hetero_mono + holocyclic + pest + pre.modern + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod8 = glmmTMB(Abundance.trend ~ hetero_mono + pest + pre.modern + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod9 = glmmTMB(Abundance.trend ~ hetero_mono + holocyclic + pest + tmp.modern + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
mod10 = glmmTMB(Abundance.trend ~ hetero_mono + holocyclic + pest + PLAND + tmp.trend + (1 | Genus/Species) + (1 | grid_id) + (1 | Dataset) + exp(pos + 0 | ID), data=aphids2, REML=F)
modav = model.avg(list(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10))
check_collinearity(mod1) #check VIF in top model
summary(modav)
confint(modav, level = 0.95) #95% confidence intervals


# Histogram of trends
png('./plots/abundance trends hist.png',res=300,height=480*6,width=480*3)
par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(5,5,1,1))
hist(aphids2$Abundance.trend[which(aphids2$Dataset=='Idaho')],breaks=seq(-60,70,5),cex.lab=2,cex.main=2,lty='blank',main='Idaho',xlab='Abundance trend %/yr',ylab='No. species*grid',xaxt='n',yaxt='n')
abline(v=0,col='black',lty=2); axis(1,lwd=3,cex.axis=1.5); axis(2,lwd=3,cex.axis=1.5)
hist(aphids2$Abundance.trend[which(aphids2$Dataset=='Midwest')],breaks=seq(-60,70,5),cex.lab=2,cex.main=2,lty='blank',main='Midwest',xlab='Abundance trend %/yr',ylab='No. species*grid',xaxt='n',yaxt='n')
abline(v=0,col='black',lty=2); axis(1,lwd=3,cex.axis=1.5); axis(2,lwd=3,cex.axis=1.5)
hist(aphids2$Abundance.trend[which(aphids2$Dataset=='UK')],breaks=seq(-60,70,5),cex.lab=2,cex.main=2,lty='blank',main='UK',xlab='Abundance trend %/yr',ylab='No. species*grid',xaxt='n',yaxt='n')
abline(v=0,col='black',lty=2); axis(1,lwd=3,cex.axis=1.5); axis(2,lwd=3,cex.axis=1.5)
dev.off()

png('./plots/abundance trends hist2.png',res=300,height=480*1.5,width=480*6)
par(mfrow=c(1,3),oma=c(0,0,0,0),mar=c(5,5,1,1))
hist(aphids2$Abundance.trend[which(aphids2$Dataset=='Idaho')],breaks=seq(-60,70,5),cex.lab=2,cex.main=2,lty='blank',main='Idaho',xlab='Abundance trend %/yr',ylab='No. species*grid',xaxt='n',yaxt='n')
abline(v=0,col='black',lty=2); axis(1,lwd=3,cex.axis=1.5); axis(2,lwd=3,cex.axis=1.5)
hist(aphids2$Abundance.trend[which(aphids2$Dataset=='Midwest')],breaks=seq(-60,70,5),cex.lab=2,cex.main=2,lty='blank',main='Midwest',xlab='Abundance trend %/yr',ylab='No. species*grid',xaxt='n',yaxt='n')
abline(v=0,col='black',lty=2); axis(1,lwd=3,cex.axis=1.5); axis(2,lwd=3,cex.axis=1.5)
hist(aphids2$Abundance.trend[which(aphids2$Dataset=='UK')],breaks=seq(-60,70,5),cex.lab=2,cex.main=2,lty='blank',main='UK',xlab='Abundance trend %/yr',ylab='No. species*grid',xaxt='n',yaxt='n')
abline(v=0,col='black',lty=2); axis(1,lwd=3,cex.axis=1.5); axis(2,lwd=3,cex.axis=1.5)
dev.off()


# Plot significant fixed effects using visreg
theme1 = theme_classic() +
	theme(axis.line = element_line(size=2)) +
	theme(axis.ticks = element_line(size=2)) + 
	theme(axis.ticks.length=unit(.25, "cm")) + 
	theme(axis.text = element_text(colour='black',size=20)) +
	theme(axis.title = element_text(colour='black',size=24)) + 
	theme(plot.margin = margin(t=15,r=5,b=5,l=5, unit="pt")) + 
	theme(plot.title = element_blank())
	
v1 = visreg(mod1,'hetero_mono',gg=T,points.par=list(pch=16,size=0.5,col=alpha('black',0.1)),line.par=list(col='black',size=2)) + 
	xlab('Host alternation') + ylab('Abundance trend %/yr') + theme1 + ylim(c(-20,20))
v2 = visreg(mod1,'pest',gg=T,points.par=list(pch=16,size=0.5,col=alpha('black',0.1)),line.par=list(col='black',size=2)) + 
	xlab('Pest status') + ylab('Abundance trend %/yr') + theme1 + ylim(c(-20,20))
v3 = visreg(mod1,'tmp.trend',gg=T,points.par=list(pch=16,size=0.5,col=alpha('black',0.1)),line.par=list(col='black',size=2)) + 
	xlab('Temp. trend sd/yr') + ylab('Abundance trend %/yr') + theme1 + ylim(c(-20,20))
v123 = grid.arrange(v1,v2,v3,nrow=1)
	ggsave(plot=v123,filename='abundance_all3.png',path='./plots/glmmTMB/',width=30,height=10,units='cm',dpi=600)

# Plot significant fixed effects using sjPlot
v0 = plot_model(mod1,line.size=2,dot.size=5,vline.color='pink',colors='black') + theme1 + xlab('Predictor') + ylab('Effect')
	ggsave(plot=v0,filename='all.png',path='./plots/glmmTMB/sjplot',width=15,height=15,units='cm',dpi=600)
v1 = plot_model(mod1,type='pred',terms=c('hetero_mono'),line.size=2,dot.size=5) + theme1 + ylim(-20,10) + xlab('Host alternation') + ylab('Abundance trend %/yr')
v2 = plot_model(mod1,type='pred',terms=c('pest'),line.size=2,dot.size=5) + theme1 + ylim(-20,10) + xlab('Pest status') + ylab('Abundance trend %/yr')
v3 = plot_model(mod1,type='pred',terms=c('tmp.trend'),line.size=2) + theme1 + ylim(-20,10) + xlab('Temp. trend sd/yr') + ylab('Abundance trend %/yr')
v123 = grid.arrange(v1,v2,v3,nrow=1)
	ggsave(plot=v123,filename='all_pred.png',path='./plots/glmmTMB/sjplot',width=35,height=15,units='cm',dpi=600)


# Plot raw relationships using violin plots
theme2 = theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')

dp1 = ggplot(aphids2, aes(x=hetero_mono, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	ylim(-60,60) +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Host alternation") + theme2
dp2 = ggplot(aphids2, aes(x=pest, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	ylim(-60,60) +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Pest status") + theme2
dps = grid.arrange(dp1,dp2,nrow=1)
ggsave(filename='./plots/glmmTMB/violin_all3.png',plot=dps,dpi=600,unit='cm',height=10,width=20)



# Traits by dataset
theme2 = theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.ticks.x = element_line(size=1),
		axis.line.x = element_line(size=1),
		axis.text.x = element_text(size=16,colour='black',angle=-45,hjust=0),
		axis.title.x = element_text(colour='black',size=20),
		axis.ticks.length.x = unit(.15, "cm")) +
	theme(legend.position='none')

aphids3 = aphids2[which(aphids2$Dataset=='Idaho'),]
dp1 = ggplot(aphids3, aes(x=hetero_mono, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	ylim(-60,60) +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Host alternation") + theme2
dp2 = ggplot(aphids3, aes(x=pest, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	ylim(-60,60) +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Pest status") + theme2

aphids3 = aphids2[which(aphids2$Dataset=='Midwest'),]
dp3 = ggplot(aphids3, aes(x=hetero_mono, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	ylim(-60,60) +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Host alternation") + theme2
dp4 = ggplot(aphids3, aes(x=pest, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	ylim(-60,60) +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Pest status") + theme2

aphids3 = aphids2[which(aphids2$Dataset=='UK'),]
dp5 = ggplot(aphids3, aes(x=hetero_mono, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	ylim(-60,60) +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Host alternation") + theme2
dp6 = ggplot(aphids3, aes(x=pest, y=Abundance.trend)) + 
	geom_violin(trim=T,size=1,scale='count',fill='grey80', color='grey80') +
	ylim(-60,60) +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=0,size=1,colour='black') +
	labs(y="Abundance trend %/yr") +
	labs(x="Pest status") + theme2

dps = grid.arrange(dp1,dp3,dp5,dp2,dp4,dp6,nrow=2)
ggsave(filename='./plots/glmmTMB/violin_all3_dataset.png',plot=dps,dpi=600,unit='cm',height=20,width=20)

