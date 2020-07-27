
#   -----------------------------------------------------------------------
# Macroalgae, Urchins, and Heat Stress Magnify Coral Loss After Coral Bleaching
# Code and analysis by Mary Donovan
# 2_summaries_and_visualizations
#   -----------------------------------------------------------------------

# This code reads in model output from 1_run_model.R and creates summaries and visualizations 

# This code was written and executed with R version 3.6.1


# Initialization ----------------------------------------------------------
library(rjags) # v.4.9 linked to JAGS 4.3.0
library(dplyr) # v.0.8.3
library(wesanderson)
library(rgdal) # v.1.4-4 linked to GDAL 2.4.2
library(ggplot2) # v.3.3.0

# inputs ------------------------------------------------------------------
# read in files created by 1_run_model.R
zmPrp <- readRDS('outputs/change_1yr_out_coda.Rdata')
mod_data <- read.csv('outputs/change_1yr_out_inputdata.csv')
locs <- read.csv('outputs/change_1yr_out_locs_all.csv')
coords <- read.csv('outputs/coord_export.csv')
delta <- read.csv('outputs/change_1yr_out_delta.csv')

# Model checks ------------------------------------------------------------
grepgo <- grep('\\bpval.mean\\b',colnames(zmPrp[[1]]))
summary(zmPrp[,grepgo])
grepgo <- grep('\\bpval.sd\\b',colnames(zmPrp[[1]]))
summary(zmPrp[,grepgo])
grepgo <- grep('\\bR2\\b',colnames(zmPrp[[1]]))
summary(zmPrp[,grepgo])

# plot predicted versus observed
grepgo <- grep('y.new',colnames(zmPrp[[1]]))
y.new_sum <- summary(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))
y.new_out <- data.frame(ynew=y.new_sum$quantiles[,'50%'])

png(file='outputs/FigureS4.png',height=1800,width=2000,res=300)
plot(y.new_out$ynew,mod_data$y,xlab='Predicted',ylab='Observed',ylim=c(0,1),xlim=c(0,1)); abline(a=0,b=1,col='red',lwd=2)
dev.off()

# plot residuals
mod_data$resid <- mod_data$y - y.new_out$ynew
hist(mod_data$resid)

# gelman diagnostics for chains
grepgo <- grep('\\bbeta_d\\b',colnames(zmPrp[[1]]))
gelman.diag(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))
gelman_supp_table <- data.frame(var=c('DHW Max','Macroalgae','Parrotfish','SSTA Freq SD','Temp Max','Depth','Exposure','Cyclones','Turbidity','Urchins','DHW Max*Macroalgae','DHW Max*Parrotfish','DHW Max*SSTA Freq SD','DHW Max*Temp Max','DHW Max*Depth','DHW Max*Exposure','DHW Max*Cyclones','DHW Max*Turbidity','DHW Max*Urchins'),
                                gel=gelman.diag(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))[[1]][,1])
grepgo <- grep('\\br\\b',colnames(zmPrp[[1]]))
gelman_supp_table <- rbind(gelman_supp_table,data.frame(var='r',
                                    gel=gelman.diag(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))[[1]][,1]))
grepgo <- grep('\\brho\\b',colnames(zmPrp[[1]]))
gelman_supp_table <- rbind(gelman_supp_table,data.frame(var='rho',
                                     gel=gelman.diag(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))[[1]][,1]))
grepgo <- grep('sigma',colnames(zmPrp[[1]]))
gelman_supp_table <- rbind(gelman_supp_table,data.frame(var=c('sigma 1','sigma 2','sigma b0'),
                                    gel=gelman.diag(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))[[1]][,1]))
write.csv(gelman_supp_table, 'outputs/gelman_table_for_supp.csv',row.names=F)

# correlation among predictors
cor(delta[c('dhw.s','macro.s','parrot.s','ssta_freq_sd.s','temp_max.s','depth.s','exposure_fac','cyclone_mag.s','turbidity.s','urchin.s')])
temp <- delta[c('dhw.s','macro.s','parrot.s','ssta_freq_sd.s','temp_max.s','depth.s','exposure_fac','cyclone_mag.s','turbidity.s','urchin.s')]
colnames(temp) <- c('DHW Max','Macroalgae','Parrotfish','SSTA Freq SD','Temp Max','Depth','Wave Exposure','Cyclones','Turbidity','Urchins')
png(file='outputs/FigureS2.png',height=3200,width=3400,res=300)
GGally::ggpairs(temp) + theme_classic() +
  theme(
    axis.ticks = element_line(colour = "black", size = 0.5),
    axis.ticks.length=unit(.25, "cm"),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.text = element_text(size=10), axis.title = element_text(size=10),
    strip.text = element_text(size=10), strip.background = element_blank()
  )
dev.off()  

# Figure 1 ----------------------------------------------------------------
# map of study sites

wlrd.p <- readOGR('data/shapefiles/TM_WORLD_BORDERS_SIMPL_PC150.shp') # from NASA (https://github.com/nasa/World-Wind-Java) and reprojected with Pacific in center

png(file='outputs/Figure1.png',height=1400,width=4300,res=300)
par(mfrow=c(1,1),mgp=c(0.5,0.6,0), mar=c(1,1,1,1))
plot(wlrd.p,ylim=c(-4400000,4400000),xlim=c(1600000,2400000), col='grey90',border='grey70')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°E','180°','60°W'),las=1,tcl=0.45,mgp=c(-1,-1.3,0),cex.axis=1.3)
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°N','0°','23°S'),las=2,tcl=0.45,mgp=c(-1,-0.6,0),hadj=0,cex.axis=1.3)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.45,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.45,mgp=c(-1,-0.6,0),hadj=0)
box()
xy <- coords
xy <- SpatialPointsDataFrame(data=xy,coords=xy[c('Long','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy,pch=21,bg='dodgerblue',cex=1.6)
dev.off()


# Figure 2 ----------------------------------------------------------------
# calculate confidence intervals
grepgo <- grep('\\bbeta_d\\b',colnames(zmPrp[[1]]))
beta_d_80 <- data.frame(
  beta = rep(NA,length(grepgo)),
  beta_up = rep(NA,length(grepgo)),
  beta_down = rep(NA,length(grepgo)),
  beta_up80 = rep(NA,length(grepgo)),
  beta_down80 = rep(NA,length(grepgo)),
  beta_up90 = rep(NA,length(grepgo)),
  beta_down90 = rep(NA,length(grepgo))
)
for(i in 1:length(grepgo)){
  beta_d_80$beta[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                  zmPrp[[2]][,grepgo[i]],
                                  zmPrp[[3]][,grepgo[i]]),0.5)
  beta_d_80$beta_up[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                     zmPrp[[2]][,grepgo[i]],
                                     zmPrp[[3]][,grepgo[i]]),0.975)
  beta_d_80$beta_down[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                       zmPrp[[2]][,grepgo[i]],
                                       zmPrp[[3]][,grepgo[i]]),0.025)
  beta_d_80$beta_up80[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                       zmPrp[[2]][,grepgo[i]],
                                       zmPrp[[3]][,grepgo[i]]),0.90)
  beta_d_80$beta_down80[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                         zmPrp[[2]][,grepgo[i]],
                                         zmPrp[[3]][,grepgo[i]]),0.10)
  # beta_d_80$beta_down60[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
  #                                        zmPrp[[2]][,grepgo[i]],
  #                                        zmPrp[[3]][,grepgo[i]]),0.15)
  # beta_d_80$beta_up60[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
  #                                      zmPrp[[2]][,grepgo[i]],
  #                                      zmPrp[[3]][,grepgo[i]]),0.85)
  beta_d_80$beta_up90[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                       zmPrp[[2]][,grepgo[i]],
                                       zmPrp[[3]][,grepgo[i]]),0.95)
  beta_d_80$beta_down90[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],
                                         zmPrp[[2]][,grepgo[i]],
                                         zmPrp[[3]][,grepgo[i]]),0.05)
}

# add label and order
beta_d_80_ordered <- beta_d_80 %>% mutate(driver_name=c('DHW Max','Macroalgae','Parrotfish','SSTA Freq SD','Temp Max','Depth','Wave Exposure','Cyclones','Turbidity','Urchins','DHW Max*Macroalgae','DHW Max*Parrotfish','DHW Max*SSTA Freq SD','DHW Max*Temp Max','DHW Max*Depth','DHW Max*Wave Exposure','DHW Max*Cyclones','DHW Max*Turbidity','DHW Max*Urchins')) %>% arrange(desc(beta))

# move dhw_max*temp_max down one so red bars are together
beta_d_80_ordered <- beta_d_80_ordered[c(1:11,14,12,17,13,15,18,16,19),]

# assign color
beta_d_80_ordered$col_use <- ifelse(
  beta_d_80_ordered$beta < 0 & round(beta_d_80_ordered$beta_up90,2) < 0, wes_palette("Zissou1")[5],
  ifelse(
    beta_d_80_ordered$beta < 0 & round(beta_d_80_ordered$beta_up80,2) < 0, rgb(207,135,127,max=255),'grey'))
beta_d_80_ordered$col_use <- ifelse(
  beta_d_80_ordered$beta > 0 & round(beta_d_80_ordered$beta_down90,2) > 0, wes_palette("Zissou1")[1],ifelse(
    beta_d_80_ordered$beta > 0 & round(beta_d_80_ordered$beta_down80,2) > 0, wes_palette("Zissou1")[2],beta_d_80_ordered$col_use))

# plot
png(file='outputs/Figure2.png',height=2000,width=2000,res=300)
par(mar=c(4,12,1,2),mfrow=c(1,1))
plot(beta_d_80_ordered$beta,seq(1:nrow(beta_d_80)),type='n',yaxt='n',ylab='',xlab=expression(gamma),xlim=c(-0.27,0.27),cex.lab=1.5)
abline(v=0,lty=2,lwd=1.5)
plotrix::plotCI(beta_d_80_ordered$beta,seq(1:nrow(beta_d_80)),ui=beta_d_80_ordered$beta_up,li=beta_d_80_ordered$beta_down,
                pch=NA,err='x',sfrac=0,col=beta_d_80_ordered$col_use,lwd=0.7,add=T)
plotrix::plotCI(beta_d_80_ordered$beta,seq(1:nrow(beta_d_80)),ui=beta_d_80_ordered$beta_up90,li=beta_d_80_ordered$beta_down90,
                err='x',yaxt='n',ylab='',pch=NA,lwd=3,add=T,sfrac=0,col=beta_d_80_ordered$col_use)
plotrix::plotCI(beta_d_80_ordered$beta,seq(1:nrow(beta_d_80)),ui=beta_d_80_ordered$beta_up80,li=beta_d_80_ordered$beta_down80,
                err='x',yaxt='n',ylab='',pch=NA,lwd=6,add=T,sfrac=0,col=beta_d_80_ordered$col_use)
axis(2,at=seq(1:nrow(beta_d_80)),labels=beta_d_80_ordered$driver_name,las=2)
dev.off()

# parameter posteriors
driver_name <- c('DHW Max','Macroalgae','Parrotfish','SSTA Freq SD','Temp Max','Depth','Wave Exposure','Cyclones','Turbidity','Urchins','DHW Max*Macroalgae','DHW Max*Parrotfish','DHW Max*SSTA Freq SD','DHW Max*Temp Max','DHW Max*Depth','DHW Max*Wave Exposure','DHW Max*Cyclones','DHW Max*Turbidity','DHW Max*Urchins')
crit <- 0.90
grepgo <- grep('\\bbeta_d\\b',colnames(zmPrp[[1]]))

png(file='outputs/FigureS1.png',height=2000,width=2700,res=300)
par(mfrow=c(4,5),mgp=c(1.7,.7,0),mar=c(3,2,1,1),oma=c(0,2,0,0))
for(i in 1:length(grepgo)){
  poster <- c(zmPrp[[1]][,grepgo[i]],zmPrp[[2]][,grepgo[i]],zmPrp[[3]][,grepgo[i]])
  plot(density(poster),xlab=driver_name[i],type='n',main='',ylab='',xlim=c(-0.4,0.4))
  abline(v=0,lty=2)
  p <- density(poster)
  polygon(p,col=rgb(211,211,211,100,max=255),border=NA)
  p_df <- data.frame(x=p$x,y=p$y)
  p_df_80 <- p_df %>% filter(x > rethinking::HPDI(poster,prob=crit)[1] & x < rethinking::HPDI(poster,prob=crit)[2])
  p_df_80 <- rbind(data.frame(x=min(p_df_80$x),y=0),p_df_80); p_df_80 <- rbind(p_df_80,data.frame(x=max(p_df_80$x),y=0))
  polygon(p_df_80,col=rgb(211,211,211,200,max=255),border=NA)
  points(p,col='darkgrey',lwd=1,type='l')
}
mtext('Posterior',side=2,outer=T,cex=1.1)
dev.off()

# Figure 3 ----------------------------------------------------------------

# Estimates of beta_hat
grepgo <- grep('B.hat',colnames(zmPrp[[1]]))
B_hat_sum <- summary(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))
B_hat_out <- data.frame(Bhat=B_hat_sum$quantiles[,'50%'],
                        Bhat_up=B_hat_sum$quantiles[,'97.5%'],
                        Bhat_down=B_hat_sum$quantiles[,'2.5%'])

# partial predictions
pred_out <- function(pred_name){
  grepgo <- grep(paste0('\\b',pred_name,'\\b'),colnames(zmPrp[[1]]))
  out <- data.frame(mean=rep(NA,length(grepgo))); out$up <- NA; out$down <- NA; out$up80 <- NA; out$down80 <- NA
  for(i in 1:length(grepgo)){
    out$mean[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],zmPrp[[2]][,grepgo[i]],zmPrp[[3]][,grepgo[i]]),0.5)
    out$up[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],zmPrp[[2]][,grepgo[i]],zmPrp[[3]][,grepgo[i]]),0.975)
    out$down[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],zmPrp[[2]][,grepgo[i]],zmPrp[[3]][,grepgo[i]]),0.025)
    out$up80[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],zmPrp[[2]][,grepgo[i]],zmPrp[[3]][,grepgo[i]]),0.90)
    out$down80[i] <- quantile(c(zmPrp[[1]][,grepgo[i]],zmPrp[[2]][,grepgo[i]],zmPrp[[3]][,grepgo[i]]),0.10)
  }
  return(out)
}

pred_DHW_out <- pred_out(pred_name='pred_dhw')
pred_DHW_out$x <- seq(from=min(delta$dhw.s),to=max(delta$dhw.s),length=100)

pred_MA_lD_out <- pred_out(pred_name='pred_MA_lD')
pred_MA_lD_out$x <- seq(from=min(delta$macro.s),to=max(delta$macro.s),length=100)

pred_MA_hD_out <- pred_out(pred_name='pred_MA_hD')
pred_MA_hD_out$x <- seq(from=min(delta$macro.s),to=max(delta$macro.s),length=100)

pred_depth_lD_out <- pred_out(pred_name='pred_depth_lD')
pred_depth_lD_out$x <- seq(from=min(delta$depth.s),to=max(delta$depth.s),length=100)

pred_depth_hD_out <- pred_out(pred_name='pred_depth_hD')
pred_depth_hD_out$x <- seq(from=min(delta$depth.s),to=max(delta$depth.s),length=100)

pred_urchin_out <- pred_out(pred_name='pred_urchin')
pred_urchin_out$x <- seq(from=min(delta$urchin.s),to=max(delta$urchin.s),length=100)

pred_exposure_out <- pred_out(pred_name='pred_exposure')

dat <- data.frame(
  B_2 = B_hat_out$Bhat[(nrow(delta)+1):nrow(B_hat_out)],
  B_2_up = B_hat_out$Bhat_up[(nrow(delta)+1):nrow(B_hat_out)],
  B_2_down = B_hat_out$Bhat_down[(nrow(delta)+1):nrow(B_hat_out)],
  delta = delta$delta,
  delta_up = delta$delta_up,
  delta_down = delta$delta_down,
  id_fac = delta$id_fac
)
delta_lm <- lm(delta~B_2,data=dat)
ax_backtrans <- data.frame(B_2 = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6))
ax_backtrans$delta <- round(predict(delta_lm, newdata=ax_backtrans),2)

delta_lm_back <- lm(B_2~delta,data=dat)
ax_forwtrans <- data.frame(delta = c(-.15,-0.1,-0.05,0,0.05,0.1))
ax_forwtrans$B_2 <- round(predict(delta_lm_back, newdata=ax_forwtrans),2)

dhw_lm <- lm(dhw.s~dhw_max,data=delta)
dhw_ax <- data.frame(dhw_max=c(0,4,8,12,16))
dhw_ax$dhw.s <- predict(dhw_lm,newdata=dhw_ax)

png(file='outputs/Figure3.png',height=3450,width=1500,res=300)
par(mfrow=c(3,1),mar=c(4.1,2,1,1),oma=c(0,3,0,0),mgp=c(2.8,1,0))

# Macro
plot(dat$B_2~delta$macro.s,xlab='Macroalgal Cover (%)',ylim=c(-0.8,0.6)
     ,ylab='',type='n',xaxt='n',yaxt='n',cex.axis=1.7,cex.lab=1.7,bty='l')
axis(1,at=c(-0.58,0.22,1.01,1.80,2.59,3.38,4.17),labels=c(0,10,20,30,40,50,60),cex.axis=1.7)
axis(2,at=ax_forwtrans$B_2,labels=ax_forwtrans$delta,cex.axis=1.7)
abline(h=0,lty=2)
polygon(c(pred_MA_lD_out$x,rev(pred_MA_lD_out$x)),c(pred_MA_lD_out$up80,rev(pred_MA_lD_out$down80)),col=rgb(91,188,214,105,max=255),border=NA)
polygon(c(pred_MA_hD_out$x,rev(pred_MA_hD_out$x)),c(pred_MA_hD_out$up80,rev(pred_MA_hD_out$down80)),col=rgb(249,132,0,105,max=255),border=NA)
points(B_hat_out$Bhat[(nrow(delta)+1):nrow(B_hat_out)]~delta$macro.s,pch=21,bg=rgb(190,190,190,125,max=255),cex=1.2)
text(-0.20,.6,'A',cex=2,pos=2)

# Urchin
plot(dat$B_2~delta$urchin.s,ylim=c(-0.8,0.6),xlab=expression("Urchin abundance"~~bgroup("(",'100 '*m^{-2},")")),
     ylab='',type='n',xaxt='n',yaxt='n',cex.axis=1.7,cex.lab=1.7,bty='l',xlim=c(-1,3.3),mgp=c(3.2,1,0))
axis(1,at=c(-0.89,-0.48,0.54,1.87,3.23),labels=c(0,1,10,100,1000),cex.axis=1.7)
axis(2,at=ax_forwtrans$B_2,labels=ax_forwtrans$delta,cex.axis=1.7)
abline(h=0,lty=2)
polygon(c(pred_urchin_out$x,rev(pred_urchin_out$x)),c(pred_urchin_out$up80,rev(pred_urchin_out$down80)),col=rgb(190,190,190,155,max=255),border=NA)
points(B_hat_out$Bhat[(nrow(delta)+1):nrow(B_hat_out)]~delta$urchin.s,pch=21,bg=rgb(190,190,190,125,max=255),cex=1.2)
text(-0.65,.6,'B',cex=2,pos=2)

# Depth
plot(dat$B_2~delta$depth.s,xlab='Depth (m)',ylim=c(-0.8,0.6)
     ,ylab='',type='n',xaxt='n',yaxt='n',cex.axis=1.7,cex.lab=1.7,bty='l',xlim=c(-2,3.7))
axis(1,at=c(-1.88,-0.55,0.77,2.10,3.42),labels=c(0,5,10,15,20),cex.axis=1.7)
axis(2,at=ax_forwtrans$B_2,labels=ax_forwtrans$delta,cex.axis=1.7)
abline(h=0,lty=2)
polygon(c(pred_depth_lD_out$x,rev(pred_depth_lD_out$x)),c(pred_depth_lD_out$up80,rev(pred_depth_lD_out$down80)),col=rgb(91,188,214,105,max=255),border=NA)
polygon(c(pred_depth_hD_out$x,rev(pred_depth_hD_out$x)),c(pred_depth_hD_out$up80,rev(pred_depth_hD_out$down80)),col=rgb(249,132,0,105,max=255),border=NA)
points(B_hat_out$Bhat[(nrow(delta)+1):nrow(B_hat_out)]~delta$depth.s,pch=21,bg=rgb(190,190,190,125,max=255),cex=1.2)
text(-1.58,.6,'C',cex=2,pos=2)

mtext(expression(paste('Change in coral cover (',Delta,')')), outer=T,side=2,cex=1.5,line=0.5)
dev.off()


# supplemental figure DHW x Macro -----------------------------------------
temp <- data.frame(macro.s = delta$macro.s, macro = delta$macro) 
temp <- temp %>% mutate(macro_bins = cut(macro, breaks = c(-Inf,0.05,0.25,Inf)))
temp %>% group_by(macro_bins) %>% summarise(min(macro),max(macro),min(macro.s),max(macro.s),n=length(macro),median(macro.s),mean(macro.s),median(macro),mean(macro))
macro_bin <- temp %>% group_by(macro_bins) %>% summarise(binz = median(macro.s))

three_levels <- data.frame(
  macro.s = rep(macro_bin$binz[1],length(unique(delta$dhw.s[delta$macro <= 0.05]))),
  dhw.s = unique(delta$dhw.s[delta$macro <= 0.05])
)
three_levels <- rbind(three_levels,data.frame(
  macro.s = rep(macro_bin$binz[2],length(unique(delta$dhw.s[(delta$macro > 0.05 & delta$macro <= 0.25)]))),
  dhw.s = unique(delta$dhw.s[(delta$macro > 0.05 & delta$macro <= 0.25)])
))
three_levels <- rbind(three_levels,data.frame(
  macro.s = rep(macro_bin$binz[3],length(unique(delta$dhw.s[(delta$macro > 0.25)]))),
  dhw.s = unique(delta$dhw.s[(delta$macro > 0.25)])
))
three_levels$pred <- NA; three_levels$pred_up80 <- NA; three_levels$pred_down80 <- NA


# save model output mcmc iterations
grepgo <- grep('\\bmu_b0\\b',colnames(zmPrp[[1]]))
mu_b0 <- c(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo])

grepgo <- grep('\\bbeta_d\\b',colnames(zmPrp[[1]]))
beta_d <- rbind(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo])
beta_matrix <- as.matrix(beta_d)

# loop through each value and predict for the macroalgae and associated reduction, for each value of dhw

for(j in 1:nrow(three_levels)){ # for each value of dhw
  X_new <- data.frame(dhw.s = three_levels$dhw.s[j],
                      macro.s = three_levels$macro.s[j],
                      parrot.s = rep(mean(delta$parrot.s),1),
                      ssta_freq_sd.s = rep(mean(delta$ssta_freq_sd.s),1),
                      temp_max.s = rep(mean(delta$temp_max.s),1),
                      depth.s = rep(mean(delta$depth.s),1),
                      exposure_fac = rep(0,1),
                      cyclone_mag.s = rep(mean(delta$cyclone_mag.s),1),
                      turbidity.s = rep(mean(delta$turbidity.s),1),
                      urchin.s = rep(mean(delta$urchin.s),1))
  # convert to model matrix
  X_new <- as.matrix(model.matrix(~ X_new$dhw.s + X_new$macro.s + X_new$parrot.s + X_new$ssta_freq_sd.s + 
                                    X_new$temp_max.s + X_new$depth.s + X_new$exposure_fac +  X_new$cyclone_mag.s + 
                                    X_new$turbidity.s + X_new$urchin.s + X_new$dhw.s*X_new$macro.s + X_new$dhw.s*X_new$parrot.s + 
                                    X_new$dhw.s*X_new$ssta_freq_sd.s + X_new$dhw.s*X_new$temp_max.s + X_new$dhw.s*X_new$depth.s + 
                                    X_new$dhw.s*X_new$exposure_fac + X_new$dhw.s*X_new$cyclone_mag.s + X_new$dhw.s*X_new$turbidity.s + 
                                    X_new$dhw.s*X_new$urchin.s))
  X_new <- as.matrix(X_new[,2:ncol(X_new)]) # remove intercept
  mcmc_go <- mu_b0 + beta_matrix %*% X_new # predict from model output
  
  # calculate 'improvement' and calculate intervals across posteriors
  three_levels$pred[j] <- median(mcmc_go)
  three_levels$pred_up80[j] <- quantile(mcmc_go,0.90)
  three_levels$pred_down80[j] <- quantile(mcmc_go,0.10)
}



png(file='outputs/Figure_S1.png',height=3450,width=1500,res=300)
par(mfrow=c(3,1),mar=c(4.1,2,1.5,1),oma=c(1,3,1,0),mgp=c(2.8,1,0))

plot(dat$B_2[delta$macro <= 0.05] ~ delta$dhw.s[delta$macro <= 0.05],ylim=c(-1.1,0.75),xlim=c(-0.6,6),xaxt='n',type='p',col='dodgerblue',yaxt='n',xlab="",cex=1.3,cex.axis=1.7,cex.lab=1.7,bty='l',ylab='')
title('Macroalgae < 5%',adj=0.05,cex.main=2)
abline(h=0,lty=2,col='grey')
axis(1,at=dhw_ax$dhw.s,labels=dhw_ax$dhw_max,cex.axis=1.7)
axis(2,at=ax_forwtrans$B_2,labels=ax_forwtrans$delta,cex.axis=1.7)
temp <- three_levels[three_levels$macro.s==macro_bin$binz[1],]; temp <- temp %>% arrange(temp$macro.s)
plotrix::plotCI(temp$dhw.s,temp$pred,ui=temp$pred_up80,li=temp$pred_down80,err='y',sfrac=0,add=T)

plot(dat$B_2[(delta$macro > 0.05 & delta$macro <= 0.25)] ~ delta$dhw.s[(delta$macro > 0.05 & delta$macro <= 0.25)],ylim=c(-1.1,0.75),xlim=c(-0.6,6),xaxt='n',type='p',col='dodgerblue',yaxt='n',xlab="",cex=1.3,cex.axis=1.7,cex.lab=1.7,bty='l',ylab='')
title('Macroalgae 5-25% ',adj=0.05,cex.main=2)
abline(h=0,lty=2,col='grey')
axis(1,at=dhw_ax$dhw.s,labels=dhw_ax$dhw_max,cex.axis=1.7)
axis(2,at=ax_forwtrans$B_2,labels=ax_forwtrans$delta,cex.axis=1.7)
temp <- three_levels[three_levels$macro.s==macro_bin$binz[2],]; temp <- temp %>% arrange(temp$macro.s)
plotrix::plotCI(temp$dhw.s,temp$pred,ui=temp$pred_up80,li=temp$pred_down80,err='y',sfrac=0,add=T)

plot(dat$B_2[delta$macro > 0.25] ~ delta$dhw.s[delta$macro > 0.25],ylim=c(-1.1,0.75),xlim=c(-0.6,6),xaxt='n',type='p',col='dodgerblue',yaxt='n',xlab="",cex=1.3,cex.axis=1.7,cex.lab=1.7,bty='l',ylab='')
title('Macroalgae > 25% ',adj=0.05,cex.main=2)
abline(h=0,lty=2,col='grey')
axis(1,at=dhw_ax$dhw.s,labels=dhw_ax$dhw_max,cex.axis=1.7)
axis(2,at=ax_forwtrans$B_2,labels=ax_forwtrans$delta,cex.axis=1.7)
temp <- three_levels[three_levels$macro.s==macro_bin$binz[3],]; temp <- temp %>% arrange(temp$macro.s)
plotrix::plotCI(temp$dhw.s,temp$pred,ui=temp$pred_up80,li=temp$pred_down80,err='y',sfrac=0,add=T)

mtext(expression(paste('Change in coral cover (',Delta,')')), outer=T,side=2,cex=1.5,line=0.5)
mtext('DHW Max', outer=T,side=1,cex=1.5,line=-0.5)
dev.off()


# values for the text

#Even at similar levels of heat stress (e.g., 4 DHW), reefs with more macroalgae saw over 10X higher mortality of corals 
new_values <- data.frame(
  dhw.s = as.numeric(rep(predict(dhw_lm,newdata=data.frame(dhw_max=4))[1],2)),
  macro.s = c(macro_bin$binz[1],macro_bin$binz[3])
)
new_values$pred <- NA; new_values$pred_up80 <- NA; new_values$pred_down80 <- NA; new_values$pred_up95 <- NA; new_values$pred_down95 <- NA

for(j in 1:nrow(new_values)){ # for each value of dhw
  X_new <- data.frame(dhw.s = new_values$dhw.s[j],
                      macro.s = new_values$macro.s[j],
                      parrot.s = rep(mean(delta$parrot.s),1),
                      ssta_freq_sd.s = rep(mean(delta$ssta_freq_sd.s),1),
                      temp_max.s = rep(mean(delta$temp_max.s),1),
                      depth.s = rep(mean(delta$depth.s),1),
                      exposure_fac = rep(0,1),
                      cyclone_mag.s = rep(mean(delta$cyclone_mag.s),1),
                      turbidity.s = rep(mean(delta$turbidity.s),1),
                      urchin.s = rep(mean(delta$urchin.s),1))
  # convert to model matrix
  X_new <- as.matrix(model.matrix(~ X_new$dhw.s + X_new$macro.s + X_new$parrot.s + X_new$ssta_freq_sd.s + 
                                    X_new$temp_max.s + X_new$depth.s + X_new$exposure_fac +  X_new$cyclone_mag.s + 
                                    X_new$turbidity.s + X_new$urchin.s + X_new$dhw.s*X_new$macro.s + X_new$dhw.s*X_new$parrot.s + 
                                    X_new$dhw.s*X_new$ssta_freq_sd.s + X_new$dhw.s*X_new$temp_max.s + X_new$dhw.s*X_new$depth.s + 
                                    X_new$dhw.s*X_new$exposure_fac + X_new$dhw.s*X_new$cyclone_mag.s + X_new$dhw.s*X_new$turbidity.s + 
                                    X_new$dhw.s*X_new$urchin.s))
  X_new <- as.matrix(X_new[,2:ncol(X_new)]) # remove intercept
  mcmc_go <- mu_b0 + beta_matrix %*% X_new # predict from model output
  
  # calculate 'improvement' and calculate intervals across posteriors
  new_values$pred[j] <- median(mcmc_go)
  new_values$pred_up80[j] <- quantile(mcmc_go,0.90)
  new_values$pred_down80[j] <- quantile(mcmc_go,0.10)
  new_values$pred_up95[j] <- quantile(mcmc_go,0.975)
  new_values$pred_down95[j] <- quantile(mcmc_go,0.025)
}
new_values

round(new_values$pred[1]/new_values$pred[2],2)
round(new_values$pred_up80[1]/new_values$pred_up80[2],2)
round(new_values$pred_down80[1]/new_values$pred_down80[2],2)

# For example, reductions in fishing of herbivorous fishes or nutrient pollution can lead to declines in macroalgal abundance from 10-40% (34–36). Our analyses suggest that such reductions in macroalgal abundance would reduce coral mortality by XX% under a moderate heat stress event. 
macro_lm_for <- lm(macro.s~macro,data=delta)

new_values <- data.frame(
  dhw.s = as.numeric(rep(predict(dhw_lm,newdata=data.frame(dhw_max=4))[1],3)),
  macro.s = as.numeric(predict(macro_lm_for,newdata=data.frame(macro=c(0,.10,.40))))
)
new_values$pred <- NA; new_values$pred_up80 <- NA; new_values$pred_down80 <- NA; new_values$pred_up95 <- NA; new_values$pred_down95 <- NA

X_new <- data.frame(dhw.s = new_values$dhw.s[1],
                    macro.s = new_values$macro.s[1],
                    parrot.s = rep(mean(delta$parrot.s),1),
                    ssta_freq_sd.s = rep(mean(delta$ssta_freq_sd.s),1),
                    temp_max.s = rep(mean(delta$temp_max.s),1),
                    depth.s = rep(mean(delta$depth.s),1),
                    exposure_fac = rep(0,1),
                    cyclone_mag.s = rep(mean(delta$cyclone_mag.s),1),
                    turbidity.s = rep(mean(delta$turbidity.s),1),
                    urchin.s = rep(mean(delta$urchin.s),1))
# convert to model matrix
X_new <- as.matrix(model.matrix(~ X_new$dhw.s + X_new$macro.s + X_new$parrot.s + X_new$ssta_freq_sd.s + 
                                  X_new$temp_max.s + X_new$depth.s + X_new$exposure_fac +  X_new$cyclone_mag.s + 
                                  X_new$turbidity.s + X_new$urchin.s + X_new$dhw.s*X_new$macro.s + X_new$dhw.s*X_new$parrot.s + 
                                  X_new$dhw.s*X_new$ssta_freq_sd.s + X_new$dhw.s*X_new$temp_max.s + X_new$dhw.s*X_new$depth.s + 
                                  X_new$dhw.s*X_new$exposure_fac + X_new$dhw.s*X_new$cyclone_mag.s + X_new$dhw.s*X_new$turbidity.s + 
                                  X_new$dhw.s*X_new$urchin.s))
X_new <- as.matrix(X_new[,2:ncol(X_new)]) # remove intercept
mcmc_go <- mu_b0 + beta_matrix %*% X_new # predict from model output

X_new <- data.frame(dhw.s = new_values$dhw.s[2],
                    macro.s = new_values$macro.s[2],
                    parrot.s = rep(mean(delta$parrot.s),1),
                    ssta_freq_sd.s = rep(mean(delta$ssta_freq_sd.s),1),
                    temp_max.s = rep(mean(delta$temp_max.s),1),
                    depth.s = rep(mean(delta$depth.s),1),
                    exposure_fac = rep(0,1),
                    cyclone_mag.s = rep(mean(delta$cyclone_mag.s),1),
                    turbidity.s = rep(mean(delta$turbidity.s),1),
                    urchin.s = rep(mean(delta$urchin.s),1))
# convert to model matrix
X_new <- as.matrix(model.matrix(~ X_new$dhw.s + X_new$macro.s + X_new$parrot.s + X_new$ssta_freq_sd.s + 
                                  X_new$temp_max.s + X_new$depth.s + X_new$exposure_fac +  X_new$cyclone_mag.s + 
                                  X_new$turbidity.s + X_new$urchin.s + X_new$dhw.s*X_new$macro.s + X_new$dhw.s*X_new$parrot.s + 
                                  X_new$dhw.s*X_new$ssta_freq_sd.s + X_new$dhw.s*X_new$temp_max.s + X_new$dhw.s*X_new$depth.s + 
                                  X_new$dhw.s*X_new$exposure_fac + X_new$dhw.s*X_new$cyclone_mag.s + X_new$dhw.s*X_new$turbidity.s + 
                                  X_new$dhw.s*X_new$urchin.s))
X_new <- as.matrix(X_new[,2:ncol(X_new)]) # remove intercept
mcmc_inc <- mu_b0 + beta_matrix %*% X_new # predict from model output

# calculate 'improvement' and calculate intervals across posteriors
median(mcmc_inc - mcmc_go)
predict(delta_lm, newdata=data.frame(B_2=median(mcmc_inc - mcmc_go)))
predict(delta_lm, newdata=data.frame(B_2=quantile(mcmc_inc - mcmc_go,0.1)))
predict(delta_lm, newdata=data.frame(B_2=quantile(mcmc_inc - mcmc_go,0.9)))


X_new <- data.frame(dhw.s = new_values$dhw.s[3],
                    macro.s = new_values$macro.s[3],
                    parrot.s = rep(mean(delta$parrot.s),1),
                    ssta_freq_sd.s = rep(mean(delta$ssta_freq_sd.s),1),
                    temp_max.s = rep(mean(delta$temp_max.s),1),
                    depth.s = rep(mean(delta$depth.s),1),
                    exposure_fac = rep(0,1),
                    cyclone_mag.s = rep(mean(delta$cyclone_mag.s),1),
                    turbidity.s = rep(mean(delta$turbidity.s),1),
                    urchin.s = rep(mean(delta$urchin.s),1))
# convert to model matrix
X_new <- as.matrix(model.matrix(~ X_new$dhw.s + X_new$macro.s + X_new$parrot.s + X_new$ssta_freq_sd.s + 
                                  X_new$temp_max.s + X_new$depth.s + X_new$exposure_fac +  X_new$cyclone_mag.s + 
                                  X_new$turbidity.s + X_new$urchin.s + X_new$dhw.s*X_new$macro.s + X_new$dhw.s*X_new$parrot.s + 
                                  X_new$dhw.s*X_new$ssta_freq_sd.s + X_new$dhw.s*X_new$temp_max.s + X_new$dhw.s*X_new$depth.s + 
                                  X_new$dhw.s*X_new$exposure_fac + X_new$dhw.s*X_new$cyclone_mag.s + X_new$dhw.s*X_new$turbidity.s + 
                                  X_new$dhw.s*X_new$urchin.s))
X_new <- as.matrix(X_new[,2:ncol(X_new)]) # remove intercept
mcmc_inc <- mu_b0 + beta_matrix %*% X_new # predict from model output

# calculate 'improvement' and calculate intervals across posteriors
median(mcmc_inc - mcmc_go)
predict(delta_lm, newdata=data.frame(B_2=median(mcmc_inc - mcmc_go)))
predict(delta_lm, newdata=data.frame(B_2=quantile(mcmc_inc - mcmc_go,0.1)))
predict(delta_lm, newdata=data.frame(B_2=quantile(mcmc_inc - mcmc_go,0.9)))
