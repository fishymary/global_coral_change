
#   -----------------------------------------------------------------------
# Macroalgae, Urchins, and Heat Stress Magnify Coral Loss After Coral Bleaching
# Code and analysis by Mary Donovan
# 1_run_model
#   -----------------------------------------------------------------------

# This code reads in global reef monitoring data by the Reef Check program, formats the data for input into the model, runs the model in JAGS, and exports model outputs. The model is a Bayesian hierarchical model with coral cover as a response variable and calculates change in coral cover as a latent variable that is a function of a series of abiotic and biotic predictors. 

# This code was written and executed with R version 3.6.1. It is set up to be executed on a High Performance Computer and run across 3 cores. 

# Initialization ----------------------------------------------------------
library(dplyr) # v.0.8.3
library(rjags) # v.4.9 linked to JAGS 4.3.0
library(parallel) # v.3.6.1
library(lubridate) # v.1.7.4

# Import and format data --------------------------------------------------
# read in data
reef_check <- read.csv('data/reef_check_all.csv')
reef_check$Date <- ymd(reef_check$Date)
preds <- read.csv('data/predictors_reef_check_timeseries.csv')
preds$Date <- ymd(preds$Date)

# set up data for the model where for each site and each focal year is categorized 'before' and the year following 'after'
change_1yr <- data.frame(new_id=NA,grp=NA,coral=NA,Year=NA,Date=lubridate::ymd(NA)) # create an empty data frame to fill
for(i in unique(reef_check$new_id)){
  temp <- reef_check[reef_check$new_id==i,]
  temp <- temp[!is.na(temp$coral),]
  nyr <- temp %>% group_by(Year) %>% tally()
  
  for(k in unique(temp$Year)){
    if(nrow(nyr[c(nyr$Year == k),]) >= 1 & nrow(nyr[c(nyr$Year == k+1),]) >= 1){
      before <- temp[c(temp$Year==k),]
      after <- temp[c(temp$Year==k+1),]
      date <- max(temp$Date[c(temp$Year==k)])
      change_1yr <- rbind(change_1yr, data.frame(new_id=i,grp='before',coral=before$coral,Year=k,Date=date))
      change_1yr <- rbind(change_1yr, data.frame(new_id=i,grp='after',coral=after$coral,Year=k,Date=date))
    }
  }
}
change_1yr <- change_1yr[!is.na(change_1yr$new_id),] # remove first empty row
length(unique(change_1yr$new_id)) # check number of unique sites

# join back in bleaching data and subset for bleaching years only
bleaching <- reef_check %>% group_by(new_id, Date) %>% summarise(bleaching=mean(bleaching,na.rm=T)) %>% ungroup()
change_1yr_complete <- change_1yr %>% left_join(bleaching, by=c('new_id','Date')) %>% na.omit() %>% ungroup()
change_1yr_complete <- change_1yr_complete %>% filter(bleaching > 0) %>% ungroup()

# format predictors of delta - join and remove data with missing predictors
change_1yr_complete <- left_join(change_1yr_complete,preds,by=c('new_id','Date','Year')) %>% ungroup()
change_1yr_complete <- change_1yr_complete %>% na.omit() %>% ungroup()

# remove points with no coral and positive bleaching
change_1yr_complete <- change_1yr_complete %>% filter(!(bleaching > 0 & coral==0))

# add a location column
locs <- reef_check %>% distinct(new_id,Location)
change_1yr_complete <- left_join(change_1yr_complete,locs,by='new_id')

# remove lesser antilles (see supplemental material for analyses showing high leverage when included)
change_1yr_complete <- change_1yr_complete %>% filter(!Location=='Lesser Antilles')

# create grouping variable for 'before' and 'after'
change_1yr_complete$grp_fact <- ifelse(change_1yr_complete$grp=='before',0,1)
X_coral <- model.matrix(~ as.factor(change_1yr_complete[,c('grp_fact')]))

# create id vector for 'replicates'
temp <- change_1yr_complete %>% distinct(new_id,Date)
temp$id_fac <- seq(1:nrow(temp))
change_1yr_complete <- left_join(change_1yr_complete,temp,by=c('new_id','Date'))

# create response variable vector. Note beta defined between zero and 1 so need to 'nudge' 0 and 1 data
y <- change_1yr_complete$coral
y <- ifelse(y == 0, 0.001, y)
y <- ifelse(y == 1, 0.999, y)

# create location vector
loc_go <- change_1yr_complete %>% distinct(new_id,Date,id_fac,Location)
loc_go$location <- as.numeric(as.factor(as.character(loc_go$Location)))

# create predictor matrix
X_delta <- change_1yr_complete %>% distinct(new_id,Date,id_fac,dhw_max,parrot,ssta_freq_sd,temp_max,depth,exposure_fac,cyclone_mag,macro,turbidity,urchin_sum_all)

X_delta$dhw.s <- scale(X_delta$dhw_max)[,1]
X_delta$macro.s <- scale(X_delta$macro)[,1]
X_delta$parrot.s <- scale(X_delta$parrot)[,1]
X_delta$ssta_freq_sd.s <- scale(X_delta$ssta_freq_sd)[,1]
X_delta$depth.s <- scale(X_delta$depth)[,1]
X_delta$temp_max.s <- scale(X_delta$temp_max)[,1]
X_delta$cyclone_mag.s <- scale(X_delta$cyclone_mag)[,1]
X_delta$turbidity.s <- scale(X_delta$turbidity)[,1]
X_delta$urchin.s <- scale(log(X_delta$urchin_sum_all+1))[,1]
X_save <- X_delta

X_delta <- as.matrix(model.matrix(~ X_delta$dhw.s + 
                                    X_delta$macro.s +
                                    X_delta$parrot.s +
                                    X_delta$ssta_freq_sd.s + 
                                    X_delta$temp_max.s + 
                                    X_delta$depth.s + 
                                    X_delta$exposure_fac + 
                                    X_delta$cyclone_mag.s +
                                    X_delta$turbidity.s +
                                    X_delta$urchin.s +
                                    X_delta$dhw.s*X_delta$macro.s +
                                    X_delta$dhw.s*X_delta$parrot.s +
                                    X_delta$dhw.s*X_delta$ssta_freq_sd.s +
                                    X_delta$dhw.s*X_delta$temp_max.s +
                                    X_delta$dhw.s*X_delta$depth.s +
                                    X_delta$dhw.s*X_delta$exposure_fac +
                                    X_delta$dhw.s*X_delta$cyclone_mag.s +
                                    X_delta$dhw.s*X_delta$turbidity.s +
                                    X_delta$dhw.s*X_delta$urchin.s 
))
X_delta <- as.matrix(X_delta[,2:ncol(X_delta)]) # remove intercept

# create predictor matrices for marginal effects plotting
X_new_DHW <- data.frame(dhw.s = seq(from=min(X_save$dhw.s),to=max(X_save$dhw.s),length=100),
                        macro.s = rep(mean(X_save$macro.s),100),
                        parrot.s = rep(mean(X_save$parrot.s),100),
                        ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                        temp_max.s = rep(mean(X_save$temp_max.s),100),
                        depth.s = rep(mean(X_save$depth.s),100),
                        exposure_fac = rep(0,100),
                        cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                        turbidity.s = rep(mean(X_save$turbidity.s),100),
                        urchin.s = rep(mean(X_save$urchin.s),100))

X_new_DHW <- as.matrix(model.matrix(~ X_new_DHW$dhw.s + 
                                      X_new_DHW$macro.s +
                                      X_new_DHW$parrot.s +
                                      X_new_DHW$ssta_freq_sd.s + 
                                      X_new_DHW$temp_max.s + 
                                      X_new_DHW$depth.s + 
                                      X_new_DHW$exposure_fac + 
                                      X_new_DHW$cyclone_mag.s +
                                      X_new_DHW$turbidity.s +
                                      X_new_DHW$urchin.s +
                                      X_new_DHW$dhw.s*X_new_DHW$macro.s +
                                      X_new_DHW$dhw.s*X_new_DHW$parrot.s +
                                      X_new_DHW$dhw.s*X_new_DHW$ssta_freq_sd.s +
                                      X_new_DHW$dhw.s*X_new_DHW$temp_max.s +
                                      X_new_DHW$dhw.s*X_new_DHW$depth.s +
                                      X_new_DHW$dhw.s*X_new_DHW$exposure_fac +
                                      X_new_DHW$dhw.s*X_new_DHW$cyclone_mag.s +
                                      X_new_DHW$dhw.s*X_new_DHW$turbidity.s +
                                      X_new_DHW$dhw.s*X_new_DHW$urchin.s 
))
X_new_DHW <- as.matrix(X_new_DHW[,2:ncol(X_new_DHW)])

X_new_MA_lD <- data.frame(dhw.s = rep(-0.2765587,100),
                          macro.s = seq(from=min(X_save$macro.s),to=max(X_save$macro.s),length=100),
                          parrot.s = rep(mean(X_save$parrot.s),100),
                          ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                          temp_max.s = rep(mean(X_save$temp_max.s),100),
                          depth.s = rep(mean(X_save$depth.s),100),
                          exposure_fac = rep(0,100),
                          cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                          turbidity.s = rep(mean(X_save$turbidity.s),100),
                          urchin.s = rep(mean(X_save$urchin.s),100))

X_new_MA_lD <- as.matrix(model.matrix(~ X_new_MA_lD$dhw.s + 
                                        X_new_MA_lD$macro.s +
                                        X_new_MA_lD$parrot.s +
                                        X_new_MA_lD$ssta_freq_sd.s + 
                                        X_new_MA_lD$temp_max.s + 
                                        X_new_MA_lD$depth.s + 
                                        X_new_MA_lD$exposure_fac + 
                                        X_new_MA_lD$cyclone_mag.s +
                                        X_new_MA_lD$turbidity.s +
                                        X_new_MA_lD$urchin.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$macro.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$parrot.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$ssta_freq_sd.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$temp_max.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$depth.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$exposure_fac +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$cyclone_mag.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$turbidity.s +
                                        X_new_MA_lD$dhw.s*X_new_MA_lD$urchin.s 
))
X_new_MA_lD <- as.matrix(X_new_MA_lD[,2:ncol(X_new_MA_lD)])

X_new_MA_hD <- data.frame(dhw.s = rep(2.0934591,100),
                          macro.s = seq(from=min(X_save$macro.s),to=max(X_save$macro.s),length=100),
                          parrot.s = rep(mean(X_save$parrot.s),100),
                          ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                          temp_max.s = rep(mean(X_save$temp_max.s),100),
                          depth.s = rep(mean(X_save$depth.s),100),
                          exposure_fac = rep(0,100),
                          cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                          turbidity.s = rep(mean(X_save$turbidity.s),100),
                          urchin.s = rep(mean(X_save$urchin.s),100))

X_new_MA_hD <- as.matrix(model.matrix(~ X_new_MA_hD$dhw.s + 
                                        X_new_MA_hD$macro.s +
                                        X_new_MA_hD$parrot.s +
                                        X_new_MA_hD$ssta_freq_sd.s + 
                                        X_new_MA_hD$temp_max.s + 
                                        X_new_MA_hD$depth.s + 
                                        X_new_MA_hD$exposure_fac + 
                                        X_new_MA_hD$cyclone_mag.s +
                                        X_new_MA_hD$turbidity.s +
                                        X_new_MA_hD$urchin.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$macro.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$parrot.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$ssta_freq_sd.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$temp_max.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$depth.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$exposure_fac +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$cyclone_mag.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$turbidity.s +
                                        X_new_MA_hD$dhw.s*X_new_MA_hD$urchin.s 
))
X_new_MA_hD <- as.matrix(X_new_MA_hD[,2:ncol(X_new_MA_hD)])

X_new_depth_lD <- data.frame(dhw.s = rep(-0.2765587,100),
                             macro.s = rep(mean(X_save$macro.s),100),
                             parrot.s = rep(mean(X_save$parrot.s),100),
                             ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                             temp_max.s = rep(mean(X_save$temp_max.s),100),
                             depth.s = seq(from=min(X_save$depth.s),to=max(X_save$depth.s),length=100),
                             exposure_fac = rep(0,100),
                             cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                             turbidity.s = rep(mean(X_save$turbidity.s),100),
                             urchin.s = rep(mean(X_save$urchin.s),100))

X_new_depth_lD <- as.matrix(model.matrix(~ X_new_depth_lD$dhw.s + 
                                           X_new_depth_lD$macro.s +
                                           X_new_depth_lD$parrot.s +
                                           X_new_depth_lD$ssta_freq_sd.s + 
                                           X_new_depth_lD$temp_max.s + 
                                           X_new_depth_lD$depth.s + 
                                           X_new_depth_lD$exposure_fac + 
                                           X_new_depth_lD$cyclone_mag.s +
                                           X_new_depth_lD$turbidity.s +
                                           X_new_depth_lD$urchin.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$macro.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$parrot.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$ssta_freq_sd.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$temp_max.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$depth.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$exposure_fac +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$cyclone_mag.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$turbidity.s +
                                           X_new_depth_lD$dhw.s*X_new_depth_lD$urchin.s 
))
X_new_depth_lD <- as.matrix(X_new_depth_lD[,2:ncol(X_new_depth_lD)])

X_new_depth_hD <- data.frame(dhw.s = rep(2.0934591,100),
                             macro.s = rep(mean(X_save$macro.s),100),
                             parrot.s = rep(mean(X_save$parrot.s),100),
                             ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                             temp_max.s = rep(mean(X_save$temp_max.s),100),
                             depth.s = seq(from=min(X_save$depth.s),to=max(X_save$depth.s),length=100),
                             exposure_fac = rep(0,100),
                             cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                             turbidity.s = rep(mean(X_save$turbidity.s),100),
                             urchin.s = rep(mean(X_save$urchin.s),100))

X_new_depth_hD <- as.matrix(model.matrix(~ X_new_depth_hD$dhw.s + 
                                           X_new_depth_hD$macro.s +
                                           X_new_depth_hD$parrot.s +
                                           X_new_depth_hD$ssta_freq_sd.s + 
                                           X_new_depth_hD$temp_max.s + 
                                           X_new_depth_hD$depth.s + 
                                           X_new_depth_hD$exposure_fac + 
                                           X_new_depth_hD$cyclone_mag.s +
                                           X_new_depth_hD$turbidity.s +
                                           X_new_depth_hD$urchin.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$macro.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$parrot.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$ssta_freq_sd.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$temp_max.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$depth.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$exposure_fac +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$cyclone_mag.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$turbidity.s +
                                           X_new_depth_hD$dhw.s*X_new_depth_hD$urchin.s 
))
X_new_depth_hD <- as.matrix(X_new_depth_hD[,2:ncol(X_new_depth_hD)])

X_new_exposure <- data.frame(dhw.s = c(rep(-0.2765587,2),rep(2.0934591,2),rep(mean(X_save$dhw.s),2)),
                             macro.s = rep(mean(X_save$macro.s),6),
                             parrot.s = rep(mean(X_save$parrot.s),6),
                             ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),6),
                             temp_max.s = rep(mean(X_save$temp_max.s),6),
                             depth.s = rep(mean(X_save$depth.s),6),
                             exposure_fac = c(rep(0,1),rep(1,1),rep(0,1),rep(1,1),rep(0,1),rep(1,1)),
                             cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),6),
                             turbidity.s = rep(mean(X_save$turbidity.s),6),
                             urchin.s = rep(mean(X_save$urchin.s),6))

X_new_exposure <- as.matrix(model.matrix(~ X_new_exposure$dhw.s + 
                                           X_new_exposure$macro.s +
                                           X_new_exposure$parrot.s +
                                           X_new_exposure$ssta_freq_sd.s + 
                                           X_new_exposure$temp_max.s + 
                                           X_new_exposure$depth.s + 
                                           X_new_exposure$exposure_fac + 
                                           X_new_exposure$cyclone_mag.s +
                                           X_new_exposure$turbidity.s +
                                           X_new_exposure$urchin.s +
                                           X_new_exposure$dhw.s*X_new_exposure$macro.s +
                                           X_new_exposure$dhw.s*X_new_exposure$parrot.s +
                                           X_new_exposure$dhw.s*X_new_exposure$ssta_freq_sd.s +
                                           X_new_exposure$dhw.s*X_new_exposure$temp_max.s +
                                           X_new_exposure$dhw.s*X_new_exposure$depth.s +
                                           X_new_exposure$dhw.s*X_new_exposure$exposure_fac +
                                           X_new_exposure$dhw.s*X_new_exposure$cyclone_mag.s +
                                           X_new_exposure$dhw.s*X_new_exposure$turbidity.s +
                                           X_new_exposure$dhw.s*X_new_exposure$urchin.s 
))
X_new_exposure <- as.matrix(X_new_exposure[,2:ncol(X_new_exposure)])

X_new_urchin <- data.frame(dhw.s = rep(mean(X_save$dhw.s),100),
                           macro.s = rep(mean(X_save$macro.s),100),
                           parrot.s = rep(mean(X_save$parrot.s),100),
                           ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                           temp_max.s = rep(mean(X_save$temp_max.s),100),
                           depth.s = rep(mean(X_save$depth.s),100),
                           exposure_fac = rep(0,100),
                           cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                           turbidity.s = rep(mean(X_save$turbidity.s),100),
                           urchin.s = seq(from=min(X_save$urchin.s),to=max(X_save$urchin.s),length=100))

X_new_urchin <- as.matrix(model.matrix(~ X_new_urchin$dhw.s + 
                                         X_new_urchin$macro.s +
                                         X_new_urchin$parrot.s +
                                         X_new_urchin$ssta_freq_sd.s + 
                                         X_new_urchin$temp_max.s + 
                                         X_new_urchin$depth.s + 
                                         X_new_urchin$exposure_fac + 
                                         X_new_urchin$cyclone_mag.s +
                                         X_new_urchin$turbidity.s +
                                         X_new_urchin$urchin.s +
                                         X_new_urchin$dhw.s*X_new_urchin$macro.s +
                                         X_new_urchin$dhw.s*X_new_urchin$parrot.s +
                                         X_new_urchin$dhw.s*X_new_urchin$ssta_freq_sd.s +
                                         X_new_urchin$dhw.s*X_new_urchin$temp_max.s +
                                         X_new_urchin$dhw.s*X_new_urchin$depth.s +
                                         X_new_urchin$dhw.s*X_new_urchin$exposure_fac +
                                         X_new_urchin$dhw.s*X_new_urchin$cyclone_mag.s +
                                         X_new_urchin$dhw.s*X_new_urchin$turbidity.s +
                                         X_new_urchin$dhw.s*X_new_urchin$urchin.s 
))
X_new_urchin <- as.matrix(X_new_urchin[,2:ncol(X_new_urchin)])

X_new_MA <- data.frame(dhw.s = rep(mean(X_save$dhw.s),100),
                       macro.s = seq(from=min(X_save$macro.s),to=max(X_save$macro.s),length=100),
                       parrot.s = rep(mean(X_save$parrot.s),100),
                       ssta_freq_sd.s = rep(mean(X_save$ssta_freq_sd.s),100),
                       temp_max.s = rep(mean(X_save$temp_max.s),100),
                       depth.s = rep(mean(X_save$depth.s),100),
                       exposure_fac = rep(0,100),
                       cyclone_mag.s = rep(mean(X_save$cyclone_mag.s),100),
                       turbidity.s = rep(mean(X_save$turbidity.s),100),
                       urchin.s = rep(mean(X_save$urchin.s),100))

X_new_MA <- as.matrix(model.matrix(~ X_new_MA$dhw.s + 
                                     X_new_MA$macro.s +
                                     X_new_MA$parrot.s +
                                     X_new_MA$ssta_freq_sd.s + 
                                     X_new_MA$temp_max.s + 
                                     X_new_MA$depth.s + 
                                     X_new_MA$exposure_fac + 
                                     X_new_MA$cyclone_mag.s +
                                     X_new_MA$turbidity.s +
                                     X_new_MA$urchin.s +
                                     X_new_MA$dhw.s*X_new_MA$macro.s +
                                     X_new_MA$dhw.s*X_new_MA$parrot.s +
                                     X_new_MA$dhw.s*X_new_MA$ssta_freq_sd.s +
                                     X_new_MA$dhw.s*X_new_MA$temp_max.s +
                                     X_new_MA$dhw.s*X_new_MA$depth.s +
                                     X_new_MA$dhw.s*X_new_MA$exposure_fac +
                                     X_new_MA$dhw.s*X_new_MA$cyclone_mag.s +
                                     X_new_MA$dhw.s*X_new_MA$turbidity.s +
                                     X_new_MA$dhw.s*X_new_MA$urchin.s 
))
X_new_MA <- as.matrix(X_new_MA[,2:ncol(X_new_MA)])

# calculate sample size for model loop
n_rep <- length(unique(change_1yr_complete$id_fac))

# save data to use in map
coord_export <- reef_check %>% distinct(new_id,Lat,Long,Location)
coord_export <- coord_export[coord_export$new_id %in% change_1yr_complete$new_id,]
write.csv(coord_export,'outputs/coord_export.csv',row.names=F)

# define model input data
data.in <- list(y=y,
                N=nrow(change_1yr_complete),
                replicate=change_1yr_complete$id_fac,
                X_coral=X_coral,
                J=n_rep,
                X_delta=X_delta,
                C=ncol(X_delta),
                location=loc_go$location,
                K=length(unique(loc_go$location)),
                X_new_DHW=X_new_DHW,
                X_new_MA_lD=X_new_MA_lD,X_new_MA_hD=X_new_MA_hD,X_new_MA=X_new_MA,
                X_new_depth_lD=X_new_depth_lD,X_new_depth_hD=X_new_depth_hD,
                X_new_exposure=X_new_exposure,
                X_new_urchin=X_new_urchin,
                prior.scale=10
)

# Define model ------------------------------------------------------------
cat("model{
    # likelihood
    for(i in 1:N){
      y[i] ~ dbeta(r*mu[i],r*(1-mu[i]))
      logit(mu[i]) <- a[replicate[i]] + inprod(b[replicate[i]],X_coral[i,2])
      y.new[i] ~ dbeta(r*mu[i],r*(1-mu[i]))
    }
    r ~ dgamma(0.1,0.1)
    
    for (j in 1:J){
      a[j] <- B[j,1]
      b[j] <- B[j,2]
      B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,]) 
      B.hat[j,1] <- mu.a    
      B.hat[j,2] <- b0[location[j]] + inprod(beta_d,X_delta[j,])
      
      delta[j] <- ilogit(a[j] + b[j]) - ilogit(a[j])
    }
    mu.a ~ dnorm(0,0.001)

    for(k in 1:K){
      b0[k] ~ dnorm(mu_b0,tau_b0)
    }
    tau_b0 ~ dgamma(0.5,0.5)
    sigma_b0 <- abs(xi)/sqrt(tau_b0)
    xi ~ dnorm(0, tau_xi)I(0,)
    tau_xi <- pow(prior.scale, -2)
    mu_b0 ~ dnorm(0,0.001)

    for(c in 1:C){
        beta_d[c] ~ dnorm(0,0.001)
    }

    for (k in 1:2) {
      sigma[k] ~ dunif(0,10) 
    }

    Tau.B[1:2,1:2] <- inverse(Sigma.B[,]) 
    Sigma.B[1,1] <- sigma[1]^2
    Sigma.B[2,2] <- sigma[2]^2
    rho ~ dunif(-1,0)
    Sigma.B[1,2] <- rho*sigma[1]*sigma[2]
    Sigma.B[2,1] <- Sigma.B[1,2]
    
    # posterior checks
    y.mean <- mean(y)
    ynew.mean <- mean(y.new)
    pval.mean <- step(ynew.mean-y.mean)
    
    sd.y <- sd(y)
    sd.y.new <- sd(y.new)
    pval.sd <- step(sd.y.new-sd.y)
    
    # R-squared
    varF <- sd(y)^2
    varE <- sd(y - y.new)^2
    R2 <- varF/(varF+varE)
    
    # predictions
    for(z in 1:100){
        pred_dhw[z] <- mu_b0 + inprod(beta_d,X_new_DHW[z,])
        pred_MA[z] <- mu_b0 + inprod(beta_d,X_new_MA[z,])
        pred_MA_lD[z] <- mu_b0 + inprod(beta_d,X_new_MA_lD[z,])
        pred_MA_hD[z] <- mu_b0 + inprod(beta_d,X_new_MA_hD[z,])
        pred_depth_lD[z] <- mu_b0 + inprod(beta_d,X_new_depth_lD[z,])
        pred_depth_hD[z] <- mu_b0 + inprod(beta_d,X_new_depth_hD[z,])
        pred_urchin[z] <- mu_b0 + inprod(beta_d,X_new_urchin[z,])
    }
    for(z in 1:6){
        pred_exposure[z] <- mu_b0 + inprod(beta_d,X_new_exposure[z,])
    }
    
    }",file='change_bayes_wMultX.jags')


# Run model ---------------------------------------------------------------
n.adapt <- 1000; n.update <- 900000; n.iter <- 100000

initfunc <- function(){return(list(
  # B=array(rnorm(2*(n_rep)),c((n_rep),2)),
  r = runif(1,30,40),
  rho = runif(1,-0.8,-0.5),
  mu.a=-1,
  mu_b0=0,
  # sigma_b0=0.1,
  sigma=c(1,0.15),
  beta_d=runif(19,-0.1,0.1)
))}

# run chains in parallel
cl <- makeCluster(3) # this determines the number of chains, must be less than the number of cores on computer
clusterExport(cl, c('data.in','n.adapt','n.update','n.iter','initfunc','n_rep'))

out <- clusterEvalQ(cl,{
  library(rjags)
  m <- jags.model("change_bayes_wMultX.jags", data = data.in, n.chains=1, n.adapt=n.adapt,inits = initfunc())
  update(m, n.iter = n.update) 
  zmCore <- coda.samples(m, variable.names = c('delta','b0','r','rho','B.hat','mu.a','sigma','xi',
                                               'y.new','pval.mean','pval.sd','R2', 'a','b',
                                               'beta_d','B','mu_b0','sigma_b0','beta_d_loc','sigma_beta','sigma',
                                               'pred_dhw','pred_MA','pred_MA_lD','pred_MA_hD','pred_depth_lD','pred_depth_hD'
                                               ,'pred_exposure','pred_urchin'), 
                         n.iter=n.iter, n.thin=10)
  return(as.mcmc(zmCore))
})
zmPrp <- mcmc.list(out)
stopCluster(cl)

saveRDS(zmPrp,file='outputs/change_1yr_out_coda.Rdata')



# export model outputs ----------------------------------------------------
grepgo <- grep('delta',colnames(zmPrp[[1]]))
delta_sum_1y <- summary(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))
delta_out_1y <- data.frame(delta=delta_sum_1y$quantiles[,'50%'],
                           delta_up=delta_sum_1y$quantiles[,'97.5%'],
                           delta_down=delta_sum_1y$quantiles[,'2.5%'])

data_out <- cbind(temp,X_save,delta_out_1y)

write.csv(data_out, 'outputs/change_1yr_out_delta.csv',row.names = F)

write.csv(data.frame(y=y,id_fac=change_1yr_complete$id_fac,change_1yr_complete$coral,change_1yr_complete$new_id,change_1yr_complete$grp,change_1yr_complete$Date), 'outputs/change_1yr_out_inputdata.csv',row.names=F)

write.csv(X_delta, 'outputs/change_1yr_out_Xdata.csv',row.names=F)

write.csv(loc_go %>% distinct(Location,location), 'outputs/change_1yr_out_locs.csv',row.names=F)
write.csv(loc_go, 'outputs/change_1yr_out_locs_all.csv',row.names=F)


