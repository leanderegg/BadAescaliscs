

#+++++++++++++++++++++++++++++++++++++++++++++++++++++
#          STEP 3: LOAD gasex data and combine lwp, aba, gs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++

require(tidyr)
require(dplyr)
require(stringr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(zoo)
library(lme4)
library(lmerTest)

####### START: Load Data ############################
#rm(list=ls())
# set working directory
#setwd("~/Documents/RESEARCH/California Buckeye/")

#se <- function(x,na.rm=T) sqrt(var(x))/sum(!is.na(x))


# set our color palette cause default be hella gnarly
mypal <- brewer.pal(n=8, "Set2")
palette(mypal)


# load raw gas exchnage data
gasexchange<-read.csv("gasexdata_AECA.csv", header = T)


####. Add in date and individual Tag, and treatment info #######
# 4th of july data doesn't have year
gasexchange$Date <- as.character(gasexchange$Date)
gasexchange$Date[which(gasexchange$Date=="4-Jul")] <- "4-Jul-18" # add a date into that shit

gasexchange.av <- gasexchange %>% group_by(Date) %>% summarise(mAmax = mean(Amax,na.rm=T), sdAmax= sd(Amax, na.rm=T))

gasexchange$Date <- as.Date(gasexchange$Date, "%d-%b-%y")
gasexchange <- gasexchange[order(gasexchange$Date),]
gasexchange$tag <- as.numeric(str_split(gasexchange$ID, pattern = "_",simplify = T)[,2])
# add in treatment
gasexchange$treatment <- rep(NA, length(gasexchange[,1]))

control <- c(589, 498, 10954,159)
mild <- c(592, 10999, 600)
severe <- c(590, 588, 591)
gasexchange$treatment[gasexchange$tag %in% control] <- "control"
gasexchange$treatment[gasexchange$tag %in% mild] <- "mild"
gasexchange$treatment[gasexchange$tag %in% severe] <- "severe"



####### Combine ABA data w/ gasex #############

## aba from ABA.R (Step 2)

# quick rename "999" and "954" to "10###"
aba$tag[which(aba$tag == 999)] <- 10999
aba$tag[which(aba$tag == 954)] <- 10954
# average any repeated measurements to the day (11 repeated measurements)
aba.clean <- aba %>% group_by(tag, treatment, date) %>% summarise(ABAFWngg = mean(ABAFWngg, na.rm=T))
# aba$tag <- as.numeric(aba$tag) # make numeric to match gasexchange, which as a numeric 'tag' column
gsaba <- left_join(gasexchange, aba.clean, by=c("tag"="tag","Date"="date", "treatment"="treatment"))
# gsaba$time.since.start <- NULL # kill time.since.start - but now killed in the combination to aba.clean step
gsaba$lwp <- NULL
gsaba$treatment <- factor(gsaba$treatment)
levels(gsaba$treatment) <- list(control="control",mwd="mild",swd="severe")
head(gsaba)




####### Combine LWP data w/ gasex #############

## wps.clean from Aescal_Greenhous_WP_Select.R (Step 1)

# with Rob's old janky ass dataframe
# enchilada <- left_join(gsaba, lwppd.all %>% select(doy.yr, tag, lwppd_MPa), by = c("Date"="doy.yr","tag"="tag"))


# with Lee's new method
# load in most recent version of wps.clean so don't have to run WP_Select code
wps.clean <- read.csv("WP_data_processed20200429.csv")[,-1]
wps.clean$DOY.yr <- as.Date(wps.clean$DOY.yr) # make sure DOY.yr is in Date format


# merge gs & ABA data with wp data
enchilada <- full_join(gsaba, wps.clean %>% select(DOY.yr, tag,treatment, lwp.m, lwp.n, flags), by = c("Date"="DOY.yr","tag"="tag", "treatment"="treatment")) %>% arrange(tag, Date)

enchilada$treatment <- factor(enchilada$treatment)
enchilada$tag <- factor(enchilada$tag)
enchilada$WUE <- enchilada$Amax/enchilada$gs
# kill one obvious outlier
enchilada$Amax[which(enchilada$WUE>0.4)] <- NA
enchilada$gs[which(enchilada$WUE>0.4)] <- NA
enchilada$WUE[which(enchilada$WUE>0.4)] <- NA
# and kill some useless columns that we'll add back in eventually
enchilada$time.since.start <- NULL
enchilada$yr <- NULL
enchilada$DOY <- NULL
# make a column to destinguish when things started getting rewatered (pch=16 -> watered pch=1 -> not watered)
#enchilada$watered <- 16

# make a day of year column
enchilada$DOY <- as.numeric(strftime(enchilada$Date, format = "%j"))
# subract out days since June 28th 
enchilada$time.since.start <- as.numeric(enchilada$DOY) - 179
# and make everything pre-drought ==0
enchilada$time.since.start[which(enchilada$time.since.start<0)] <- 0

# make a time since rewatering column
enchilada$time.since.rewater <- NA
# cacluate rewatering for each tree seperately
enchilada$time.since.rewater[which(enchilada$tag=="591")] <- enchilada$DOY[which(enchilada$tag=="591")] - 198
enchilada$time.since.rewater[which(enchilada$tag=="10999")] <- enchilada$DOY[which(enchilada$tag=="10999")] - 200
enchilada$time.since.rewater[which(enchilada$tag=="592")] <- enchilada$DOY[which(enchilada$tag=="592")] - 200
enchilada$time.since.rewater[which(enchilada$tag=="600")] <- enchilada$DOY[which(enchilada$tag=="600")] - 200
enchilada$time.since.rewater[which(enchilada$tag=="590")] <- enchilada$DOY[which(enchilada$tag=="590")] - 201
enchilada$time.since.rewater[which(enchilada$tag=="588")] <- enchilada$DOY[which(enchilada$tag=="588")] - 202
enchilada$time.since.rewater[which(enchilada$time.since.rewater<0)] <- 0

# make a variable for watering status (0=watered, 1=droughted, 2=rewatered)
enchilada$watered <- 0
enchilada$watered[which(enchilada$tag %in% c("591","10999","592","600","590","588") & enchilada$time.since.start>0)] <- 1
enchilada$watered[which(enchilada$time.since.rewater>0)] <- 2

# also making a factor for the water status
enchilada$water.status <- factor(enchilada$watered)
levels(enchilada$water.status) <- list(watered=0, drought=1, rewatered=2)

# it's getting annoying to clean the lwp by filtering for 'flags==0', so just going to make lwp.clean column
enchilada$lwp.clean <- enchilada$lwp.m
enchilada$lwp.clean[which(enchilada$flags>0)] <- NA

#### Smoothing lwp with a loess smoother:

## quick visual of what this code is doing (for random individual 589):
# lwp.smoothed <- loess(lwp.clean~DOY, data=enchilada[which(enchilada$tag==589),])
# sm <- predict(lwp.smoothed)
# plot(lwp.clean~DOY, enchilada[which(enchilada$tag==589),])
# lines(sm~enchilada$DOY[which(enchilada$tag==589 & enchilada$lwp.clean<1)])

# now loop through individuals (except 159, which has no data), fit the smoother, use predict to add to df
enchilada$lwp.smooth <- NA
for(i in unique(enchilada$tag)[-1]){
  smoothed.lwp <- loess(lwp.clean~DOY, data=enchilada[which(enchilada$tag==i),])
  enchilada$lwp.smooth[which(enchilada$tag==i & enchilada$lwp.clean<1)] <- predict(smoothed.lwp)
}

enchilada$time.period <- NA
enchilada$time.period[which(enchilada$DOY<179)] <- "pre"
enchilada$time.period[which(enchilada$DOY>190 & enchilada$DOY<202 & enchilada$watered<2)] <- "drought"
enchilada$time.period[which(enchilada$DOY>203)] <- "recovery"

 

# ####### Adding some derived metrics ############

# set a minimum water potential, gs, Amax and max ABA for each individual
  # only using period between start of drydown and the start of rewatering
ench <- enchilada %>% group_by(tag, treatment) %>% filter(DOY>179 & watered<2) %>% summarise(minlwp.flag = min(lwp.m, na.rm=T),minlwp=min(lwp.clean, na.rm=T), minlwp.smooth=min(lwp.smooth, na.rm=T), maxABA = max(ABAFWngg, na.rm=T), mings=min(gs, na.rm=T), minAmax = min(Amax, na.rm=T))

# calculate the maximum post-experiment rates for recovery
recov <- enchilada %>% group_by(tag) %>% 
  filter(time.since.rewater>0) %>% 
  summarise(maxlwp.recov = max(lwp.clean, na.rm=T), maxlwp.flag.recov= max(lwp.m, na.rm=T), maxlwp.smooth.recov=max(lwp.smooth, na.rm=T),
            minABA.recov = min(ABAFWngg, na.rm=T), maxAmax.recov = max(Amax, na.rm=T), maxgs.recov=max(gs, na.rm=T))


enc <- left_join(ench, recov)


#### Quick ANOVAs to see whether treatments actually differed:

summary(aov(minlwp~treatment, ench[-1,])) # p=0.241
summary(aov(minlwp.flag~treatment, ench[-1,])) # p=0.0745 .
summary(aov(minlwp.smooth~treatment, ench[-1,])) # p=0.394
summary(aov(maxABA~treatment, ench[-1,])) # p=0.151 (looks obvious that controls are lower max ABA, but low sample size) 
summary(aov(mings~treatment, ench[-1,])) # p=0.00173**
summary(aov(minAmax~treatment, ench[-1,])) # p=0.0253*
  # so gas ex does, ABA is close, but filtered lwp doesn't and unfiltered is close


# pull out the date of the lowest/highest values
ench$date.min.lwp <- as.Date(NA)
ench$date.max.ABA <- as.Date(NA)
ench$date.min.gs <- as.Date(NA)
ench$date.min.Amax <- as.Date(NA)
for(i in unique(enchilada$tag)[-1]){
  tmp <- enchilada[which(enchilada$tag==i & enchilada$DOY>179 & enchilada$watered<2),]
  ench$date.min.lwp[which(ench$tag==i)] <- as.Date(tmp$Date[which(tmp$lwp.m == min(tmp$lwp.m, na.rm=T))])
  ench$date.max.ABA[which(ench$tag==i)] <- as.Date(tmp$Date[which(tmp$ABAFWngg == max(tmp$ABAFWngg, na.rm=T))])
  ench$date.min.gs[which(ench$tag==i)] <- as.Date(tmp$Date[which(tmp$gs == min(tmp$gs, na.rm=T))])
  ench$date.min.Amax[which(ench$tag==i)] <- as.Date(tmp$Date[which(tmp$Amax == min(tmp$Amax, na.rm=T))])
  
}

# add in to enchilada a column indicating min wp, and date of min so we can analyze recovery
enchilada$minlwp <- ench$minlwp[match(enchilada$tag, ench$tag)]
enchilada$date.minlwp <- ench$date.min.lwp[match(enchilada$tag, ench$tag)]
enchilada$minlwp.smooth <- ench$minlwp.smooth[match(enchilada$tag, ench$tag)]
enchilada$maxABA <- ench$maxABA[match(enchilada$tag, ench$tag)]
enchilada$date.maxABA <- ench$date.max.ABA[match(enchilada$tag, ench$tag)]
enchilada$mings <- ench$mings[match(enchilada$tag, ench$tag)]
enchilada$date.mings <- ench$date.min.gs[match(enchilada$tag, ench$tag)]
enchilada$minAmax <- ench$minAmax[match(enchilada$tag, ench$tag)]
enchilada$date.minAmax <- ench$date.min.Amax[match(enchilada$tag, ench$tag)]

enchilada$maxlwp <- enc$maxlwp.recov[match(enchilada$tag, enc$tag)]
enchilada$maxlwp.smooth <- enc$maxlwp.smooth.recov[match(enchilada$tag, enc$tag)]
enchilada$minABA <- enc$minABA.recov[match(enchilada$tag, enc$tag)]
enchilada$maxgs <- enc$maxgs.recov[match(enchilada$tag, enc$tag)]
enchilada$maxAmax <- enc$maxAmax.recov[match(enchilada$tag, enc$tag)]




enchilada$Amax[which(enchilada$tag==590 & enchilada$time.since.rewater>0)] <- 0.4
enchilada$gs[which(enchilada$tag==590 & enchilada$time.since.rewater>0)] <- 27

### calculate treatment means ###
treat.means <- enchilada %>% group_by(DOY, Date, treatment) %>% summarise(lwp=mean(lwp.m, na.rm=T), ABA = mean(ABAFWngg, na.rm=T), gs = mean(gs, na.rm=T), E=mean(E, na.rm=T), Amax=mean(Amax, na.rm=T))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########## BEGIN: FIGURES #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#_________________________________________________________
####### TIMECOURSE FIGURES #########
#_________________________________________________________


#### .. Old line figures #####
quartz(width=5, height=6)
par(mfrow=c(3,1), mar=c(0,4,0,1), oma=c(4,0,1,0))
plot(gs~Date, enchilada, col=treatment, pch=16, type="n")
for(i in unique(enchilada$tag)){
  lines(gs~Date, enchilada[which(enchilada$tag==i),], col=treatment, lwd=2)
}
abline(v=as.Date("2018-06-28", "%F"))
plot(lwp.m~Date, enchilada, col=treatment, pch=16, type="n")
for(i in unique(enchilada$tag)){
  lines(lwp.m~Date, enchilada[which(enchilada$tag==i & enchilada$lwp.m<0),], col=treatment, lwd=2)
}
abline(v=as.Date("2018-06-28", "%F"))
plot(ABAFWngg~Date, enchilada, col=treatment, pch=16, type="n")
for(i in unique(enchilada$tag)){
  lines(ABAFWngg~Date, enchilada[which(enchilada$tag==i & enchilada$ABAFWngg>0),], col=treatment, lwd=2)
}
abline(v=as.Date("2018-06-28", "%F"))

#quartz.save("gs_lwp_ABA_throughtime_lines_v1.pdf")


### .. Newer point figures (all inds one plot) ####

# *** ALL DATA ***
  # note: to add a line for the treatment means, add: 
  #           + geom_line(data=treat.means[which(!is.na(treat.means$gs)),], aes(col=treatment, shape=NULL), size=1.5)
  #       to add a line for the loess smothing of each treatment, add:   
  #           + geom_smooth(se=F, aes(pch=NULL))
# #*** Gs
# ggplot(enchilada, aes(x=DOY, y=gs, col=factor(treatment), shape=water.status)) + geom_point() + geom_vline(xintercept=179) +
#   geom_line(data=treat.means[which(!is.na(treat.means$gs)),], aes(col=treatment, shape=NULL), size=1.5)
# 
# #*** E
# ggplot(enchilada, aes(x=DOY, y=E, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)
# #*** Amax
# ggplot(enchilada, aes(x=DOY, y=Amax, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)
# #*** lwp (without removing flags)
# # ggplot(enchilada, aes(x=DOY, y=lwp.m, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)
# #*** lwp (flags removed)
# ggplot(enchilada, aes(x=DOY, y=lwp.clean, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)
# #*** ABA
# ggplot(enchilada, aes(x=DOY, y=gs, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)

## make special dataset for statistics, which only has the pre, drought, and recovery periods, no in betweens
ench.for.stats <- enchilada[which(!is.na(enchilada$time.period)),]

# *** Experiment Only ***
    # removed data pre-DOY=170
#*** Gs
tmp <- enchilada[-which(enchilada$tag==590 & enchilada$time.since.rewater>0), ]
ggplot(tmp[which(enchilada$DOY>170),], aes(x=DOY, y=gs, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179) + geom_smooth(se=F, aes(pch=NULL)) + geom_point(aes(x=201, y=27, col=NULL))
# gsmod <- lmer(gs~treatment*time.period + (1|tag), ench.for.stats )
# summary(gsmod)
# 
# gsaov <- aov(gs~treatment*time.period, ench.for.stats)
# gshsd <- TukeyHSD(gsaov)

#*** E
ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=E, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))+ geom_point(aes(x=201, y=.084, col=NULL))
#*** Amax
ggplot(tmp[which(tmp$DOY>170),], aes(x=DOY, y=Amax, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))+ geom_point(aes(x=201, y=.4, col=NULL))
#*** lwp (without removing flags)
# ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=lwp.m, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))
#*** lwp (flags removed)
ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=lwp.clean, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))
# ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=lwp.smooth, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))
#*** ABA
ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=ABAFWngg, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))
#*** WUE
ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=WUE, col=treatment, pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))


# *** Digging in to lwp ***

# Each individual with loess smoother (visualizes my smoothing method)
# ggplot(enchilada, aes(x=DOY, y=lwp.clean, col=treatment)) + geom_point()+geom_smooth() + facet_wrap(facets = ~tag)
# 
# # all individuals on one plot, smoother color=treatment
#   # uncleaned data (6s left in)
# # ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=lwp.m, col=treatment, shape=tag)) + geom_point()+geom_smooth(se=F) + geom_vline(xintercept=200)
#   # cleaned data (no 6s)
# # with ggplot smoother
# ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=lwp.clean, col=treatment, shape=tag)) + geom_point()+geom_smooth(se=F) + geom_vline(xintercept=200)
# # just with smoothed data from df
# ggplot(enchilada[which(enchilada$lwp.clean<1),], aes(x=DOY, y=lwp.clean, col=treatment, shape=tag)) + geom_point() + geom_vline(xintercept=200) + geom_line(aes(y=lwp.smooth))
# # smoothed by treatment (from smoothed lwp)
# ggplot(enchilada[which(enchilada$DOY>170),], aes(x=DOY, y=lwp.smooth, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)  + geom_smooth(se=F, aes(pch=NULL))


#_________________________________________________________
####### PhysVariable SCATTERPLOT FIGURES #########
#_________________________________________________________


### ABA vs GS (with pre-treatment excluded)
ggplot(enchilada[which(gsaba$Date>"2018-06-28"),], aes(x=ABAFWngg, y=gs, col=treatment, shape=water.status)) + geom_point() #+ geom_smooth(se=F,method = "loess",span=1)
ggplot(enchilada[which(gsaba$Date>"2018-06-28"),], aes(x=log(ABAFWngg), y=log(gs), col=treatment, shape=water.status)) + geom_point() #+ geom_smooth(se=F,method = "loess",span=1)

# do a quick general additative model showing the significance of this nonlinear trend:
gs <- enchilada$gs[which(enchilada$ABAFWngg>0 & enchilada$gs>0)]
aba <- enchilada$ABAFWngg[which(enchilada$ABAFWngg>0 & enchilada$gs>0)]
treat <- enchilada$treatment[which(enchilada$ABAFWngg>0 & enchilada$gs>0)]
waterstat <- enchilada$water.status[which(enchilada$ABAFWngg>0 & enchilada$gs>0)]
gsaba.gam <- gam(I(log(gs))~s(I(log(aba)), fx=FALSE, k=-1, bs="cr")) #
summary(gsaba.gam) # very significant (p<2e-16)
# testing for effect of treatment, ns as main effect (on intercept)
gsaba.gam2 <- gam(I(log(gs))~s(I(log(aba)), fx=FALSE, k=-1, bs="cr")+treat)
# testing for different relationships by treatment
gsaba.gam3 <- gam(I(log(gs))~s(I(log(aba)), fx=FALSE, k=-1, bs="cr", by=as.numeric(treat=="control")) +
                    s(I(log(aba)), fx=FALSE, k=-1, bs="cr", by=as.numeric(treat=="swd")) + 
                    s(I(log(aba)), fx=FALSE, k=-1, bs="cr", by=as.numeric(treat=="mwd")))
# all treats have significant smoothers, but are they needed?
AIC(gsaba.gam, gsaba.gam2, gsaba.gam3) # NOPE! one smoother == golden

### Plot the relationship
gsaba.preds <- predict(gsaba.gam, se=T, type="response")
I1 <- order(aba) # order points by ABA, so plotting doesn't go crazy

# plot(log(gs)~log(ABAFWngg), enchilada, col=treatment, pch=as.numeric(enchilada$watered[which(enchilada$gs>0&enchilada$ABAFWngg>0)])+16)
# lines(gsaba.preds$fit[I1]~log(aba[I1]))
# mtext("p<0.001***",side=3, line=-1, adj=.8)

### or plot with ggplot
ggplot() + geom_point(data=enchilada, aes(x=log(ABAFWngg), y=log(gs), col=treatment, pch=water.status)) + geom_line(data=data.frame(fit = gsaba.preds$fit[I1],ABA=aba[I1] ), aes(x=log(ABA), y=fit))




# Gs vs LWP and ABA vs LWP
ggplot(enchilada[which(enchilada$ABAFWngg>0),], aes(x=lwp.smooth, y=gs, col=treatment, pch=water.status)) + geom_point() 
gslwp.gam <- gam(gs~s(lwp.smooth, fx=FALSE, k=-1, bs="cr"), data=enchilada) #
gslwp.gam2 <- gam(gs~s(lwp.smooth, fx=FALSE, k=-1, bs="cr")+ treatment, data=enchilada) #
gslwp.gam3 <- gam(gs~s(lwp.smooth, fx=FALSE, k=-1, bs="cr", by=as.numeric(treatment=="control")) +
                    s(lwp.smooth, fx=FALSE, k=-1, bs="cr", by=as.numeric(treatment=="swd")) + 
                    s(lwp.smooth, fx=FALSE, k=-1, bs="cr", by=as.numeric(treatment=="mwd")), data=enchilada) #
summary(gslwp.gam3) # somehow super significant, but must be meaningless?
AIC(gslwp.gam, gslwp.gam2, gslwp.gam3)

ggplot(enchilada[which(enchilada$ABAFWngg>0 & enchilada$Date>"2018-06-28" ),], aes(y=ABAFWngg, x=lwp.smooth, col=treatment, shape=water.status)) + geom_point()
ABAlwp.gam <- gam(ABAFWngg~s(lwp.smooth, fx=FALSE, k=-1, bs="cr"), data=enchilada) #
summary(ABAlwp.gam) # not significant


# Gs and ABA vs LWP.smooth
ggplot(enchilada[which(enchilada$ABAFWngg>0),], aes(x=lwp.smooth, y=log(gs), col=treatment)) + geom_point() 
ggplot(enchilada[which(enchilada$ABAFWngg>0 & enchilada$Date>"2018-06-28" ),], aes(x=ABAFWngg, y=lwp.smooth, col=treatment, shape=water.status)) + geom_point()
  # ABA actually starts to look ok against the smoothed wp. Gs still looks like shit though
ABAFWngglwp.gam <- gam(ABAFWngg~s(lwp.smooth, fx=FALSE, k=-1, bs="cr"), data=enchilada) #
ABAFWngglwp.gam2 <- gam(ABAFWngg~s(lwp.smooth, fx=FALSE, k=-1, bs="cr")+ treatment, data=enchilada) #
ABAFWngglwp.gam3 <- gam(ABAFWngg~s(lwp.smooth, fx=FALSE, k=-1, bs="cr", by=as.numeric(treatment=="control")) +
                    s(lwp.smooth, fx=FALSE, k=-1, bs="cr", by=as.numeric(treatment=="swd")) + 
                    s(lwp.smooth, fx=FALSE, k=-1, bs="cr", by=as.numeric(treatment=="mwd")), data=enchilada) #
summary(ABAFWngglwp.gam3) # somehow super significant, but must be meaningless?
AIC(ABAFWngglwp.gam, ABAFWngglwp.gam2, ABAFWngglwp.gam3)

# GS 




########## Predictors of recovery ###########

## Gs recovery as f(min lwp)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minlwp.smooth, y=gs, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minlwp.smooth, y=maxgs, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## Amax recovery as f(min lwp)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minlwp.smooth, y=Amax, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minlwp.smooth, y=maxAmax, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## ABA recovery as f(min lwp)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minlwp.smooth, y=ABAFWngg, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minlwp.smooth, y=minABA, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## Gs recovery as f(min Amax)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minAmax, y=gs, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minAmax, y=maxgs, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## Amax recovery as f(min gs)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=mings, y=Amax, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=mings, y=maxAmax, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## ABA recovery as f(min Amax)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minAmax, y=ABAFWngg, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=minAmax, y=minABA, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## Gs recovery as f(max ABA)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=maxABA, y=gs, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=maxABA, y=maxgs, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))


## Amax recovery as f(max ABA)
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=maxABA, y=Amax, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))
ggplot(enchilada[which(enchilada$water.status=="rewatered"),], aes(x=maxABA, y=maxAmax, col=treatment))+ geom_point() + geom_point(aes(x=-1.34, y=0.4, col=NULL))






### Old plots looking at gas echange #######
ggplot(enchilada, aes(x=gs, y=Amax, col=water.status, pch=water.status)) + geom_point()+facet_wrap(facets=~tag)


ggplot(enchilada, aes(x=gs, y=WUE, col=water.status, pch=water.status)) + geom_point()+facet_wrap(facets=~tag)

ggplot(enchilada[which(enchilada$flags==0),], aes(x=Date, y=lwp.m, col=water.status, pch=water.status)) + geom_point()+facet_wrap(facets=~tag)

ggplot(enchilada[which(enchilada$flags==0),], aes(x=Date, y=lwp.m, col=factor(treatment), pch=water.status)) + geom_point() + geom_smooth(se=F, col=treatment)


ggplot(enchilada, aes(x=DOY, y=gs, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179)


ggplot(enchilada, aes(x=DOY, y=E, col=factor(treatment), pch=water.status)) + geom_point() + geom_vline(xintercept=179) + geom_vline(xintercept=200) + geom_smooth(aes(col=tag))


# visualizing individual WPs in different panels
# quartz(width=5, height=6) #lwp
# par(mfcol=c(5,2), oma=c(4,4,3,0.1),mar=c(1,1,0.1,0.1))
# for(i in unique(enchilada$tag)){
#   plot(lwp.m~Date, enchilada[which(enchilada$tag==i & enchilada$lwp.m<0),], col=treatment, ylim=c(-6,0), xlim=as.Date(c("2018-04-27","2018-07-30")))
# } 
# mtext(outer=TRUE, side = 3, "lwp")
# 
# 
# quartz(width=5, height=6) #aba

# par(mfcol=c(5,2), oma=c(4,4,3,0.1),mar=c(1,1,0.1,0.1))
# for(i in unique(enchilada$tag)){
#   plot(ABAFWngg~Date, enchilada[which(enchilada$tag==i),], col=treatment, xlim=as.Date(c("2018-04-27","2018-07-30")))
# }

# mtext(outer=TRUE, side = 3, "ABA")

quartz(width=5, height=6)
ggplot(enchilada, aes(x=Date, y=ABAFWngg,col=treatment)) + geom_point() + facet_wrap(facets = ~tag)

quartz(width=5, height=6)
ggplot(enchilada, aes(x=Date, y=lwp.m,col=treatment)) + geom_point() + facet_wrap(facets = ~tag)



### Plotting scatterplots of GS, ABA and lwp
quartz(width=5, height=6)
ggplot(gsaba, aes(x=log(ABAFWngg), y=log(gs), col=treatment)) + geom_point() #+ geom_smooth(se=F,method = "loess",span=1)

# same plot, minus pre-treatment
quartz(width=5, height=6)
ggplot(gsaba[which(gsaba$Date>"2018-06-28"),], aes(x=log(ABAFWngg), y=log(gs), col=treatment)) + geom_point() #+ geom_smooth(se=F,method = "loess",span=1)


# Gs vs LWP
quartz(width=5, height=6)
ggplot(enchilada[which(enchilada$ABAFWngg>0),], aes(x=lwp.m, y=log(gs), col=log(ABAFWngg) )) + geom_point() 

#quartz(width=5, height=6)
#ggplot(enchilada[which(enchilada$ABAFWngg>0),], aes(x=minlwp, y=log(gs), col=log(ABAFWngg) )) + geom_point()


quartz(width=5, height=6)
#plot(lwp.m~ABAFWngg, enchilada, col=treatment)
ggplot(enchilada[which(enchilada$ABAFWngg>0 & enchilada$Date>"2018-06-28" ),], aes(x=ABAFWngg, y=lwp.m, col=treatment, size=time.since.rewater+1 )) + geom_point()


















# ++++++++++++++++++++++++++++++++++++++++++
######## average to treatment, rob's code #########

quartz(width=5, height=6)
ench <- enchilada %>% group_by(Date, DOY, treatment) %>% summarise(E.m = mean(E, na.rm=T), gs.m=mean(gs,na.rm=T), lwp.m = mean(lwp.m, na.rm=T), ABA.m=mean(ABAFWngg, na.rm=T))


quartz(width=5, height=6)
par(mfrow=c(3,1), mar=c(0,5,0,1), oma=c(3,2,1,0), mgp=c(3,1,0))
inds <- unique(gasexchange$ID)
length(inds)
cols <- c("red","brown","darkgreen","violet","blue","cornflowerblue","darkorange","black","pink")
plot(Amax~Date, gasexchange, type="n", xaxt="n", xlab="",ylim=c(0,30))
k <- 0
for(i in unique(gasexchange$ID)){
  k <- k + 1
  lines(Amax~Date, gasexchange[which(gasexchange$ID==i),], col=cols[k])
}
plot(gs~Date, gasexchange, type="n", xaxt="n", xlab="",ylim=c(0,1000))
k <- 0
for(i in unique(gasexchange$ID)){
  k <- k + 1
  lines(gs~Date, gasexchange[which(gasexchange$ID==i),], col=cols[k])
}
plot(E~Date, gasexchange, type="n", xlab="", col = 1,ylim=c(0,10))
k <- 0
for(i in unique(gasexchange$ID)){
  k <- k + 1
  lines(E~Date, gasexchange[which(gasexchange$ID==i),], col=cols[k])
}
# group by treatment
quartz(width=5, height=6)
par(mfrow=c(3,1), mar=c(0,5,0,1), oma=c(4,1,1,0), mgp=c(3,1,0))
inds <- unique(gasexchange$ID)
treats <- c("control","mwd","swd")
control <- c(2,7,8)
mwd <- c(3,4,9)
swd <- c(5,1,6)
cols <- c("red","brown","darkgreen","violet","blue","cornflowerblue","darkorange","black","pink")
treatment <- 1:length(gasexchange$ID)
for (i in c(1:length(gasexchange$ID))){
  treatment[i] <- ifelse(gasexchange$ID[i] %in% inds[control],"control",
                         ifelse(gasexchange$ID[i] %in% inds[mwd], "mwd","swd"))
}
gasexchange <- cbind(gasexchange,treatment)
gasexchange <- as.data.frame(gasexchange)
gasexchange[,1] <- as.Date(gasexchange[,1])
gasexchange[,5] <- as.numeric(gasexchange[,5])
gasexchange[,6] <- as.numeric(gasexchange[,6])
gasexchange[,7] <- as.numeric(gasexchange[,7])
gasexchange[,8] <- as.numeric(gasexchange[,8])
gasexchange[,9] <- as.factor(gasexchange[,9])
ge.data <- matrix(NA,ncol=8,nrow=3*length(unique(gasexchange[,1])))
ge.data <- as.data.frame(ge.data)
ge.data[,1] <- as.character(ge.data[,1])
ge.data[,2] <- as.Date(ge.data[,2])
ge.data[,3] <- as.numeric(ge.data[,3])
ge.data[,4] <- as.numeric(ge.data[,4])
ge.data[,5] <- as.numeric(ge.data[,5])
ge.data[,6] <- as.numeric(ge.data[,6])
ge.data[,7] <- as.numeric(ge.data[,7])
ge.data[,8] <- as.numeric(ge.data[,8])
dates <- unique(gasexchange[,1])
j <- 0
ge.data
for (i in c(1:length(unique(gasexchange$Date)))){
  vals <- tapply(gasexchange$Amax[gasexchange[,1] == dates[i]],gasexchange$treatment[gasexchange[,1] == dates[i]],mean,na.rm=T)
  ge.data[c((j+i):c(j+i+2)),1]  <- names(vals)
  dat <- rep(dates[i],3)
  ge.data[c((j+i):c(j+i+2)),2] <- dat
  ge.data[c((j+i):c(j+i+2)),3] <- vals[1:3]
  vals <- tapply(gasexchange$Amax[gasexchange[,1] == dates[i]],gasexchange$treatment[gasexchange[,1] == dates[i]],se)
  ge.data[c((j+i):c(j+i+2)),4] <- vals[1:3]
  vals <- tapply(gasexchange$gs[gasexchange[,1] == dates[i]],gasexchange$treatment[gasexchange[,1] == dates[i]],mean,na.rm=T)
  ge.data[c((j+i):c(j+i+2)),5] <- vals[1:3]
  vals <- tapply(gasexchange$gs[gasexchange[,1] == dates[i]],gasexchange$treatment[gasexchange[,1] == dates[i]],se)
  ge.data[c((j+i):c(j+i+2)),6] <- vals[1:3]
  vals <- tapply(gasexchange$E[gasexchange[,1] == dates[i]],gasexchange$treatment[gasexchange[,1] == dates[i]],mean,na.rm=T)
  ge.data[c((j+i):c(j+i+2)),7] <- vals[1:3]
  vals <- tapply(gasexchange$E[gasexchange[,1] == dates[i]],gasexchange$treatment[gasexchange[,1] == dates[i]],se)
  ge.data[c((j+i):c(j+i+2)),8] <- vals[1:3]
  j <- j + 2
}
ge.data
colnames(ge.data) <- c("Treatment","Date","Amax","A.se","gs","gs.se","E","E.se")
names(ge.data)
par(mfrow=c(3,1), mar=c(0,5,0,0.1), oma=c(7,1,1,0), mgp=c(3,1,0))
cols <- c("deepskyblue4","darkgoldenrod3","red")
pchs <- c(16,15,1)
plot(Amax~Date,ge.data, type="n", xaxt="n", xlab="",ylim=c(0,20),col=cols,las=1,xlim=c(as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[36]),2],na.rm=T)),max(as.numeric(ge.data[,2]),na.rm=T)),ylab=expression(paste("Amax (", mu,"mol ",m^-2," ", s^-1,")")))
abline(v=dates[14] + 0.5,lty=2,col="black",lwd=2)
legend(dates[13] + 1,20,pch=pchs,col=cols,c(unique(ge.data$Treatment)))
k <- 0
for(i in unique(ge.data$Treatment)){
  k <- k + 1
  points(ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),3], col=cols[k],pch=pchs[k])
  arrows(ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),3]-ge.data[which(ge.data$Treatment==i),4],ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),3]+ge.data[which(ge.data$Treatment==i),4],code=3,length=0.01, col=cols[k],angle=90)
}

rw <- seq(from=as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[79]),2],na.rm=T))-1,to=max(as.numeric(ge.data[90,2]),na.rm=T), by=1)

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
rect(rw[1],-6,rw[4],30,col=addTrans("cornflowerblue",80),border = NA)
plot(gs~Date, ge.data, type="n", xaxt="n", xlab="",ylim=c(0,500),las=1,xlim=c(as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[36]),2],na.rm=T)),max(as.numeric(ge.data[,2]),na.rm=T)),ylab = expression(paste(g[s], " (mmol ",m^-2," ", s^-1,")")))
abline(v=dates[14] + 0.5,lty=2,col="black",lwd=2)
k <- 0
for(i in unique(ge.data$Treatment)){
  k <- k + 1
  points(ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),5], col=cols[k],pch=pchs[k])
  arrows(ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),5]-ge.data[which(ge.data$Treatment==i),6],ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),5]+ge.data[which(ge.data$Treatment==i),6],code=3,length=0.01, col=cols[k],angle=90)
}
rw <- seq(from=as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[79]),2],na.rm=T))-1,to=max(as.numeric(ge.data[90,2]),na.rm=T), by=1)

rect(rw[1],-6,rw[4],550,col=addTrans("cornflowerblue",80),border = NA)
plot(E~Date, ge.data, type="n", xlab="", col = 1,ylim=c(0,10),xlim=c(as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[36]),2],na.rm=T)),max(as.numeric(ge.data[,2]),na.rm=T)),xaxt="n",las=1,ylab=expression(paste("E (mol ",m^-2," ",s^-1,")")))
abline(v=dates[14] + 0.5,lty=2,col="black",lwd=2)
k <- 0
for(i in unique(ge.data$Treatment)){
  k <- k + 1
  points(ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),7], col=cols[k],pch=pchs[k])
  arrows(ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),7]-ge.data[which(ge.data$Treatment==i),8],ge.data[which(ge.data$Treatment==i),2],ge.data[which(ge.data$Treatment==i),7]+ge.data[which(ge.data$Treatment==i),8],code=3,length=0.01, col=cols[k],angle=90)
}
doydate <- read.csv("doydate.csv",header=T)
h <- strptime(as.character(doydate$date), "%d-%b-%y")
h <- format(h, format="%Y-%m-%d")
h <- as.character(h)
g <- as.character(ge.data[which(ge.data$Date > ge.data$Date[36]),2])
g <- unique(g)
lab.dates <- doydate[match(g,h),2]
vals.h.min <- match(g[1],h) 
these <- seq(from=as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[36]),2],na.rm=T)),to=max(as.numeric(ge.data[,2]),na.rm=T), by=1)
vals.h.max <- vals.h.min + length(these)
gnew <- doydate[vals.h.min:vals.h.max,2]

at.h <- seq(from=min(these),to=max(these)+2,by=3)
at.l <- seq(from=1,to=length(these),by=3)
those <- which(these %in% at.h)
axis(side=1,labels=gnew[at.l], at=these[at.l],tick=T,las=2)
mtext(side=1,line=6,"Date",at=0.5,outer=T)

rw <- seq(from=as.numeric(min(ge.data[which(ge.data$Date > ge.data$Date[79]),2],na.rm=T))-1,to=max(as.numeric(ge.data[90,2]),na.rm=T), by=1)
rect(rw[1],-6,rw[4],12,col=addTrans("cornflowerblue",80),border = NA)

