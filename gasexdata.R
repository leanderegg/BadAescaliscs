

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



####### Combine ABA and LWP data w/ gasex #############3

## aba from ABA.R (Step 2)

# quick rename "999" and "954" to "10###"
aba$tag[which(aba$tag == 999)] <- 10999
aba$tag[which(aba$tag == 954)] <- 10954
# aba$tag <- as.numeric(aba$tag) # make numeric to match gasexchange, which as a numeric 'tag' column
gsaba <- left_join(gasexchange, aba, by=c("tag"="tag","Date"="date", "treatment"="treatment"))
gsaba$time.since.start <- NULL
head(gsaba)





## wps.clean from Aescal_Greenhous_WP_Select.R (Step 1)

# with Rob's old janky ass dataframe
# enchilada <- left_join(gsaba, lwppd.all %>% select(doy.yr, tag, lwppd_MPa), by = c("Date"="doy.yr","tag"="tag"))


# with Lee's new method
# load in most recent version of wps.clean so don't have to run WP_Select code
wps.clean <- read.csv("WP_data_processed20200327.csv")[,-1]
wps.clean$DOY.yr <- as.Date(wps.clean$DOY.yr) # make sure DOY.yr is in Date format


# merge gs & ABA data with wp data
enchilada <- left_join(gsaba, wps.clean %>% select(DOY.yr, tag, lwp.m), by = c("Date"="DOY.yr","tag"="tag"))

enchilada$treatment <- factor(enchilada$treatment)
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






# 
# 
# 
# ####### Adding some derived metrics ############
# 
# # set a minimum water potential for each individual
# ench <- enchilada %>% group_by(ID, tag, treatment) %>% summarise(minlwp = min(lwp.m, na.rm=T), maxABA = max(ABAFWngg, na.rm=T), mings=min(gs, na.rm=T))
# 
# ench$date.min <- as.Date(NA)
# for(i in unique(enchilada$tag)){
#   tmp <- enchilada[which(enchilada$tag==i),]
#   ench$date.min[which(ench$tag==i)] <- as.Date(tmp$Date[which(tmp$lwp.m == min(tmp$lwp.m, na.rm=T))])
# }
# 
# # add in to enchilada a column indicating min wp, and date of min so we can analyze recover
# enchilada$minlwp <- ench$minlwp[match(enchilada$tag, ench$tag)]
# enchilada$date.min <- ench$date.min[match(enchilada$tag, ench$tag)]


#_________________________________________________________
####### Plot timecourses of Gs, LWP, and ABA #########
#_________________________________________________________

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
ggplot(enchilada[which(enchilada$ABAFWngg>0 & enchilada$Date>"2018-06-28"),], aes(x=ABAFWngg, y=lwp.m, col=treatment )) + geom_point()






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

