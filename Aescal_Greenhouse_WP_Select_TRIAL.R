#+++++++++++++++++++++++++++++++++++++++++++++++++++++
#          STEP 1: LOAD AND CLEAN WP DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++









#rm(list=ls())
# set working directory
#setwd("~/Documents/RESEARCH/California Buckeye/")
require(tidyr)
require(dplyr)
require(stringr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(ggplot2)



se <- function(x,na.rm=T) sqrt(var(x))/sum(!is.na(x))


#____________________________________________________________________________________________________________
####### START: Loak raw psychrometer data, qc flag data, and calibrate to get raw wp data #######################

wpdata<-read.csv("GREENHOUSEWPselect.csv", header = FALSE)
head(wpdata)
iteration.number <- 1:length(wpdata[,1])
iteration.number[1] <- 1
for (i in c(2:length(iteration.number))){
  iteration.number[i] <- ifelse(wpdata[i,1] == 106,iteration.number[i-1]+1,iteration.number[i-1])
}
wpdata.n <- cbind(iteration.number,wpdata)
numberofiterations <- max(iteration.number)/2

complete.data <- matrix(NA, ncol=14, nrow=numberofiterations)
complete.data <- as.data.frame(complete.data)
complete.data[,1] <- as.numeric(complete.data[,1])
complete.data[,2] <- as.numeric(complete.data[,2])
complete.data[,3] <- as.numeric(complete.data[,3])
complete.data[,4] <- as.numeric(complete.data[,4])
complete.data[,5] <- as.numeric(complete.data[,5])
complete.data[,6] <- as.numeric(complete.data[,6])
complete.data[,7] <- as.numeric(complete.data[,7])
complete.data[,8] <- as.numeric(complete.data[,8])
complete.data[,9] <- as.numeric(complete.data[,9])


colnames(complete.data) <- c("Yr","DOY","Time","Temp","PSY1","PSY2","PSY3","PSY4","PSY5","PSY.one1","PSY.one2","PSY.one3","PSY.one4","PSY.one5")

#make an additional matrix for storing QC flags for each psy for each measurement
qc1 <- data.frame(matrix(NA, ncol=5, nrow=numberofiterations))
names(qc1) <- c("PSY1.qc","PSY2.qc","PSY3.qc", "PSY4.qc", "PSY5.qc")


complete.data.2 <- matrix(NA, ncol=12, nrow=numberofiterations)
complete.data.2 <- as.data.frame(complete.data.2)
complete.data.2[,1] <- as.numeric(complete.data.2[,1])
complete.data.2[,2] <- as.numeric(complete.data.2[,2])
complete.data.2[,3] <- as.numeric(complete.data.2[,3])
complete.data.2[,4] <- as.numeric(complete.data.2[,4])
complete.data.2[,5] <- as.numeric(complete.data.2[,5])
complete.data.2[,6] <- as.numeric(complete.data.2[,6])
complete.data.2[,7] <- as.numeric(complete.data.2[,7])
complete.data.2[,8] <- as.numeric(complete.data.2[,8])
colnames(complete.data.2) <- c("Yr","DOY","Time","Temp","PSY6","PSY7","PSY8","PSY9","PSY.one6","PSY.one7","PSY.one8","PSY.one9")

#make an additional matrix for storing QC flags for each psy for each measurement
qc2 <- data.frame(matrix(NA, ncol=4, nrow=numberofiterations))
names(qc2) <- c("PSY6.qc","PSY7.qc","PSY8.qc", "PSY9.qc")


ninitial <- seq(from=1,to=max(iteration.number),by=1)
nadd <- seq(from=0,to=max(iteration.number)-1,by=1)
ncombine <- matrix(NA,ncol=3,nrow=max(iteration.number))
ncombine <- as.data.frame(ncombine)
ncombine[,1] <- ninitial
ncombine[,2] <- nadd
ncombine[,3] <- ncombine[,1] + ncombine[,2]



# run through, pull out each psy measurement, and average values from the peaks
QC <- F # T/F -run quality control script?
qc.ver <- "20180310" #quality control version

these <- which(ncombine[,3] < max(iteration.number))
those <- seq(from=2,to=max(iteration.number),by=2)
psy <- 1:5
psy.one <- 1:5
qc.flag <- rep("fail", times=5)
l <- 0
for (i in c(ncombine[these,3])){
  l <- l + 1
  yr <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),3] #pulling out yr from header row indicated by 106
  doy <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),4] # pulling out doy from header rows
  time <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),5] # pulling out time stamp from header rows (in hhmm format)
  temp <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),6] # pulling out ambient Temp from header row
  j <- 0
  for (k in c(109,111,113,115,117)){
    j <- j + 1
    all <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(6:8)]
    all <- as.numeric(all)
    psy[j] <- mean(all,na.rm=T)
    all <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(6)]
    all <- as.numeric(all)
    psy.one[j] <- mean(all,na.rm=T)
    if(QC==T){
      if(is.na(wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(5)])){
        qc.flag[j] <- NA 
      }
      else{
        plot(c(unlist(wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(5:14)])), typ="b", lwd=3)
        points(c(unlist(wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(6:8)]))~c(2,3,4),pch=16, col="red")
        print(paste("i=",i,"k=",k, "qc.row=",l,"qc.col=",j))
        qc.flag[j] <- askYesNo("Quality of trace?(1=good,2=kinda bad, 3=DELETE)", default=F, prompts=getOption("askYesNo", gettext(c("1","2","3"))))
      }
      
    }
  }
  complete.data[l,1] <- yr
  complete.data[l,2] <- doy
  complete.data[l,3] <- time
  complete.data[l,4] <- temp
  complete.data[l,c(5:9)] <- psy #store average of 3 value plateau
  complete.data[l,c(10:14)] <- psy.one # store value of first point of plateau (third measurement)
  if(QC==T){
    qc1[l, c(1:5)] <- qc.flag 
  }
}

# write qc output so you don't have to rerun:
# write.csv(qc1, paste0("Psychrometer_QC1_", qc.ver,".csv"))
# only run if QC ==T and you made a new qc file

qc1 <- read.csv("Psychrometer_QC1_20180309.csv")


### repeat for psychrometers 5-9
l <- 0
psy <- 1:4
psy.one <- 1:4

# restarting if interrupted in middle
# start right before you failed:
#l <- 
for (i in c(those)){
  l <- l + 1
  yr <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),3]
  doy <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),4]
  time <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),5]
  temp <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == 106),6]
  j <- 0
  qc.flag <- rep("fail", times=4)
  for (k in c(109,111,113,115)){
    j <- j + 1
    all <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(6:8)]
    all <- as.numeric(all)
    psy[j] <- mean(all,na.rm=T)
    all <- wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(6)]
    all <- as.numeric(all)
    psy.one[j] <- mean(all,na.rm=T)
    if(QC==T){
      if(is.na(wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(5)])){
        qc.flag[j] <- NA 
      }
      else{
        plot(c(unlist(wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(5:14)])), typ="b", lwd=3)
        points(c(unlist(wpdata.n[which(wpdata.n[,1] == i & wpdata.n[,2] == k),c(6:8)]))~c(2,3,4),pch=16, col="red")
        print(paste("i=",i,"k=",k, "qc.row=",l,"qc.col=",j))
        qc.flag[j] <- askYesNo("Quality of trace?(1=good,2=kinda bad, 3=DELETE)", default=F, prompts=getOption("askYesNo", gettext(c("1","2","3"))))
      }
      
    }
    
  }
  complete.data.2[l,1] <- yr
  complete.data.2[l,2] <- doy
  complete.data.2[l,3] <- time
  complete.data.2[l,4] <- temp
  complete.data.2[l,c(5:8)] <- psy
  complete.data.2[l,c(9:12)] <- psy.one
  if(QC==T){
    qc2[l, c(1:4)] <- qc.flag 
  }
}

# write.csv(qc2, paste0("Psychrometer_QC2_", qc.ver,".csv"))
qc2 <- read.csv("Psychrometer_QC2_20180309.csv")


complete.data[,c(5:9)] <- complete.data[,c(5:9)]/(0.325+0.027*complete.data[,4]) 
complete.data[,c(10:14)] <- complete.data[,c(10:14)]/(0.325+0.027*complete.data[,4]) 

complete.data.2[,c(5:8)] <- complete.data.2[,c(5:8)]/(0.325+0.027*complete.data.2[,4]) 
complete.data.2[,c(9:12)] <- complete.data.2[,c(9:12)]/(0.325+0.027*complete.data.2[,4]) 

calib <- read.table("PsyCalibrationSLOPEINT.txt",header=T)
# calibrate mean values for Psychrometers 1:9
complete.data[,5] <- (complete.data[,5]-calib[1,3])/calib[1,2]
complete.data[,6] <- (complete.data[,6]-calib[2,3])/calib[2,2]
complete.data[,7] <- (complete.data[,7]-calib[3,3])/calib[3,2]
complete.data[,8] <- (complete.data[,8]-calib[4,3])/calib[4,2]
complete.data[,9] <- (complete.data[,9]-calib[5,3])/calib[5,2]
complete.data.2[,5] <- (complete.data.2[,5]-calib[6,3])/calib[6,2]
complete.data.2[,6] <- (complete.data.2[,6]-calib[7,3])/calib[7,2]
complete.data.2[,7] <- (complete.data.2[,7]-calib[8,3])/calib[8,2]
complete.data.2[,8] <- (complete.data.2[,8]-calib[9,3])/calib[9,2]
# Values based on a single value
complete.data[,10] <- (complete.data[,10]-calib[1,3])/calib[1,2]
complete.data[,11] <- (complete.data[,11]-calib[2,3])/calib[2,2]
complete.data[,12] <- (complete.data[,12]-calib[3,3])/calib[3,2]
complete.data[,13] <- (complete.data[,13]-calib[4,3])/calib[4,2]
complete.data[,14] <- (complete.data[,14]-calib[5,3])/calib[5,2]
complete.data.2[,9] <- (complete.data.2[,9]-calib[6,3])/calib[6,2]
complete.data.2[,10] <- (complete.data.2[,10]-calib[7,3])/calib[7,2]
complete.data.2[,11] <- (complete.data.2[,11]-calib[8,3])/calib[8,2]
complete.data.2[,12] <- (complete.data.2[,12]-calib[9,3])/calib[9,2]




####### Quick visualizing raw output ###########################3
#raw predawn traces (for all psychrometers)
par(mfcol=c(5,2), oma=c(4,4,0.1,0.1),mar=c(1,1,0.1,0.1))
cols <- c("black","red","blue","darkgreen","brown","turquoise","violet","darkorange","cornflowerblue")
for (i in c(5:9)){
  plot(complete.data[complete.data[,3] < 500 & complete.data[,3] > 400,i]~complete.data[complete.data[,3] < 500 & complete.data[,3] > 400,2],ylim=c(-6,0),las=1, col=cols[i-4])
}
for (i in c(5:8)){
  plot(complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,i]~complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,2],ylim=c(-6,0),las=1, col=cols[i+1])
}

# #raw predawn traces (for all psychrometers) based on a single value (max of 3 points)
# par(mfcol=c(5,2), oma=c(4,4,0.1,0.1),mar=c(1,1,0.1,0.1))
# cols <- c("black","red","blue","darkgreen","brown","turquoise","violet","darkorange","cornflowerblue")
# for (i in c(10:14)){
#   plot(complete.data[complete.data[,3] < 500 & complete.data[,3] > 400,i]~complete.data[complete.data[,3] < 500 & complete.data[,3] > 400,2],ylim=c(-6,0),las=1, col=cols[i-9])
# }
# for (i in c(9:12)){
#   plot(complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,i]~complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,2],ylim=c(-6,0),las=1, col=cols[i-3])
# }
# head(complete.data)



####### END: load raw psychrometer data and quality flags #############



#________________________________________________________________________
#### START: Lee's new way of processing wp data to long form ########



# step 1: take wide data and make it long
# step 2: merge in quality control flags
# step 3: aggregate things by date, filtering for bad values


# take complete.data and make it long form (take 5 psy columns and make them 1 column)
wpslong1 <- gather(select(.data = complete.data, Yr,DOY,Time,Temp,PSY1:PSY5) , key = "Psy_num", value="lwppd", PSY1:PSY5)

# take qc1 and also make it long
names(qc1)<- str_replace(names(qc1), pattern = ".qc","") # pull out the ".qc" part of the quality control data flag column names so they will match up with wpslong psychrometer names
qualcont1 <- data.frame(complete.data[,1:4], qc1[,-1]) # tacking on date, time, temp info from complete.data and killig extra column named "X" that came in when we read.csv (was rownmanes)
qclong1 <- gather(qualcont1,key = "Psy_num", value="QC_flag", PSY1:PSY5 ) # take wide form and make long

#merge wps and qc flags for first 5 psychrometers
wps1 <- full_join(wpslong1, qclong1)



# repeat the above steps for psychrometers 6-9
# take complete.data and make it long form (take 5 psy columns and make them 1 column)
wpslong2 <- gather(select(.data = complete.data.2, Yr,DOY,Time,Temp,PSY6:PSY9) , key = "Psy_num", value="lwppd", PSY6:PSY9)

# take qc2 and also make it long
names(qc2)<- str_replace(names(qc2), pattern = ".qc","")
qualcont2 <- data.frame(complete.data.2[,1:4], qc2[,-1])
qclong2 <- gather(qualcont2,key = "Psy_num", value="QC_flag", PSY6:PSY9 )

#merge wps and qc flags
wps2 <- full_join(wpslong2, qclong2)


# row bind psychrometers 1-5 and 6-9 together into one dataframe
wps.all <- rbind(wps1, wps2)



# now summarize all repeat measurments to get a mean per day. also filter out bad psychrmenter readings
wps.clean <- wps.all %>% filter(QC_flag==TRUE | QC_flag==FALSE) %>% group_by(Yr, DOY, Psy_num) %>% summarise(lwp.m = mean(lwppd, na.rm=T), lwp.sd = sd(lwppd, na.rm=T), lwp.n = n(), Temp.m = mean(Temp), Temp.sd = sd(Temp))



###### . Adding tag and treatment info ############
# extract individuals for matching up with tags
wps.clean$ind <- as.numeric(str_replace(wps.clean$Psy_num,"PSY", ""))



### add in a 'tag' column to match rest of dataset:
tag <- c(590,592,589,498,10999,588,591,10954,600)
ind <- c(1:9)
treatment <- c("swd", "mwd", "control", "control","mwd","swd","swd","control","mwd")
tag.mapper <- data.frame(tag, ind, treatment)

# add final column to wp data with individuals tag name and treatment
wps.clean$tag <- tag.mapper$tag[match(wps.clean$ind, tag.mapper$ind)]
wps.clean$treatment <- tag.mapper$treatment[match(wps.clean$ind, tag.mapper$ind)]


# convert DOY to dates so we can match
wps.clean$doy.yr <- paste(wps.clean$DOY, "18", sep="-")
  # make a new column that R recognizes as a date
wps.clean$DOY.yr <- as.Date(wps.clean$doy.yr, "%j-%y")



###### END: average wps and turn long form ###################
#____________________________________________________________________________




####### ^^ RUN EVERYTHING ABOVE ^^ ######################









#______________________________________________________________________________________

################## Rob's way of cleaning and processing data ##########################

# DONT RUN !!!!


######## taking wide form data and making it long form (and combining psys 1-5 &6-9) ##########

# make empty dataframe to put averaged values into
lwppd.all <- matrix(NA,ncol=4,nrow=9*length(unique(complete.data[complete.data[,3] < 500 & complete.data[,3] > 400,2])))
lwppd.all <- as.data.frame(lwppd.all)
colnames(lwppd.all) <- c("doy","ind","lwppd_MPa","treatment")



#  # of days with predawns
ks <- length(unique(complete.data[complete.data[,3] < 500 & complete.data[,3] > 400,2]))

k <- 1
clean <- c(1,3,5,7,9)
check <- c(6)
alsocheck <- c(8)
# loop through psychrometer columns, and pull out the max wp
for (i in c(5:9)){
  # take the max lwp from each day of measurments for a single psychrometer
  max.wps <- tapply(complete.data[complete.data$Time < 500 & complete.data$Time > 400,i],complete.data[complete.data$Time < 500 & complete.data$Time > 400,2],max,na.rm=T) # any time we have no good readings for a day, returns -Inf
  max.wps <- as.vector(max.wps)
  # pull out day of year for each measurement
  sampling.doys <- tapply(complete.data[complete.data$Time < 500 & complete.data$Time > 400,2],complete.data[complete.data$Time < 500 & complete.data$Time > 400,2],max,na.rm=T)
  
  # excluding wps more negative than -1.75 in pre-treatment period
  include <- ifelse((i-4) %in% clean & sampling.doys < 188,which(max.wps > -1.75),c(1:length(max.wps)))
    # if psychrometer is 1,
  max.wps <- max.wps[include]
  lwppd.all[c(k:c(k+(ks-1))),3] <- max.wps
  lwppd.all[c(k:(k+ (ks-1))),1] <- sampling.doys
  lwppd.all[c(k:(k + (ks-1))),2] <- rep(i-4,length(max.wps))
  k <- k + ks
}
lwppd.all[which(lwppd.all[,3] > 0),3] <- NA # get rid of bad wp that are positive
lwppd.all$lwppd_MPa[which(lwppd.all$lwppd_MPa< -10)] <- NA # put NAs in places that returned -Infs

for (i in c(5:8)){
  max.wps <- tapply(complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,i],complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,2],mean,na.rm=T)
  max.wps <- as.vector(max.wps)
  sampling.doys <- tapply(complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,2],complete.data.2[complete.data.2[,3] < 500 & complete.data.2[,3] > 400,2],mean,na.rm=T)
  include <- ifelse((i+1) %in% clean & sampling.doys < 188,which(max.wps > -1.75),c(1:length(max.wps)))
  include <- ifelse((i+1) %in% check & sampling.doys < 188,which(max.wps > -1.75),include)
  include <- ifelse((i+1) %in% alsocheck & sampling.doys < 189,which(max.wps > -1.75),include) 
  max.wps <- max.wps[include]
  lwppd.all[c(k:(k+(ks-1))),3] <- max.wps
  lwppd.all[c(k:(k+(ks-1))),1] <- sampling.doys
  lwppd.all[c(k:(k+(ks-1))),2] <- rep(i+1,length(max.wps))
  k <- k + ks
}
lwppd.all[which(lwppd.all[,3] > 0),3] <- NA


# group by treatment
treats <- c("control","mwd","swd")
control <- c(3,4,8)
mwd <- c(2,5,9)
swd <- c(1,6,7)
for (i in c(1:length(lwppd.all$ind))){
  lwppd.all[i,4] <- ifelse(lwppd.all[i,2] %in% control,"control",
                           ifelse(lwppd.all[i,2] %in% mwd, "mwd","swd"))
}
lwppd.all[which(lwppd.all[,3] == -Inf),3] <- NA




### add in a 'tag' column to match rest of dataset:
tag <- c(590,592,589,498,10999,588,591,10954,600)
ind <- c(1:9)
treatment <- c("swd", "mwd", "control", "control","mwd","swd","swd","control","mwd")
tag.mapper <- data.frame(tag, ind, treatment)

### turning DOY into a date that plays well with other datasets

lwppd.all$tag <- tag.mapper$tag[match(lwppd.all$ind, tag.mapper$ind)]
lwppd.all$doy.yr <- paste(lwppd.all$doy, "18", sep="-")
lwppd.all$doy.yr <- as.Date(lwppd.all$doy.yr, "%j-%y")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## Avveraging to treatment #############
mean.data <- tapply(lwppd.all[,3],list(lwppd.all[,1],lwppd.all[,4]),mean,na.rm=T)
mean.data.se <- tapply(lwppd.all[,3],list(lwppd.all[,1],lwppd.all[,4]),se)


# plotting...some shit?

par(mfrow=c(3,3))
for (i in c(1:9)){
  plot(lwppd.all[lwppd.all$ind == i,3]~lwppd.all[lwppd.all$ind == i,1], ylim=c(-5,0),yaxt="n",xaxt="n",type="l")
  axis(side=1,tick=T,labels=F)
  axis(side=2,tick=T,labels=F)
}

# plotting different shit
pchs <- c(16,15,1)
cols <- c("deepskyblue4","darkgoldenrod3","red")
par(new=F,mfrow=c(1,1),oma=c(4,3,0.1,0.1),mar=c(3,1,0.1,0.1))
doys <- as.numeric(rownames(mean.data))
for (i in c(1:3)){
  plot(mean.data[,i]~doys,ylim=c(-5,0),las=1, pch=pchs[i],col=cols[i],xlab="",ylab="",xlim=c(167,208),xaxt="n")
  lines(mean.data[,i]~doys,ylim=c(-5,0),las=1,col=cols[i],xlab="",ylab="",xaxt="n")
  arrows(doys,mean.data[,i]+mean.data.se[,i],doys,mean.data[,i]-mean.data.se[,i],code=3,length=0.01, col=cols[i],angle=90)
  par(new=T)
}
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
rect(198,-6,202,1,col=addTrans("cornflowerblue",80),border = NA)
doydate <- read.csv("doydate.csv",header=T)
lab.dates <- seq(from=168,to=208,by=3)
these <- match(lab.dates,doydate[,1])
axis(side=1,labels=doydate[these,2], at=lab.dates,tick=T,las=2)
abline(v=178.5,lty=2,lwd=3,col="black")
legend(169,-3,pch=pchs,col=cols,c("Control","MWD","SWD"))
mtext(side=2,outer=T,"Predawn LWP (MPa)",line=1.5,at=0.55)
mtext(side=1,outer=T,"Date",line=2.5,at=0.5)

