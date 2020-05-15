
#+++++++++++++++++++++++++++++++++++++++++++++++++++++
#          STEP 2: LOAD AND CLEAN ABBA DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++





# ABA data
# Buckeye 2018

#setwd("~/Documents/RESEARCH/California Buckeye")
# require(dplyr)

aba <- read.table("ABA.csv", header=T, sep=",") %>% select(1:3)
abadate.raw <- read.table("ABAdate.csv", header=T, sep=",")
head(aba)
head(abadate.raw)
abadate <- abadate.raw %>% group_by(Date) %>% summarise (Collection_num = unique(Collection_num))
abadate$Date <- as.Date(abadate$Date,"%F")
# Load functions
se <- function(x) sqrt(var(x, na.rm=T)/sum(!is.na(x)))
treatment <- rep(NA, length(aba[,1]))
aba <- cbind(aba, treatment)
control <- c(589, 498, 954,159)
mild <- c(592, 999, 600)
severe <- c(590, 588, 591)
aba$treatment[aba$tag %in% control] <- "control"
aba$treatment[aba$tag %in% mild] <- "mild"
aba$treatment[aba$tag %in% severe] <- "severe"
aba$date <- abadate$Date[match(aba$time.since.start, abadate$Collection_num)]


######### END: load ABA data #############################




########## ^^ RUN EVERYTHING ABOVE #########################
#____________________________________________________________












####### Rob's old Code ###############3
# commented out by LDLA on 04/06/2020
# 
# ## plotting time course of ABA for each individual #######
# cols <- c("purple","red","orange","saddlebrown","black","blue","green","cornflowerblue","pink","yellow")
# tag.i <- as.vector(unique(aba$tag))
# names(aba)
# 
#  tiff(filename = "FigX ABA.tif",
#       width = 8, height = 10, units = "cm", pointsize = 6,
#       compression = c("none"),
#       bg = "white", res = 1200, family = "", restoreConsole = TRUE,
#       type = c("windows", "cairo"),antialias = "cleartype")
# 
# par(mfrow=c(4,3),oma=c(4,4,0.1,0.1),mar=c(1,1,0.1,0.1))
# for (i in c(1:length(tag.i))){
#   plot(aba[aba$tag == tag.i[i],1]~aba[aba$tag == tag.i[i],3], type="b",col=cols[i],pch=16,xlim=c(1,22), ylim=c(0,360))
#   # rect(0,-50,7,400,col=cols[i])
#   text(1,340,LETTERS[i])
#   text(13,240,tag.i[i])
# }
# 
# mtext(side=2,outer=T, line=2,expression(paste("ABA level (g ",g^-1,")")))
# mtext(side=1,outer=T, line=2,"Day of sampling")
# # dev.off()
# 
# 
# 
# 
# 
# 
# 
# # treatment means
# aba.means <- tapply(aba$ABAFWngg,list(aba$treatment,aba$time),mean, na.rm=T)
# aba.means.se <- tapply(aba$ABAFWngg,list(aba$treatment,aba$time),se)
# time <- sort(unique(aba$time))
# cols <- c("purple","orange","grey")
# #tiff(filename = "FigX ABA treats.tif",
# #    width = 8, height = 10, units = "cm", pointsize = 6,
# #     compression = c("none"),
# #     bg = "white", res = 1200, family = "", restoreConsole = TRUE,
# #     type = c("windows", "cairo"),antialias = "cleartype")
# par(new=F)
# par(mfrow=c(1,1),oma=c(4,4,0.1,0.1),mar=c(1,1,0.1,0.1))
# for (i in c(1:3)){
#   plot(aba.means[i,]~time, type="b",col=cols[i],pch=16,xlim=c(1,22), ylim=c(0,250))
#   arrows(time, aba.means[i,]-aba.means.se[i,],time,aba.means[i,]+aba.means.se[i,], code=3, angle=90, length =0.005, col=cols[i])
#   par(new=T)
# }
# rect(0,-50,6,400,col=cols[i])
# text(1,240,LETTERS[1])
# legend(10,220,pch=16, c("Control", "Mild", "Severe"), col=cols)
# mtext(side=2,outer=T, line=2,expression(paste("ABA level (g ",g^-1,")")))
# mtext(side=1,outer=T, line=2,"Day of sampling")
# #dev.off()
# 

