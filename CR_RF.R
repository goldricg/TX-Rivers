#Code to pull sub watersheds from ArcGIS and summarize information
#before running this code make sure to run 4 python codes in order written below:
#WTS_Delineation
#WTS_areacalc
#AttSum

#Written by Grace Goldrich Middaugh
#Written 11/7/2021
#Last Updated by Grace Goldrich-Middaugh
#Last Updated 06/03/2022

#Set working directories
setwd("C:/Users/grace/Box Sync/Oregonstate/TX/ColoradoSites") #windows
wname <- getwd() # set back to working directory
dname <- paste(wname,"Data",sep="/") # will open data file in working directory
Oname <- paste(wname,"Outputs",sep="/") # will open outputs file in working directory
fname<-paste(wname,"Plots/Chem",sep="/") # will open plots file in working directory
xname<-paste(wname,"Functions",sep="/") # will open functions file in working directory

#load packages
library(tidyverse)
library(ggplot2)
library(reshape2)
library(readxl)
library(GGally)
library(randomForest)
library(caTools)
library(arcgisbinding)
arc.check_product() #set up connection

#Load data
setwd(dname)
psum <- readRDS("CS_SimpleSum.RDS")
CS <- readRDS("CS.data.long.RDS")
wqpc <- readRDS("CSwqpchargebalance.RDS")
wqpc$Year <- format(as.Date(wqpc$date, format="%Y-%m-%d"), "%Y")
data <- readRDS("Allwts_chemdata041922.RDS")

#add NEAR_FID to wqp data
CS <- readRDS("CS.data.long.RDS")
CSunique <- CS %>% 
  group_by(siteno, NEAR_FID) %>% 
  filter(row_number() == 1)
CSunique <- na.omit(as.data.frame(CSunique[,c(2,9)]))
wqpc <- wqpc[!is.na(wqpc$siteno),]
wqpc$NEAR <- NA
for(i in 1:nrow(wqpc)){
  for(j in 1:nrow(CSunique)){
    if(wqpc$siteno[i] == CSunique$siteno[j]){
      wqpc$NEAR[i] <- CSunique$NEAR_FID[j]
    }
  }
}

#process and convert to compositional
wqp <- na.omit(reshape2::melt(wqpc, id.vars=c("date", "siteno", "lat", "long", "CBerror", "NEAR")))
wqp <- split(wqp, wqp$variable) #remove Phosphorus, Nitrate, Alkalinity, and Flow
wqp <- do.call(rbind,wqp)
wqp <- wqp[!(wqp$variable == "Phosphorus"|wqp$variable == "Alkalinity"|wqp$variable == "Flow"|
               wqp$variable == "Aluminum"|wqp$variable == "Bromide"|wqp$variable == "Iron"|wqp$variable == "Nitrate"),]
wqp$value <- as.numeric(wqp$value)
wqp <- reshape2::dcast(wqp, date+siteno+lat+long+NEAR~variable, value.var="value", mean)

param <- c("Bicarbonate","Calcium", "Chloride", "Potassium", "Magnesium", "Sodium", "Sulfate", "Silica")

wqp <- wqp[complete.cases(wqp[,param]),] #4863 complete records
length(unique(wqp$siteno)) #actually only 78 sites with complete records
wqcomp <- wqp[,1:5]

wqcomp[,param] <- wqp[,param]/rowSums(wqp[,param])
rowSums(wqcomp[,param])

#summary dataframe
site_name <- c("Site1", "Site2", "Site3", "Site4", "Site5", "Site6", "Site7", "Site8", "Site9",
               "Site10", "Site11", "Site12", "Site13", "Site14", "Site15")
numbers <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
site_no <- c("08117995", "08119500", "08121000", "08123850", "08126380", "08136700", "08147000", "08158000",
             "08159200", "08159500", "08160400", "08161000", "08162000", "08162500", "08162501")
segment <- c("1413", "1412", "1412", "1412", "1426", "1410", "1409",
             "1428", "1434", "1434", "1434", "1402", "1402", "1402", "1401")
lat <- c(32.62861, 32.53915, 32.39250, 32.05361, 31.71528, 31.49361, 31.21778, 30.24614,
         30.10444, 30.01250, 29.91222, 29.70611, 29.30889, 28.97389, 28.77417)
long <- c(-101.285000, -101.055222, -100.878333, -100.761667, -100.026111, -99.573611, -98.564167,
          -97.680056, -97.319167, -97.161667, -96.903611, -96.536667, -96.103611, -96.012222, -95.997500)
colodat <- data.frame(site_name, numbers,site_no, segment, lat, long)
arc.write(file.path("C:/Users/grace/Box Sync/Oregonstate/TX/Spatial/All_WTS\\All_WTS.gdb/CRsampsites"), data=colodat, coords=c('long', 'lat'), shape_info=list(type='Point', WKID=4269), overwrite=TRUE)

sites <- arc.open(file.path("C:/Users/grace/Box Sync/Oregonstate/TX/Spatial/All_WTS\\All_WTS.gdb", "old_CRsites"))
sites <- arc.select(sites, fields=c("siteno", "lat", "long", "trib", "id"))
psum <- inner_join(sites, psum, by=c("id" = "ID"))

wdata <- wqcomp[wqcomp$siteno %in% psum$siteno,] #subset chem data to match 78 sites
wdata <- left_join(wdata,psum, by=c("siteno" = "siteno"))

wdata$yday <- lubridate::yday(format(as.Date(wdata$date, format="%Y-%m-%d")))
wdata$syday <- sin(wdata$yday)
data <- data[which(data$River.x.x == "Colorado"),]
data <- distinct(data, siteno,.keep_all = TRUE)
wdata <- left_join(wdata, data, by=c("siteno" = "siteno"))

##### Random Forest prediction ####

#Split into test and train datasets
#center and standardize
param <- c("Bicarbonate.x", "Calcium.x",  "Chloride.x",  "Potassium.x", "Magnesium.x", "Sodium.x",     
            "Sulfate.x",  "Silica.x")
lithnames <- c("Sdep.x", "Conglomerate.x", "Carbonates.x", "Evaporites.x",
                         "Mudstone.x", "Igmet.x", "Sandstone.x", "MAP")
data <- as.data.frame(scale(wdata[,c(param, lithnames)], center=TRUE, scale = TRUE))
#colnames(data) <- c("Bicarbonate", "Calcium",  "Chloride",  "Potassium", "Magnesium", "Sodium",     
#                    "Sulfate",  "Silica", "Sdep", "Conglomerate", "Carbonates", "Evaporites",
#                    "Mudstone", "Igmet", "Sandstone", "MAP")
sample = sample.split(data$Bicarbonate, SplitRatio = .80)
train = subset(data, sample == TRUE)
test  = subset(data, sample == FALSE)
dim(train)
dim(test)

#run regression and save importance, predict for test dataset
i <- 1
fit <- list()
imp <- list()

for(i in 1:length(param)){
  train1 <- cbind(train[,param[i]], train[,lithnames]) #make dataset with all chem var and 1 lithology class
  colnames(train1) <- c(param[i], lithnames)
  
  forest_input<-formula(paste(param[i], "~."))
  fit[[i]] <- randomForest(forest_input,data=train1, importance=TRUE) #run regression rf
  imp[[i]] <- as.data.frame(importance(fit[[i]])) #save importance
  imp[[i]]$preds <- rownames(imp[[i]])
  imp[[i]]$param <- param[i]
  
  #use model to predict test dataset
  test1 <- cbind(test[,param[i]], test[,lithnames]) #make dataset with all chem var and 1 lithology class
  colnames(test1) <- c(param[i], lithnames)
  
  #test1[,lithnames[i]] <- factor(test1[,lithnames[i]]) #for classification
  pred = predict(fit[[i]], newdata=test1)
  
  #compare predicted to observed
  test1$pred <- pred
  fit[[i]]$cor <- cor(test1[,param[i]], test1$pred)
}

names(fit) <- param
names(imp) <- param


#extract info to dataframe
sum <- do.call(rbind, imp)
sum$spcor <- NA
i <- 1
j <- 1
for(i in 1:length(lithnames)){
  for(j in 1:length(param)){
    sum[which(sum$preds == lithnames[i] & sum$param == param[j]),]$spcor <- 
      cor(data[,param[j]], data[,lithnames[i]])
  }
}

sum <- sum[!sum$param %in% sum$preds,]

#extract performance on test data
cors <- c()
sols <- c()
x <- 1
for(i in 1:length(fit)){
  cors[x] <- fit[[i]]$cor
  sols[x] <- names(fit[i])
  x <- x+1
}
perf <- as.data.frame(cbind(round(cors,2), sols))
plot <- ggplot(data=perf, aes(x=factor(sols), y=round(cors,2)))+
  geom_point()+
  labs(x="Solute", y="Model performance")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)

#make lithology plot
plot <- ggplot(data=sum)+
  geom_point(aes(x=preds, y=param, size=`%IncMSE`, color=spcor))+
  labs(x="Lithology", y="Solutes", size="IncMSE", color="Corr")+
  scale_size(range=c(3,20))+
  scale_color_distiller(type="div", palette= "RdYlBu")+
  #scale_color_gradient2(low="red", mid="white", high="blue")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)

#plot performance of rf model for each 
####group land use classes ####
drops <- c("OpenWater.x", "Water.x")
wdata <- wdata[ , !(names(wdata) %in% drops)]
lcg <- wdata[,c("date.x","siteno","id","trib","yday", "syday","NEAR",param, lithnames)]
lulc <- c("Barren", "Open", "Forest", "CultivatedCrops", "PastureHay", "Developed")
lcg$Barren <- wdata$Barren.x
lcg$Open <- wdata$GrasslandHerbaceous+wdata$Shrubs+wdata$EmergentHerbaceousWetlands+wdata$WoodyWetlands
lcg$Forest <- wdata$DeciduousForest+wdata$EvergreenForest+wdata$MixedForest
lcg$CultivatedCrops <- wdata$CultivatedCrops
lcg$PastureHay <- wdata$PastureHay
lcg$Developed <- wdata$DevelopedHighIntenisy+wdata$DevelopedLowIntensity+wdata$DevelopedMediumIntensity+wdata$DevelopedOpenSpace

#split into test and train
lcg1 <- as.data.frame(scale(lcg[,c(param, lithnames, lulc)], center=TRUE, scale = TRUE))
sample = sample.split(lcg1$Bicarbonate, SplitRatio = .80)
train = subset(lcg1, sample == TRUE)
test  = subset(lcg1, sample == FALSE)

#do rf for grouped land use
preds <- c(lithnames, lulc)
i <- 1
fit1 <- list()
imp1 <- list()
for(i in 1:length(param)){
  train1 <- cbind(train[,preds], train[,param[i]]) #make dataset with all wts var and 1 chem
  colnames(train1) <- c(colnames(train[,preds]), param[i])
  
  forest_input <- formula(paste(param[i], "~."))
  fit1[[i]] <- randomForest(forest_input,data=train1, importance=TRUE) #run regression rf
  imp1[[i]] <- as.data.frame(importance(fit1[[i]])) #save importance
  imp1[[i]]$preds <- rownames(imp1[[i]])
  imp1[[i]]$param <- param[i]
  
  #use model to predict test dataset
  test1 <- cbind(test[,param[i]], test[,preds]) #make dataset with all wts var and 1 chem
  colnames(test1) <- c(param[i], preds)
  pred = predict(fit1[[i]], newdata=test1)
  
  #compare predicted to observed
  test1$pred <- pred
  fit1[[i]]$cor <- cor(test1[,param[i]], test1$pred)
}

names(fit1) <- param
names(imp1) <- param
z <- print(fit1[[1]])
#extract info to dataframe
sum <- do.call(rbind, imp1)
sum$spcor <- NA
i <- 1
j <- 1
for(i in 1:length(preds)){
  for(j in 1:length(param)){
    sum[which(sum$preds == preds[i] & sum$param == param[j]),]$spcor <- 
      cor(lcg1[,param[j]], lcg1[,preds[i]])
  }
}

sum <- sum[!sum$param %in% sum$preds,]

#extract performance on test data
cors <- c()
sols <- c()
mse <- c()
x <- 1
for(i in 1:length(fit1)){
  cors[x] <- fit1[[i]]$cor
  mse[x] <- mean(fit1[[i]]$mse)
  sols[x] <- names(fit1[i])
  x <- x+1
}

perf <- as.data.frame(cbind(round(cors,2),sols))
perf$mse <- mse
names(perf) <- c("Rsq", "Solute", "MSE")
str(perf)
perf$Rsq <- as.numeric(perf$Rsq)^2 #convert from correlation to R^2
plot <- ggplot(data=perf)+
  geom_point(aes(x=factor(Solute),y=Rsq), color="black")+
  geom_point(aes(x=factor(Solute),y=MSE), color="red")+
  labs(x="Solute", y="Model performance")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)

setwd(dname)
write.csv(perf, "RFperformance.csv")

#make summary plot
#identify values in top 85th percentile to label on plot
x <- quantile(sum$`%IncMSE`, probs=c(0.8))
sum$num <- ifelse(sum$`%IncMSE` > x, sum$`%IncMSE`, NA)
sum$preds <- gsub(".x", "", sum$preds)
sum$param <- gsub(".x", "", sum$param)
lithnames <- gsub(".x", "", lithnames)
sum$preds <- factor(sum$preds, levels=c(lithnames, lulc))

library(viridisLite)
setwd(fname)
jpeg("RF_outputsum.jpg", width=1000, height=800, quality=100)
plot <- ggplot(data=sum)+
  geom_point(aes(x=preds, y=param, size=`%IncMSE`, fill=spcor), pch=21, color="black")+
  labs(x=" ", y="Solute", size="% IncMSE", fill="Correlation")+
  geom_text(data=sum, aes(x=preds, y=param, label=round(num, 0)), 
            color="black", size=5)+
  scale_size_continuous(range=c(2,25), limits=range(sum$`%IncMSE`), breaks=c(5,10,20))+
  #scale_fill_viridis_c(option="magma", begin=0.22)+
  scale_color_distiller(type="div", palette= "RdYlBu", aesthetics=c("colour", "fill"))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 20))
print(plot)
dev.off()
##### Partial dependence plots #####
data(iris)
set.seed(543)
iris.rf <- randomForest(Species~., iris)
partialPlot(iris.rf, iris, Petal.Width, "versicolor")

## Looping over variables ranked by importance:
data(airquality)
airquality <- na.omit(airquality)
set.seed(131)
library(randomForest)
ozone.rf <- randomForest(Ozone ~ ., airquality, importance=TRUE)
imp <- importance(ozone.rf)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op <- par(mfrow=c(2, 3))
for (i in seq_along(impvar)) {
  partialPlot(ozone.rf, airquality, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]),
              ylim=c(30, 70))
}
par(op)


##### Analyses for supplemental info ####
#do RF for zones 2-5, 2-8, and 9-14 separately
#split into test and train
#subset to zones 2-5
lcgup <- lcg[which(lcg$trib == 1),] #1 = tributary, 0 = main stem
length(unique(lcgup$siteno))
sample = sample.split(lcgup$Bicarbonate, SplitRatio = .80)
train = subset(lcgup, sample == TRUE)
test  = subset(lcgup, sample == FALSE)

lcgnames <- c(lithnames,"Barren", "Open", "Forest", "CultivatedCrops", "PastureHay", "Developed")

#do rf for grouped land use
i <- 1
fit1 <- list()
imp1 <- list()
for(i in 1:length(lcgnames)){
  train1 <- cbind(train[,param], train[,lcgnames[i]]) #make dataset with all chem var and 1 lithology class
  colnames(train1) <- c(colnames(train[,param]), lcgnames[i])
  
  #train1[,lithnames[i]] <- factor(train1[,lithnames[i]]) #make response a factor - for classification
  fit1[[i]] <- randomForest(train1[,lcgnames[i]] ~.,data=train1, importance=TRUE) #run regression rf
  imp1[[i]] <- as.data.frame(importance(fit1[[i]])) #save importance
  imp1[[i]]$solutes <- rownames(imp1[[i]])
  imp1[[i]]$lc <- lcgnames[i]
  
  #use model to predict test dataset
  test1 <- cbind(test[,param], test[,lcgnames[i]]) #make dataset with all chem var and 1 lithology class
  colnames(test1) <- c(colnames(test[,param]), lcgnames[i])
  pred = predict(fit1[[i]], newdata=test1)
  
  #compare predicted to observed
  test1$pred <- pred
  fit1[[i]]$cor <- cor(test1[,lcgnames[i]], test1$pred)
}

names(fit1) <- lcgnames
names(imp1) <- lcgnames
#make summary dataframe
sum2 <- do.call(rbind, imp1)
t <- do.call(rbind,fit1)
solutes <- param
sum2$cor <- NA
for(i in 1:length(lcgnames)){
  for(j in 1:length(solutes)){
    sum2[which(sum2$lc == lcgnames[i] & sum2$solutes == solutes[j]),]$cor <- 
      cor(lcg[,solutes[j]], lcg[,lcgnames[i]])
  }
}

sum2 <- sum2[!sum2$solutes %in% sum2$lc,]

#make plot of subset RF output with all vars
sum2$lc <- factor(sum2$lc, levels=c(lcgnames))
x <- quantile(sum2$`%IncMSE`, probs=c(0.8))
sum2$num <- ifelse(sum2$`%IncMSE` > x, sum2$`%IncMSE`, NA)

plot <- ggplot(data=sum2)+
  geom_point(aes(x=lc, y=solutes, size=`%IncMSE`, color=cor))+
  labs(x="Land Use", y="Solutes", size="IncMSE", color="Corr")+
  geom_text(data=sum2, aes(x=lc, y=solutes, label=round(num, 0)), 
            color="black", size=5)+
  scale_size(range=c(3,25))+
  scale_color_distiller(type="div", palette= "RdYlBu")+
  #viridis::scale_color_viridis()+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)

#do random forest using seasonality and tributary as predictors
lcgup <- lcg
#look at seasonal variation
lcgup$month <- format(as.Date(lcgup$date, format="%Y-%m-%d"), "%m")
lcgup$season <- NA
lcgup$season <- ifelse(lcgup$month == "12"|lcgup$month == "01"|lcgup$month == "02","W",
                     ifelse(lcgup$month == "03"|lcgup$month == "04"|lcgup$month == "05","SP",
                            ifelse(lcgup$month == "06"|lcgup$month == "07"|lcgup$month == "08", "SU", "F")))

sample = sample.split(lcgup$Bicarbonate, SplitRatio = .80)
train = subset(lcgup, sample == TRUE)
test  = subset(lcgup, sample == FALSE)
i <- 1
fit1 <- list()
imp1 <- list()
for(i in 1:length(lcgnames)){
  train1 <- cbind(train[,param], train[,"season"]) #make dataset with all chem var and 1 lithology class
  colnames(train1) <- c(param, "season")
  
  #train1[,lithnames[i]] <- factor(train1[,lithnames[i]]) #make response a factor - for classification
  fit1[[i]] <- randomForest(train1[,"season"] ~.,data=train1, importance=TRUE,) #run regression rf
  imp1[[i]] <- as.data.frame(importance(fit1[[i]])) #save importance
  imp1[[i]]$solutes <- rownames(imp1[[i]])
  imp1[[i]]$lc <- "trib"
  
  #use model to predict test dataset
  test1 <- cbind(test[,param], test[,"trib"]) #make dataset with all chem var and 1 lithology class
  colnames(test1) <- c(colnames(test[,param]), "trib")
  pred = predict(fit1[[i]], newdata=test1)
  
  #compare predicted to observed
  test1$pred <- pred
  fit1[[i]]$cor <- cor(test1[,"trib"], test1$pred)
}

names(fit1) <- lcgnames
names(imp1) <- lcgnames
#make summary dataframe
sum2 <- do.call(rbind, imp1)
t <- do.call(rbind,fit1)
solutes <- param
sum2$cor <- NA
for(i in 1:length(lcgnames)){
  for(j in 1:length(solutes)){
    sum2[which(sum2$lc == "trib" & sum2$solutes == solutes[j]),]$cor <- 
      cor(lcg[,solutes[j]], lcg[,"trib"])
  }
}

sum2 <- sum2[!sum2$solutes %in% sum2$lc,]
setwd(dname)
write.csv(sum2[1:8,], "RF_tribtable.csv")
