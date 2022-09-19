#Code to pull sub watersheds from ArcGIS and summarize information
#before running this code make sure to run 4 python codes in order written below:
#WTS_Delineation
#WTS_areacalc
#AttSum

#Written by Grace Goldrich Middaugh
#Written 11/7/2021
#Last Updated by Grace Goldrich-Middaugh
#Last Updated 02/10/2022

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

#load original summary site list
setwd(dname)
#sites <- readRDS("Subwts_chem.RDS") #27 most measured sites
#sites <- readRDS("Chargebalancesites.RDS") #108 charge balanced sites with complete cases
#arc.write(file.path(paste("C:/Users/grace/Desktop/CR_SimpleWTS\\CR_SimpleWTS.gdb", "sites", sep="/")), data=sites, coords=c('long', 'lat'), shape_info=list(type='Point', WKID=4269), overwrite=TRUE)
sites <- arc.open(file.path("C:/Users/grace/Box Sync/TX/Spatial\\CR_SimpleWTS.gdb", "sites"))
sites <- arc.select(sites, fields=c('siteno', 'lat', 'long','id','trib'))
#saveRDS(sites,"Chargebalancesites.RDS")

#### Skip to line ~58 once watersheds are delineated ####
#write complete cases to GIS as individual points to delineate watersheds
setwd(dname)
wqp <- readRDS("CSwqpchargebalance.RDS")
wqp <- wqp[!is.na(wqp$siteno),] #read in charge balanced measurements
wqp$Year <- format(as.Date(wqp$date, format="%Y-%m-%d"), "%Y")
param <- c("Bicarbonate","Calcium", "Chloride", "Potassium",
           "Magnesium", "Sodium", "Sulfate", "Silica")
wqp <- wqp[complete.cases(wqp[,param]),]
length(unique(wqp$siteno)) #108 sites with complete records for param, 1101 records w/ flow, 60 sites

#make summary data to put into GIS

P <- paste("P", seq(1,108,1), sep="_")
i <- 1
for(i in 1:length(P)){
  arc.write(file.path(paste("C:/Users/grace/Desktop/CR_SimpleWTS\\CR_SimpleWTS.gdb", paste(P[i]), sep="/")), data=sites[i,], coords=c('long', 'lat'), shape_info=list(type='Point', WKID=4269), overwrite=TRUE)
}
####

#pull WSP (watershed polygons) and sum area attribute to get total area for each subwts
CRwts <- list()
WSP <- paste("WSP",(1:106)[-7], sep="_")
Lith <- paste("Lith",(1:106)[-7], sep="_")
LCov <- paste("LC",(1:106)[-7], sep="_")
ID <- (1:106)[-7]
area <- c()
i <- 1
for(i in 1:length(WSP)){
  CRwts[[i]] <- arc.open(file.path("C:/Users/grace/Desktop/CR_SimpleWTS\\CR_SimpleWTS.gdb", paste(WSP[i])))
  CRwts[[i]] <- arc.select(CRwts[[i]], fields=c('Id', 'WTS_area'))
  area[i] <- sum(CRwts[[i]]$WTS_area)
}
asum <- as.data.frame(cbind(ID,area))

#summarize lithology data
Rlist <- list()
Rtypes <- c("Sedimentary Deposit", "Conglomerate", "Carbonates", "Water", "Evaporites",
           "Mudstone", "Igneous-Metamorphic", "Sandstone")
asum[Rtypes] <- NA

for(i in 1:length(WSP)){
  Rlist[[i]] <- arc.open(file.path("C:/Users/grace/Desktop/CR_SimpleWTS\\CR_SimpleWTS.gdb", paste(Lith[i])))
  Rlist[[i]] <- arc.select(Rlist[[i]], fields=c('RockType', 'Area_SqKm'))
  
  for(j in 1:length(Rtypes)){
  asum[i,Rtypes[j]] <- sum(Rlist[[i]][Rlist[[i]]$RockType %in% Rtypes[j],]$Area_SqKm)

  }
}

names(asum)[names(asum) == "Sedimentary Deposit"] <- "Sdep"
names(asum)[names(asum) == "Igneous-Metamorphic"] <- "Igmet"
Rtypes <- c("Sdep", "Conglomerate", "Carbonates", "Water", "Evaporites",
            "Mudstone", "Igmet", "Sandstone")

#summarize land cover data
LClist <- list()
setwd(dname)
LULCtypes <- readRDS("LULCtypes.RDS") #read in vector of LULC types
asum[LULCtypes] <- NA #generate a column for each in asum df

for(i in 1:length(WSP)){
  LClist[[i]] <- arc.open(file.path("C:/Users/grace/Desktop/CR_SimpleWTS\\CR_SimpleWTS.gdb", paste(LCov[i])))
  LClist[[i]] <- arc.raster(LClist[[i]])
  LClist[[i]] <- LClist[[i]]$attribute_table()
  for(j in 1:length(LULCtypes)){
    asum[i,LULCtypes[j]] <- sum(LClist[[i]][LClist[[i]]$Land_Use %in% LULCtypes[j],]$Area_SqKm)
  }
}

names(asum) <- gsub(" ", "", names(asum))
names(asum) <- gsub("/", "", names(asum))
LULCtypes <- gsub(" ", "", LULCtypes)
LULCtypes <- gsub("/", "", LULCtypes)

#make summary dataframe as percents
psum <- asum[,c("ID", "area")]
psum[Rtypes] <- NA #add lithology columns
psum[LULCtypes] <- NA #add LULC columns

psum[Rtypes] <- asum[Rtypes]/psum$area
rowSums(psum[Rtypes]) #check all sum to 1 (or close enough)

psum[LULCtypes] <- asum[LULCtypes]/psum$area
rowSums(psum[LULCtypes])

#save as RDS
setwd(dname)
#saveRDS(psum, "CS_SimpleSum.RDS")
psum <- readRDS("CS_SimpleSum.RDS")

#join spatial data to chemistry data
setwd(dname)
CS <- readRDS("CS.data.long.RDS")
wqpc <- readRDS("CSwqpchargebalance.RDS")
wqpc$Year <- format(as.Date(wqpc$date, format="%Y-%m-%d"), "%Y")

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

#convert to compositional
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
sites <- arc.open(file.path("C:/Users/grace/Box Sync/Oregonstate/TX/Spatial/All_WTS\\All_WTS.gdb", "old_CRsites"))
sites <- arc.select(sites, fields=c("siteno", "lat", "long", "trib", "id"))
psum <- inner_join(sites, psum, by=c("id" = "ID"))

wdata <- wqcomp[wqcomp$siteno %in% psum$siteno,] #subset chem data to match 78 sites
wdata <- left_join(wdata,psum, by=c("siteno" = "siteno"))

wdata$yday <- lubridate::yday(format(as.Date(wdata$date, format="%Y-%m-%d")))
wdata$syday <- sin(wdata$yday)

##### Random Forest prediction ####

#Split into test and train datasets
#center and standardize
lithnames <- Rtypes[-4]
data <- as.data.frame(scale(wdata[,c(param, lithnames)], center=TRUE, scale = TRUE))
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
  
  #train1[,lithnames[i]] <- factor(train1[,lithnames[i]]) #make response a factor - for classification
  forest_input<-formula(paste(param[i], "~."))
  fit[[i]] <- randomForest(forest_input,data=train1, importance=TRUE) #run regression rf
  imp[[i]] <- as.data.frame(importance(fit[[i]])) #save importance
  imp[[i]]$preds <- rownames(imp[[i]])
  imp[[i]]$param <- param[i]
  
  #use model to predict test dataset
  test1 <- cbind(test[,param[i]], test[,lithnames]) #make dataset with all chem var and 1 lithology class
  colnames(test1) <- c(colnames(test[,param]), lithnames[i]))
  
  #test1[,lithnames[i]] <- factor(test1[,lithnames[i]]) #for classification
  pred = predict(fit[[i]], newdata=test1)
  
  #compare predicted to observed
  test1$pred <- pred
  fit[[i]]$cor <- cor(test1[,lithnames[i]], test1$pred)
}

names(fit) <- lithnames
names(imp) <- lithnames


#extract info to dataframe
sum <- do.call(rbind, imp)
solutes <- param
sum$cor <- NA
for(i in 1:length(lithnames)){
  for(j in 1:length(solutes)){
    sum[which(sum$lith == lithnames[i] & sum$solutes == solutes[j]),]$cor <- 
      cor(wdata[,solutes[j]], wdata[,lithnames[i]])
  }
}

sum <- sum[!sum$solutes %in% sum$lith,]

#make lithology plot
plot <- ggplot(data=sum)+
  geom_point(aes(x=lith, y=solutes, size=`%IncMSE`, color=cor))+
  labs(x="Lithology", y="Solutes", size="IncMSE", color="Corr")+
  scale_size(range=c(3,20))+
  scale_color_distiller(type="div", palette= "RdYlBu")+
  #scale_color_gradient2(low="red", mid="white", high="blue")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)

####group land use classes ####
lithnames <- Rtypes[-4]
drops <- c("OpenWater", "Water")
wdata <- wdata[ , !(names(wdata) %in% drops)]
lcg <- wdata[,c("date","siteno","id","trib","yday", "syday","NEAR",param, lithnames)]
#saveRDS(lcg, "CRlcg.RDS")
setwd(dname)
lcg <- readRDS("CRlcg.RDS")
lcg$Barren <- wdata$Barren
lcg$Open <- wdata$GrasslandHerbaceous+wdata$Shrubs+wdata$EmergentHerbaceousWetlands+wdata$WoodyWetlands
lcg$Forest <- wdata$DeciduousForest+wdata$EvergreenForest+wdata$MixedForest
lcg$CultivatedCrops <- wdata$CultivatedCrops
lcg$PastureHay <- wdata$PastureHay
lcg$Developed <- wdata$DevelopedHighIntenisy+wdata$DevelopedLowIntensity+wdata$DevelopedMediumIntensity+wdata$DevelopedOpenSpace

#split into test and train
sample = sample.split(lcg$Bicarbonate, SplitRatio = .80)
train = subset(lcg, sample == TRUE)
test  = subset(lcg, sample == FALSE)

lcgnames <- c("Barren", "Open", "Forest", "CultivatedCrops", "PastureHay", "Developed")

#do rf for grouped land use
i <- 1
fit1 <- list()
imp1 <- list()
for(i in 1:length(param)){
  train1 <- cbind(train[,lcgnames], train[,param[i]]) #make dataset with all chem var and 1 lithology class
  colnames(train1) <- c(colnames(train[,lcgnames]), param[i])
  
  #train1[,lithnames[i]] <- factor(train1[,lithnames[i]]) #make response a factor - for classification
  fit1[[i]] <- randomForest(train1[,param[i]]~.,data=train1, importance=TRUE) #run regression rf
  imp1[[i]] <- as.data.frame(importance(fit1[[i]])) #save importance
  imp1[[i]]$solutes <- rownames(imp1[[i]])
  imp1[[i]]$param <- param[i]
  
  #use model to predict test dataset
  test1 <- cbind(test[,param[i]], test[,lcgnames]) #make dataset with all chem var and 1 lithology class
  colnames(test1) <- c(colnames(test[,param[i]]), lcgnames)
  pred = predict(fit1[[i]], newdata=test1)
  
  #compare predicted to observed
  test1$pred <- pred
  fit1[[i]]$cor <- cor(test1[,lcgnames[i]], test1$pred)
}

names(fit1) <- lcgnames
names(imp1) <- lcgnames

x <- 1
cors <- c()
for(i in 1:length(fit1)){
cors[x] <- fit1[[i]]$cor
x <- x+1
}

lcgcors <- as.data.frame(cbind(lcgnames,cors))

#extract info to dataframe
sum2 <- do.call(rbind, imp1)
solutes <- param
sum2$cor <- NA
for(i in 1:length(lcgnames)){
  for(j in 1:length(solutes)){
    sum2[which(sum2$lc == lcgnames[i] & sum2$solutes == solutes[j]),]$cor <- 
      cor(lcg[,solutes[j]], lcg[,lcgnames[i]])
  }
}

sum2 <- sum2[!sum2$solutes %in% sum2$lc,]

#make plot of grouped land use
sum2$lc <- factor(sum2$lc, levels=c())
plot <- ggplot(data=sum2)+
  geom_point(aes(x=lc, y=solutes, size=`%IncMSE`, color=cor))+
  labs(x="Land Use", y="Solutes", size="IncMSE", color="Corr")+
  scale_size(range=c(3,25))+
  scale_color_distiller(type="div", palette= "RdYlBu")+
  #viridis::scale_color_viridis()+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)


#combine and plot lithology and grouped landuse with labels
sum <- sum[!(sum$lith == "Water"),]
colnames(sum2) <- colnames(sum)
sum <- rbind(sum, sum2)
Rtypes <- Rtypes[-4] #remove water
sum$lith <- factor(sum$lith, levels =c(Rtypes,lcgnames))

#identify values in top 85th percentile to label on plot
x <- quantile(sum$`%IncMSE`, probs=c(0.8))
sum$num <- ifelse(sum$`%IncMSE` > x, sum$`%IncMSE`, NA)

plot <- ggplot(data=sum)+
  geom_point(aes(x=lith, y=solutes, size=`%IncMSE`, color=cor))+
  labs(x="Watershed Factor", y="Solutes", size="%IncMSE", color="Corr")+
  scale_size(range=c(1,15))+
  geom_text(data=sum, aes(x=lith, y=solutes, label=round(num, 0)), 
            color="black", size=5)+
  scale_color_distiller(type="div", palette= "RdYlBu")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90))
print(plot)

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
