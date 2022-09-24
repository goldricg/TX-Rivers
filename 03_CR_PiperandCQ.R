#Created 11/9/2020
#Created by Grace Goldrich-Middaugh
#Last Updated 3/14/21

#scripting download and filtering of water chemistry data with arcgisbinding

#set working directories - might have to transfer later
setwd("C:\\GISdocs\\TX_scriptdat") #windows
wname <- getwd() # set back to working directory
gdbname <- paste(wname,"TX_scriptdat.gdb",sep="\\") # will open data file in working directory
boxname <- paste("C:/Users/grace/Box Sync/Oregonstate/TX")
colorado <- paste("C:/Users/grace/Box Sync/Oregonstate/TX/ColoradoSites")
dname <- paste(colorado, "Data",sep="/") # will open data file in working directory
fname<-paste(colorado,"Plots/Chem",sep="/") # will open plots file in working directory


#load packages
library(dataRetrieval)
library(magrittr)
library(dplyr)
library(ggplot2)
library(EnvStats)
library(MASS)

setwd(dname)
sites <- read.csv("WQPstations.csv")
sites <- sites[which(sites$MonitoringLocationTypeName == "River/Stream"),]
x <- sites[grep(c("USGS", "TCEQ", "TXSTRMTM"), sites$MonitoringLocationIdentifier),]
sites <- na.omit(unique(sites$MonitoringLocationIdentifier))

allloc <- readRDS("allrawCSWQPdat.RDS")
#filter to parameters of interest
CSdf <- allloc
z <- CSdf[which(CSdf$param == "Alkalinity, total"),]
CSdf <- CSdf[which(CSdf$param == "Alkalinity"| CSdf$param == "Alkalinity, bicarbonate"|CSdf$param == "Alkalinity, total"|CSdf$param == "Iron" |CSdf$param == "Bromide"| CSdf$param == "Total dissolved solids" | CSdf$param == "Phosphorus" |CSdf$param == "Chloride"|CSdf$param == "Sulfate"|CSdf$param == "Temperature, water"| CSdf$param == "Stream flow, instantaneous"|
                     CSdf$param == "Specific conductance"|CSdf$param == "pH"| CSdf$param == "Bicarbonate"| CSdf$param == "Stream flow, mean. daily"|CSdf$param == "Nitrate"|CSdf$param == "Calcium"|CSdf$param == "Magnesium"|CSdf$param == "Sodium"|
                     CSdf$param == "Silica"|CSdf$param == "Potassium"|CSdf$param == "Nitrogen"|CSdf$param == "Aluminum" |CSdf$param == "Flow"),]
CSdf$measureval <- as.numeric(CSdf$measureval)
#fix units (umho/cm = uS/cm)
unique(CSdf$units)
CSdf <- CSdf[!is.na(CSdf$units),]
CSdf <- CSdf[!(CSdf$param == "Flow" & CSdf$units == "None"),]
CSdf <- CSdf[!(CSdf$param == "Stream flow, instantaneous" | CSdf$units == "tons/day" | CSdf$units == "mg/kg"| CSdf$units == "mg/kg as P"|
                 CSdf$units == "nu"| CSdf$units == "ft/sec"| CSdf$units == "mgd"| CSdf$units == "tons/ac ft"),]
z <- CSdf[which(CSdf$units == "ft3/s"),] #4,393 measurements

for(i in 1:nrow(CSdf)){
  if(CSdf$units[i]== "ft3/s"|CSdf$units[i] == "ft3/sec"|CSdf$units[i] == "cfs"){
    CSdf$measureval[i] <- CSdf$measureval[i]/35.3147
    CSdf$units[i] <- "m3/sec"
  }else if(CSdf$param[i] == "Total dissolved solids"){
    CSdf$measureval[i] <- CSdf$measureval[i]*735.468 #from tons/ac ft to mg/l
    CSdf$units[i] <- "mg/l"
  }else if(CSdf$units[i] == "mg/l asNO3"){
    CSdf$measureval[i] <- CSdf$measureval[i]*0.2259
    CSdf$units[i] <- "mg/l as N"
  }else if(CSdf$units[i] == "ug/l"){
    CSdf$measureval[i] <- CSdf$measureval[i]/1000
    CSdf$units[i] <- "mg/l"
  }else if(CSdf$units[i] == "mg/l PO4"){
    CSdf$measureval[i] <- CSdf$measureval[i]*0.32613
    CSdf$units[i] <- "mg/l as P"
  }
}

unique(CSdf$units)
z <- CSdf[which(CSdf$param == "Total dissolved solids"),]
unique(z$units)
range(z$measureval, na.rm=T)
z <- CSdf[(CSdf$param == "Bromide"),] #9,021 measurements of total alkalinity, 2,112 Iron, 696 Aluminum, 385 Bromide
CSdf <- CSdf[!(CSdf$units == "ueq/L"),] #only 31 measurements so remove for now
###

#put into wide format
CScast <- reshape2::dcast(CSdf, date+siteno+lat+long~param, value.var="measureval",mean)
for(i in 1:nrow(CScast)){
  if(is.na(CScast$Flow[i])){
    CScast$Flow[i] <- CScast$'Stream flow, mean. daily'[i]
  }
}
CScast$'Stream flow, mean. daily' <- NULL
z <- CScast[!is.na(CScast$Flow),]
setwd(dname)
#saveRDS(CScast, "CSWQPfull.RDS")

  ###### Start HERE if you have already filtered and saved data #####
setwd(dname)
CScast <- readRDS("CSWQPfull.RDS")
length(unique(CScast$siteno))
#convert to meq/l
CScomp <- CScast
CScomp$Alkalinity <- CScomp$Alkalinity*1.22 #convert from CaCO3 to Ca(HCO3)2
param <- c("Alkalinity","Aluminum","Bicarbonate","Bromide", "Calcium", "Chloride", "Potassium","Magnesium","Nitrate","Phosphorus","Iron","Sodium","Sulfate", "Silica")
mass <- c(61.0168,26.981539,61.0168,79.904,40.078,35.453,
          39.0983,24.305,14.0067,30.9738,55.84,
          22.989769,96.06,28.0855)
charge <- c(1,3,1,1,2,1,1,2,2,1,2,1,2,1) #this will just keep P and Silica in mmol/L

CScomp <- dplyr::select(CScomp, select=c("Alkalinity","Aluminum","Bicarbonate","Bromide", "Calcium", "Chloride", "Potassium","Magnesium","Nitrate","Phosphorus","Iron","Sodium","Sulfate", "Silica", "pH", 
                                         "Temperature, water", "Total dissolved solids"))
names(CScomp) <- c(paste(param), "pH", "temp", "TDS")
CScomp <- cbind(CScomp, CScast[,c(1:4,13)])
CScomp <- reshape2::melt(CScomp, id.vars=c("date", "siteno","lat", "long", "Flow", "pH", "temp", "TDS"))
CScomp <- split(CScomp, CScomp$variable)
i <- 1
for(i in 1:length(CScomp)) {
  g_L <- CScomp[[param[i]]]$value/1000
  CScomp[[param[i]]]$mmol_L <- (g_L/mass[i])*1000
  CScomp[[param[i]]]$meq_L <- CScomp[[param[i]]]$mmol_L*charge[i]
}
CScomp <- do.call(rbind, CScomp)
CScomp1 <- reshape2::dcast(CScomp, date+siteno+lat+long+Flow+pH+temp+TDS~variable, value.var="meq_L")
z <- CScomp1[!is.na(CScomp1$Bicarbonate | CScomp1$Alkalinity),]
for(i in 1:nrow(CScomp1)){
  if(is.na(CScomp1$Bicarbonate[i])){
    CScomp1$Bicarbonate[i] <- CScomp1$Alkalinity[i]
  }
}
#calculate charge balance
for(i in 1:nrow(CScomp1)){
  pos <- sum(c(CScomp1$Calcium[i], CScomp1$Potassium[i], CScomp1$Magnesium[i], CScomp1$Sodium[i]), na.rm=T)
  neg <- sum(c(CScomp1$Chloride[i],CScomp1$Sulfate[i], CScomp1$Bicarbonate[i], CScomp1$Nitrate[i]), na.rm=T)
  CScomp1$CBerror[i] <- ((pos-neg)/(pos+neg))*100
}
mean(abs(CScomp1$CBerror), na.rm=T)
#only keep records with <10% error
CScomp1 <- CScomp1[!(abs(CScomp1$CBerror) > 10),]
mean(abs(CScomp1$CBerror), na.rm=T) #12,413 charge balanced samples
setwd(dname)
#saveRDS(CScomp1, "CSwqpchargebalance.RDS")

CScomp1 <- readRDS("CSwqpchargebalance.RDS")
CScomp1 <- CScomp1[!is.na(CScomp1$siteno),]
length(unique(CScomp1$siteno)) #155 sites
z <- CScomp1[!is.na(CScomp1$temp),] #7,568 measurements, 6.852 measurements with pH, 3,868 with temp

#make schoeller plot of charge balanced samples
#add NEAR FID
CS <- readRDS("CS.data.long.RDS")
CSunique <- CS %>% 
  group_by(siteno, NEAR_FID) %>% 
  filter(row_number() == 1)
CSunique <- na.omit(as.data.frame(CSunique[,c(2,9)]))
CScomp1$NEAR <- NA
for(i in 1:nrow(CScomp1)){
  for(j in 1:nrow(CSunique)){
    if(CScomp1$siteno[i] == CSunique$siteno[j]){
      CScomp1$NEAR[i] <- CSunique$NEAR_FID[j]
    }
  }
}

CSunique <- na.omit(CSunique)
CScomp1 <- CScomp1[!is.na(CScomp1$siteno),]
for(i in 1:nrow(CScomp1)){
  for(j in 1:nrow(CSunique)){
    if(CScomp1$siteno[i] == CSunique$siteno[j]){
      CScomp1$NEAR[i] <- CSunique$NEAR_FID[j]
    }
  }
}


meq.sch1 <- CScomp1[,c(1,2,8:21,23)]
meq.sch1 <- reshape2::melt(meq.sch1,id.vars=c("date","siteno","NEAR"))
meq.sch1$year <- format(as.Date(meq.sch1$date, format="%Y-%m-%d"), "%Y")

setwd(fname)
pdf("WQP_schoeller_3_14_21.pdf")
plot <- ggplot(data=meq.sch1,aes(x=variable,y=value, color=as.numeric(year)))+
  geom_hline(yintercept=1, color="red")+
  geom_point()+
  viridis::scale_color_viridis()+
  geom_line(aes(group=siteno))+
  scale_y_log10(labels=function(x) round(x,2)) +
  labs(x="Parameter", y="Concentration (meq/L)", title="Colorado River")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90))+
  facet_wrap(~NEAR)
print(plot)
dev.off()
#plot just site 8 as example
m <- meq.sch1[(meq.sch1$NEAR == "8"),]
pdf("WQPsite8_schoeller_3_14_21.pdf")
plot <- ggplot(data=meq.sch1,aes(x=variable,y=value, color=as.numeric(year)))+
  geom_hline(yintercept=1, color="red")+
  geom_point()+
  viridis::scale_color_viridis()+
  geom_line(aes(group=siteno))+
  scale_y_log10(labels=function(x) round(x,2)) +
  labs(x="Parameter", y="Concentration (meq/L)", title="Colorado River Site 8")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90))
print(plot)
dev.off()

#make piper diagram
toPercent <- function (d) {
  totalCations <- d$Calcium + d$Magnesium + d$Sodium + d$Potassium
  d$Calcium <- 100 * (d$Calcium/totalCations)
  d$Magnesium <- 100 * (d$Magnesium/totalCations)
  d$Sodium <- 100 * (d$Sodium/totalCations)
  d$Potassium <- 100 * (d$Potassium/totalCations)
  totalAnions <- d$Chloride + d$Sulfate + d$Carbonate + d$Bicarbonate
  d$Chloride <- 100 * (d$Chloride/totalAnions)
  d$Sulfate <- 100 * (d$Sulfate/totalAnions)
  d$Carbonate <- 100 * (d$Carbonate/totalAnions)
  d$Bicarbonate <- 100 * (d$Bicarbonate/totalAnions)
  return(d)
}

setwd(dname)
CScomp1 <- readRDS("CRlcg.RDS")
param <- c("Bicarbonate","Calcium", "Chloride", "Potassium", "Magnesium", "Sodium", "Sulfate")

pip <- CScomp1[complete.cases(CScomp1[,c("date","siteno","NEAR", "trib",param)]),]
pip$Carbonate <- 0
pip <- as.data.frame(cbind(toPercent(pip[,c("Carbonate", param)]),pip[,c("date", "siteno", "NEAR", "trib")]))
#check that add to 100
cation.sums <- apply(pip[, c("Calcium", "Magnesium", "Sodium", "Potassium")], 1, FUN = sum)
anion.sums  <- apply(pip[, c("Chloride", "Sulfate", "Carbonate", "Bicarbonate")], 1, FUN = sum)

#setwd(dname)
#saveRDS(pip, "WQPpiper.RDS")
transform_piper_data <- function(Mg, Ca, Cl,SO4, name=NULL){
  if(is.null(name)){
    name = rep(1:length(Mg),3)
  } else {
    name = rep(name,3)
  }
  y1 <- Mg * 0.86603
  x1 <- 100*(1-(Ca/100) - (Mg/200))
  y2 <- SO4 * 0.86603
  x2 <-120+(100*Cl/100 + 0.5 * 100*SO4/100)
  new_point <- function(x1, x2, y1, y2, grad=1.73206){
    b1 <- y1-(grad*x1)
    b2 <- y2-(-grad*x2)
    M <- matrix(c(grad, -grad, -1,-1), ncol=2)
    intercepts <- as.matrix(c(b1,b2))
    t_mat <- -solve(M) %*% intercepts
    data.frame(x=t_mat[1,1], y=t_mat[2,1])
  }
  np_list <- lapply(1:length(x1), function(i) new_point(x1[i], x2[i], y1[i], y2[i]))
  npoints <- do.call("rbind",np_list)
  data.frame(observation=name,x=c(x1, x2, npoints$x), y=c(y=y1, y2, npoints$y))
}

piper_data <- transform_piper_data(Ca   = pip$Calcium,
                                   Mg   = pip$Magnesium,
                                   Cl   = pip$Chloride,
                                   SO4  = pip$Sulfate,
                                   name = pip$siteno)
piper_data <- merge(piper_data,
                    pip[,c("date","NEAR", "siteno", "trib")],
                    by.x = "observation",
                    by.y = "siteno")

###define ggplot_piper function ####
ggplot_piper <- function() {
  library(ggplot2)
  grid1p1 <<- data.frame(x1 = c(20,40,60,80), x2= c(10,20,30,40),y1 = c(0,0,0,0), y2 = c(17.3206,34.6412,51.9618, 69.2824))
  grid1p2 <<- data.frame(x1 = c(20,40,60,80), x2= c(60,70,80,90),y1 = c(0,0,0,0), y2 = c(69.2824, 51.9618,34.6412,17.3206))
  grid1p3 <<- data.frame(x1 = c(10,20,30,40), x2= c(90,80,70,60),y1 = c(17.3206,34.6412,51.9618, 69.2824), y2 = c(17.3206,34.6412,51.9618, 69.2824))
  grid2p1 <<- grid1p1
  grid2p1$x1 <- grid2p1$x1+120
  grid2p1$x2 <- grid2p1$x2+120
  grid2p2 <<- grid1p2
  grid2p2$x1 <- grid2p2$x1+120
  grid2p2$x2 <- grid2p2$x2+120
  grid2p3 <<- grid1p3
  grid2p3$x1 <- grid2p3$x1+120
  grid2p3$x2 <- grid2p3$x2+120
  grid3p1 <<- data.frame(x1=c(100,90, 80, 70),y1=c(34.6412, 51.9618, 69.2824, 86.603), x2=c(150, 140, 130, 120), y2=c(121.2442,138.5648,155.8854,173.2060))
  grid3p2 <<- data.frame(x1=c(70, 80, 90, 100),y1=c(121.2442,138.5648,155.8854,173.2060), x2=c(120, 130, 140, 150), y2=c(34.6412, 51.9618, 69.2824, 86.603))
  
  p <- ggplot() +
    ## left hand ternary plot
    geom_segment(aes(x=0,y=0, xend=100, yend=0)) +
    geom_segment(aes(x=0,y=0, xend=50, yend=86.603)) +
    geom_segment(aes(x=50,y=86.603, xend=100, yend=0)) +
    ## right hand ternary plot
    geom_segment(aes(x=120,y=0, xend=220, yend=0)) +
    geom_segment(aes(x=120,y=0, xend=170, yend=86.603)) +
    geom_segment(aes(x=170,y=86.603, xend=220, yend=0)) +
    ## Upper diamond
    geom_segment(aes(x=110,y=190.5266, xend=60, yend=103.9236)) +
    geom_segment(aes(x=110,y=190.5266, xend=160, yend=103.9236)) +
    geom_segment(aes(x=110,y=17.3206, xend=160, yend=103.9236)) +
    geom_segment(aes(x=110,y=17.3206, xend=60, yend=103.9236)) +
    ## Add grid lines to the plots
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid1p1, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid1p2, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid1p3, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid2p1, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid2p2, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid2p3, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid3p1, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid3p2, linetype = "dashed", size = 0.25, colour = "grey50") +
    ### Labels and grid values
    #geom_text(aes(50,-10, label="Ca^2"), parse=T, size=4) + # Commented out, as parse=TRUE can cause issues
    
    geom_text(aes(c(20,40,60,80),c(-5,-5,-5,-5), label=c(80, 60, 40, 20)), size=3) +
    geom_text(aes(c(35,25,15,5),grid1p2$y2, label=c(80, 60, 40, 20)), size=3) +
    geom_text(aes(c(95,85,75,65),grid1p3$y2, label=c(80, 60, 40, 20)), size=3) +
    # geom_text(aes(17,50, label="Mg^2"), parse=T, angle=60, size=4) +
    coord_equal(ratio=1)+  
    geom_text(aes(17,50, label="Mg^2"), angle=60, size=4, parse=TRUE) +  
    geom_text(aes(82.5,50, label="Na + K"), angle=-60, size=4) +
    geom_text(aes(50,-10, label="Ca^2"), size=4, parse=TRUE) +
    
    
    geom_text(aes(170,-10, label="Cl^-phantom()"), size=4, parse=TRUE) +
    geom_text(aes(205,50, label="SO^4"), angle=-60, size=4, parse=TRUE) +
    geom_text(aes(137.5,50, label="Alkalinity~as~HCO^3"), angle=60, size=4, parse=TRUE) +
    geom_text(aes(72.5,150, label="SO^4~+~Cl^-phantom()"), angle=60, size=4, parse=TRUE) +
    geom_text(aes(147.5,150, label="Ca^2~+~Mg^2"), angle=-60, size=4, parse=TRUE) + 
    
    geom_text(aes(c(155,145,135,125),grid2p2$y2, label=c(20, 40, 60, 80)), size=3) +
    geom_text(aes(c(215,205,195,185),grid2p3$y2, label=c(20, 40, 60, 80)), size=3) +
    geom_text(aes(c(140,160,180,200),c(-5,-5,-5,-5), label=c(20, 40, 60, 80)), size=3) +
    geom_text(aes(grid3p1$x1-5,grid3p1$y1, label=c(80, 60, 40, 20)), size=3) +
    geom_text(aes(grid3p1$x2+5,grid3p1$y2, label=c(20, 40, 60, 80)), size=3) +
    geom_text(aes(grid3p2$x1-5,grid3p2$y1, label=c(20, 40, 60, 80)), size=3) +
    geom_text(aes(grid3p2$x2+5,grid3p2$y2, label=c(80, 60, 40, 20)), size=3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  return(p)
}

##end ##
unique(piper_data$NEAR)
piper_data$NEAR <- factor(piper_data$NEAR)
piper_data$triblab <- ifelse(piper_data$trib == 1, "Tributary", "Main Stem")
piper_data$month <- format(as.Date(piper_data$date, format="%Y-%m-%d"), "%m")
piper_data$season <- NA
piper_data$season <- ifelse(piper_data$month == "12"|piper_data$month == "01"|piper_data$month == "02","W",
                     ifelse(piper_data$month == "03"|piper_data$month == "04"|piper_data$month == "05","SP",
                            ifelse(piper_data$month == "06"|piper_data$month == "07"|piper_data$month == "08", "SU", "F")))

piper_data1 <- piper_data[which(piper_data$season == "SU"|piper_data$season == "W"),]

setwd(fname)
jpeg("piperszn_SUW.jpeg",width=1000, height=800, quality=100)
ggplot_piper() +
  geom_point(data=piper_data1, mapping = aes(x,y, color=season), size=4.5, alpha=0.25) +
  #scale_colour_manual(name="Zone", values=c("2"="#440154FF","3"="#481F70FF","4"="#443A83FF",
   #                                     "5"="#3B528BFF","6"="#31688EFF","7"="#287C8EFF",
    #                                    "8"="#21908CFF","9"="#20A486FF","10"="#35B779FF",
     #                                   "12"="#5DC863FF","13"="#8FD744FF","14"="#C7E020FF"),
      #                aesthetics = c("color", "fill")) +
  #scale_shape_manual(name="Position",values=c("Main Stem"=16, "Tributary"=3))+
  labs(color="Season", shape="Position")+
  theme(legend.position = c(1, 0.5), 
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size=12)))
dev.off()

setwd(dname)
CScomp1 <- readRDS("CSwqpchargebalance.RDS")


####C-Q relationships ####
CScomp1 <- CScomp1[!is.na(CScomp1$Flow),] #remove records with no flow, 3,607 complete samples

modlist <- list() #select sites with > 4 measurements
x <- 1
z <- split(CScomp1, CScomp1$siteno)
for(i in 1:length(z)){
  if(nrow(z[[i]]) > 4){
   modlist[[x]] <- z[[i]] 
   x <- x+1
  }
}

#calculate stats for number of obs
num_rows<-lapply(modlist, nrow)
num_rows_vect<-unlist(num_rows)
max(num_rows_vect)


#fit linear models by individual sites
library(MASS)
z <- do.call(rbind,modlist)
lmod <- list()
bc <- c()
slope <- c()
rsq <- c()
radj <- c()
pval <- c()
statsite <- c()
locsite <- c()
statparam <- c()
SE <- c()
lat <- c()
long <- c()
CVc <- c()
CVq <- c()
CVratio <- c()
plot <- list()
x <- 1
i <- 1
j <- 1
param <- c("Bicarbonate", "Calcium", "Chloride", "Potassium","Magnesium","Phosphorus","Nitrate","Sodium","Sulfate", "Silica")
setwd(paste(colorado, "Plots/Chem", sep='/'))
pdf("WQP_cqresiduals11_21_21.pdf")
for(i in 1:length(modlist)) {
  for(j in 1:length(param)){
    data <- modlist[[i]][paste(param[j])]
    data <- cbind(data, modlist[[i]]$Flow)
    data[,1] <- log(data[,1])
    data[,2] <- log(data[,2])
    data <- data[is.finite(data[,1]),] #remove -Inf values
    data <- data[is.finite(data[,2]),]
    data <- data[!is.na(data[,1]),] #remove remaining NA rows
    data <- data[!is.na(data[,2]),]
    if(nrow(data) >4){
    lmod[[x]] <- lm(data[,1]~data[,2], data=data) #run linear model
    #plot residuals 
    Obs <- seq(1:nrow(data))
    plot[[x]] <- ggplot(data=lmod[[x]], aes(x=Obs, y=lmod[[x]]$residuals))+
      geom_point()+
      geom_hline(yintercept=0, color="red")+
      labs(x="Observation", y="Residual", title=paste(param[j], modlist[[i]]$NEAR[1], sep=" Site "))
    print(plot[[x]])
    
    lsum <- base::summary(lmod[[x]])
    slope[x] <- lsum$coefficients[2]
    SE[x] <- lsum$coefficients[4]
    pval[x] <- lsum$coefficients[8]
    rsq[x] <- lsum$r.squared
    radj[x] <- lsum$adj.r.squared
    locsite[x] <- modlist[[i]]$siteno[1]
    lat[x] <- modlist[[i]]$lat[1]
    long[x] <- modlist[[i]]$long[1]
    statparam[x] <- paste(param[j])
    x <- x+1
  }
  }
}
dev.off()
statsum1 <- data.frame(locsite,lat,long, statparam, slope, rsq, radj, SE, pval)

#add NEAR
x <- left_join(statsum1, CSunique, by=c("locsite"="siteno"))
setwd(dname)
#saveRDS(statsum1, "WQPslopeSummary.RDS")
statsum1 <- readRDS("WQPslopeSummary.RDS")
#plot these values
#plot slopes
statsum2$statsite <- factor(statsum2$statsite, levels=c('2','3','4','5','6','7','8','9','12','13','14','15'))
statsum1 <- within(statsum1, color <- "not significant") #color = R if p>0.05
statsum1[!(statsum1$pval > 0.05), "color"] <- "significant"
statsum1$statparam <- factor(statsum1$statparam, levels=c("Bicarbonate", "Calcium", "Chloride", "Magnesium", "Sulfate", "Sodium","Potassium", "Silica", "Nitrate", "Phosphorus"))
unique(statsum1$statparam)

setwd(fname)
pdf("PosterC_Q.pdf")
plot <- ggplot(data=statsum1)+
  #geom_errorbar(aes(x=statsite,y=slope,ymin=low, ymax=high), width=.2, size=0.5)+
  geom_point(aes(x=long, y=slope,color=rsq, shape=color), size=5)+
  geom_smooth(aes(x=as.numeric(long),y=slope), alpha=0, color="black")+
  ylim(-0.8,0.5)+
  scale_color_distiller(type="div", palette= "RdYlBu")+
  scale_shape_manual(name="Significance", values=c("significant"=19, "not significant"=3))+
  labs(x="Longitude (upstream to downstream)", y="Slope (log C vs log Q)", shape="Significance", color=expression(R^2))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), legend.position = c(0.76,0.14),
        strip.text=element_text(size=16), axis.title=element_text(size=16))+
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1))+
  facet_wrap(~statparam,scales="fixed")
print(plot)

CVdat <- statsum1
for(i in 1:nrow(CVdat)){
  if(CVdat$statparam[i] == "Bicarbonate"){
    CVdat$category[i] <- paste("Bicarbonate")
  }else if(CVdat$statparam[i] == "Calcium"|CVdat$statparam[i] == "Magnesium"|CVdat$statparam[i] == "Sodium"){
    CVdat$category[i] <- paste("Cation")
  }else if(CVdat$statparam[i] == "Chloride"|CVdat$statparam[i]=="Sulfate"){
    CVdat$category[i] <- paste("Anion")
  }else if(CVdat$statparam[i]=="Potassium"|CVdat$statparam[i] == "Nitrate"|CVdat$statparam[i]=="Phosphorus"){
    CVdat$category[i] <- paste("Nutrient")
  }
}
plot2 <- ggplot(data=CVdat,aes(x=exp(CVratio), y=slope, color=statparam, fill=statparam))+
  geom_point(size=2)+
  #geom_convexhull(alpha=0.5)+
  stat_ellipse(aes(x=exp(CVratio), y=slope,color=statparam),type = "norm")+
  ylim(-1,1)+
  xlim(-1.1,4.1)+
  geom_segment(x=0,y=0, xend=1,yend=1,lwd=0.7, color="black")+
  geom_segment(x=0,y=0,xend=1,yend=-1,lwd=0.7, color="black")+
  geom_vline(xintercept=1, color="grey60")+
  geom_hline(yintercept=0, color="grey60")+
  #viridis::scale_color_viridis(discrete=T)+
  labs(x="CVc/CVq", y="Slope (log C vs log Q)")+
  theme(legend.title=element_blank())+
  facet_wrap(~category,scales="fixed")
print(plot2)
dev.off()
