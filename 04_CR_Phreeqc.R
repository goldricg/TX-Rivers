#Created 2/20/21
#Created by Grace Goldrich-Middaugh
#Last Updated 5/24/21
#Last updated by Grace Goldrich-Middaugh

#code to use the phreeqc package for TX surface and groundwater chemical data

#load packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(hydrogeo)
library(viridis)
library(reshape2)
library(phreeqc)

#set working directories
setwd("C:/Users/grace/Box Sync/Oregonstate/TX/ColoradoSites") #windows
wname <- getwd() # set back to working directory
dname <- paste(wname,"Data",sep="/") # will open data file in working directory
Oname <- paste(wname,"Outputs",sep="/") # will open outputs file in working directory
fname<-paste(wname,"Plots/Chem",sep="/") # will open plots file in working directory
xname<-paste(wname,"Functions",sep="/") # will open functions file in working directory

#load data
setwd(dname)
wqpc <- readRDS("CSwqpchargebalance.RDS") #charge balanced samples in meq/L
CS <- readRDS("CS.data.long.RDS") #data with near sites labelled

#add NEAR_FID to wqp data
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

#convert to mmol/L
param <- c("Alkalinity","Aluminum","Bicarbonate","Bromide", "Calcium", "Chloride", "Potassium","Magnesium","Nitrate","Phosphorus","Iron","Sodium","Sulfate", "Silica")
charge <- c(1,3,1,1,2,1,1,2,2,1,2,1,2,1) #this will just keep P and Silica in mmol/L
for(i in 1:length(param)){
  wqpc[,paste(param[i])] <- wqpc[,paste(param[i])]/charge[i]
}

#subset to only sites with all params
wqpc <- dplyr::select(wqpc, select=-c(Alkalinity))

#set up dataframe with complete conc measurements
param <- c("Bicarbonate","Calcium", "Chloride", "Potassium","Magnesium","Sodium","Sulfate", "Silica")

si.dat <- wqpc[complete.cases(wqpc[,c(param)]),] #5282 measurements without iron or aluminum

#look for pH outliers
plot <- ggplot(data=si.dat)+
  geom_point(aes(x=long, y=pH))
print(plot)
#remove 1 pH point of 42.85
si.dat <- si.dat[which(si.dat$pH < 20),]
z <- si.dat[max(si.dat$Potassium),]

#subset to observations with temperature measurements
si.dat <- si.dat[!is.na(si.dat$temp),] #3868 observations
mean(si.dat$temp)

#set up si columns
si.dat$aragonite <- NA
si.dat$anhydrite <- NA
si.dat$calcite <- NA
si.dat$chalcedony <- NA
si.dat$dolomite <- NA
si.dat$gypsum <- NA
si.dat$quartz <- NA
si.dat$halite <- NA
si.dat$kmica <- NA
si.dat$kfeld <- NA
si.dat$illite <- NA

# load the phreeqc.dat database
phrLoadDatabaseString(phreeqc.dat)
for(i in 1:nrow(si.dat)){
# accumulate the input
phrAccumulateLine("TITLE SI test")
phrAccumulateLine("SOLUTION 1")
phrAccumulateLine(paste("pH ", si.dat[i,"pH"]))
phrAccumulateLine(paste(" temp ", si.dat[i,"temp"]))
phrAccumulateLine(" units mmol/l")
phrAccumulateLine(paste("Alkalinity", si.dat[i,"Bicarbonate"]))
phrAccumulateLine(paste("Ca", si.dat[i,"Calcium"]))
phrAccumulateLine(paste("Cl", si.dat[i,"Chloride"]))
phrAccumulateLine(paste("K", si.dat[i,"Potassium"]))
phrAccumulateLine(paste("Mg", si.dat[i,"Magnesium"]))
#phrAccumulateLine(paste("Fe", si.dat[i,"Iron"]))
#phrAccumulateLine(paste("Al", si.dat[i,"Aluminum"]))
phrAccumulateLine(paste("Na", si.dat[i,"Sodium"]))
phrAccumulateLine(paste("Si", si.dat[i,"Silica"]))
phrAccumulateLine(paste("S(6)", si.dat[i,"Sulfate"]))
phrAccumulateLine("SELECTED_OUTPUT")
phrAccumulateLine(" -file test.sel")
phrAccumulateLine(" -si aragonite calcite chalcedony dolomite quartz halite gypsum")
phrAccumulateLine("END")
# run it and echo the name of the output file
if (is.null(phrRunAccumulated())) {
  cat(paste("see ", phrGetOutputFileName(), ".\n", sep = ""))
}

x <- phrGetSelectedOutput()
si.dat$aragonite[i] <- x[[1]]$si_aragonite[1]
#si.dat$anhydrite[i] <- x[[1]]$si_anhydrite[1]
si.dat$calcite[i] <- x[[1]]$si_calcite[1]
si.dat$chalcedony[i] <- x[[1]]$si_chalcedony[1]
si.dat$dolomite[i] <- x[[1]]$si_dolomite[1]
si.dat$gypsum[i] <- x[[1]]$si_gypsum[1]
si.dat$quartz[i] <- x[[1]]$si_quartz[1]
si.dat$halite[i] <- x[[1]]$si_halite[1]
#si.dat$kmica[i] <- x[[1]]$si_k.mica[1]
#si.dat$kfeld[i] <- x[[1]]$si_k.feldspar[1]
#si.dat$illite[i] <- x[[1]]$si_illite[1]
}

#put in long format to plot
long.si <- reshape2::melt(si.dat, id.vars=c("siteno", "lat", "long", "date", "NEAR","Flow", "temp"))
long.si <- long.si[which(long.si$variable == "aragonite"|long.si$variable == "anhydrite"|long.si$variable == "calcite"|
                           long.si$variable == "chalcedony"|long.si$variable == "gypsum"|long.si$variable == "kmica"|long.si$variable == "kfeld"|
                           long.si$variable == "quartz"|long.si$variable == "dolomite"|long.si$variable == "halite"|long.si$variable == "illite"),]

plot <- ggplot(data=long.si, aes(x=value,y=log10(Flow)))+
  geom_point(aes(color=variable))
print(plot)

#plot each mineral by flow and color by site
z <- long.si[which(long.si$variable == "calcite"|long.si$variable == "gypsum"|long.si$variable == "aragonite"|
                     long.si$variable == "chalcedony"|long.si$variable == "halite"|long.si$variable == "quartz"),]
z <- z[!is.na(z$NEAR),]
levels(z$variable)
z$NEAR <- factor(z$NEAR)
mins <- as_labeller(c("calcite" = "Calcite","aragonite" = "Aragonite", "gypsum" = "Gypsum", "chalcedony" = "Chalcedony", "halite"="Halite"))
#plot SI vs flow
plot <- ggplot(data=z, aes(x=value, y=log10(Flow)))+
  geom_point(aes(color=NEAR), size=3)+
  scale_color_brewer(type="div", palette= "RdYlBu")+
  labs(x="SI", y="log(Flow)", color="Zone")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~variable, scales="fixed", labeller=mins) #+
  #theme(legend.position = "none")
print(plot)

#plot just quartz versus flow
z <- long.si[which(long.si$variable == "quartz"),]
plot <- ggplot(data=z, aes(x=value, y=log10(Flow)))+
  geom_point(aes(color=NEAR), size=3)
print(plot)

#plot SI vs long
plot <- ggplot(data=long.si)+
  geom_point(aes(x=long, y=value, group=siteno, color=variable))+
  geom_hline(yintercept=1)+
  ylim(-4,2)
print(plot)

#make boxplots with each mineral group
z <- long.si[which(long.si$variable == "aragonite"|long.si$variable == "gypsum"|
                     long.si$variable == "chalcedony"|long.si$variable == "halite"),]
z <- z[!is.na(z$NEAR),]
z <- z[!is.na(z$variable),]
mins <- as_labeller(c("aragonite" = "Aragonite", "gypsum" = "Gypsum", "chalcedony" = "Chalcedony", "halite"="Halite"))

plot <- ggplot(data=z, aes(x=factor(NEAR), y=value, fill=variable))+
  geom_boxplot(position=position_dodge(width=0.8))+
  geom_hline(yintercept=1)+
  #stat_summary(fun.data = give.n, geom = "text", fun = median,
   #            position = position_dodge(width = 0.75))+
  labs(x="Zone", y="Calculated Saturation Index", fill="Mineral")+
  #scale_fill_manual(values=c(RColorBrewer::brewer.pal(4,"Spectral")), 
   #                  labeller=mins)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plot)

i <- 1
p <- ggplot()
for(i in 1:length(mins)){
  z <- long.si[which(long.si$variable == mins[i]),]
  p<-p+geom_boxplot(data=z, aes(x=factor(NEAR), y=value, fill=variable),
               position=position_dodge(width=0.8))
  p<-p+geom_text(data=n.sum[(n.sum$variable == mins[i]),], aes(label=value,x=NEAR, y=median(as.numeric(long.si[which(long.si$variable == mins[i]),]$value))-0.5), 
            size=6, position=position_dodge(width=0.8))
}
p<-p+labs(x="Zone", y="SI", fill="Mineral")+
  theme(axis.text.x = element_text(angle = 90))+
  theme_bw()
print(p)


#plot SI vs season for aragonite
main <- long.si[which(long.si$variable == "aragonite"),]
main$month <- format(as.Date(main$date, format="%Y-%m-%d"), "%m")
main$season <- NA
main$season <- ifelse(main$month == "12"|main$month == "01"|main$month == "02","W",
                     ifelse(main$month == "03"|main$month == "04"|main$month == "05","SP",
                            ifelse(main$month == "06"|main$month == "07"|main$month == "08", "SU", "F")))
main <- main[!is.na(main$temp),]
plot <- ggplot(data=main, aes(y=log10(Flow), x=value))+
  geom_point(aes(color=temp),size=4)+
  viridis::scale_color_viridis()
print(plot)

#### Calculate for temp of 10C ####
si.dat1 <- wqpc[complete.cases(wqpc[,c(param)]),] #5282 measurements without iron or aluminum
si.dat1 <- si.dat1[!is.na(si.dat1$temp),] #3868 observations

si.dat1 <- si.dat1[which(si.dat1$pH < 20),] #remove 1 pH point of 42.85
y <- si.dat[!is.na(si.dat$temp),]
mean(y$temp)
#set up si columns
si.dat1$aragonite <- NA
si.dat1$chalcedony <- NA
si.dat1$gypsum <- NA
si.dat1$halite <- NA

# load the phreeqc.dat database
phrLoadDatabaseString(phreeqc.dat)
for(i in 1:nrow(si.dat1)){
  # accumulate the input
  phrAccumulateLine("TITLE SI test")
  phrAccumulateLine("SOLUTION 1")
  phrAccumulateLine(paste("pH ", si.dat1[i,"pH"]))
  phrAccumulateLine(" temp 25.0")
  phrAccumulateLine(" units mmol/l")
  phrAccumulateLine(paste("Alkalinity", si.dat1[i,"Bicarbonate"]))
  phrAccumulateLine(paste("Ca", si.dat1[i,"Calcium"]))
  phrAccumulateLine(paste("Cl", si.dat1[i,"Chloride"]))
  phrAccumulateLine(paste("K", si.dat1[i,"Potassium"]))
  phrAccumulateLine(paste("Mg", si.dat1[i,"Magnesium"]))
  phrAccumulateLine(paste("Na", si.dat1[i,"Sodium"]))
  phrAccumulateLine(paste("Si", si.dat1[i,"Silica"]))
  phrAccumulateLine(paste("S(6)", si.dat1[i,"Sulfate"]))
  phrAccumulateLine("SELECTED_OUTPUT")
  phrAccumulateLine(" -file test.sel")
  phrAccumulateLine(" -si aragonite calcite chalcedony dolomite quartz halite gypsum")
  phrAccumulateLine("END")
  # run it and echo the name of the output file
  if (is.null(phrRunAccumulated())) {
    cat(paste("see ", phrGetOutputFileName(), ".\n", sep = ""))
  }
  
  x <- phrGetSelectedOutput()
  si.dat1$aragonite[i] <- x[[1]]$si_aragonite[1]
  si.dat1$chalcedony[i] <- x[[1]]$si_chalcedony[1]
  si.dat1$gypsum[i] <- x[[1]]$si_gypsum[1]
  si.dat1$halite[i] <- x[[1]]$si_halite[1]
}

#put in long format to plot
long.si <- reshape2::melt(si.dat1, id.vars=c("siteno", "lat", "long", "date", "NEAR","Flow", "temp"))

z1 <- long.si[which(long.si$variable == "aragonite"|long.si$variable == "gypsum"|
                     long.si$variable == "chalcedony"|long.si$variable == "halite"),]
z1 <- z1[!is.na(z1$NEAR),]
z1 <- z1[!is.na(z1$variable),]
mins <- as_labeller(c("aragonite" = "Aragonite", "gypsum" = "Gypsum", "chalcedony" = "Chalcedony", "halite"="Halite"))

plot <- ggplot(data=z1, aes(x=factor(NEAR), y=value, fill=variable))+
  geom_boxplot(position=position_dodge(width=0.8))+
  geom_hline(yintercept=1)+
  #stat_summary(fun.data = give.n, geom = "text", fun = median,
  #            position = position_dodge(width = 0.75))+
  labs(x="Zone", y="Calculated Saturation Index", fill="Mineral")+
  #scale_fill_manual(values=c(RColorBrewer::brewer.pal(4,"Spectral")), 
  #                  labeller=mins)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plot)

z$pchange <- ((z$value-z1$value)/z$value)*100

