#Test PCA for Colorado river based on PCA_example.R code in R_info folder
#transformation of water quality portal prefiltered data
#Grace Goldrich-Middaugh
#10-13-21

#set working directories
setwd("C:/Users/grace/Box Sync/Oregonstate/TX/ColoradoSites") #windows
wname <- getwd() # set back to working directory
dname <- paste(wname,"Data",sep="/") # will open data file in working directory
Oname <- paste(wname,"Outputs",sep="/") # will open outputs file in working directory
fname<-paste(wname,"Plots/PCA",sep="/") # will open plots file in working directory
xname<-paste(wname,"Functions",sep="/") # will open functions file in working directory

#colorado surface chem PCA
#install_github("vqv/ggbiplot")
#library('ggbiplot')
library('ggplot2')
library('tidyr')
library('WVPlots') #devtools::install_github('WinVector/WVPlots',build_vignettes=TRUE)
barbell_plot = function(frame, xvar, ymin, ymax, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar))
  } else {
    gplot = ggplot(frame, aes_string(x=xvar, color=colorvar))
  }
  
  gplot + geom_point(aes_string(y=ymin)) + 
    geom_point(aes_string(y=ymax)) +
    geom_linerange(aes_string(ymin=ymin, ymax=ymax)) +
    ylab("value")
}
dotplot_identity = function(frame, xvar, yvar, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar))
  } else {
    gplot = ggplot(frame, 
                   aes_string(x=xvar, y=yvar, ymax=yvar, 
                              color=colorvar))
  }
  gplot + geom_point() + geom_linerange(aes(ymin=0))
}
extractProjection <- function(ndim,princ) {
  # pull off the rotation.  
  proj <- princ$rotation[,1:ndim] 
  # sign was arbitrary, so flip in convenient form
  for(i in seq_len(ndim)) {
    si <- sign(mean(proj[,i]))
    if(si!=0) {
      proj[,i] <- proj[,i]*si
    }
  }
  proj
}
rsq <- function(x,y) {
  1 - sum((y-x)^2)/sum((y-mean(y))^2)
}

library(ggfortify)
library(factoextra)
library('vtreat')
library(dplyr)
library(arcgisbinding)
arc.check_product()

# load Colorado all chem, gw chem, and precipitation data
setwd("C:/Users/grace/Box Sync/Oregonstate/TX/ColoradoSites/Data")
CS <- readRDS("CS.data.long.RDS")
wqpc <- readRDS("CSwqpchargebalance.RDS")
length(unique(wqpc$siteno))
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

# format data with columns = chem parameters, rows = observations with complete cases
z <- split(CS, CS$param) #P, K, and Bco have the fewest measurements, rest close to 200 obs each, pH has 365
wqp <- na.omit(reshape2::melt(wqpc, id.vars=c("date", "siteno", "lat", "long", "CBerror", "NEAR")))
wqp <- split(wqp, wqp$variable) #remove Phosphorus, Alkalinity
wqp <- do.call(rbind,wqp)
wqp <- wqp[!(wqp$variable == "Phosphorus"|wqp$variable == "Alkalinity"|wqp$variable == "Flow"|
             wqp$variable == "Aluminum"|wqp$variable == "Bromide"|wqp$variable == "Iron"|wqp$variable == "Nitrate"),]
wqp$value <- as.numeric(wqp$value)
wqp <- reshape2::dcast(wqp, date+siteno+lat+long+NEAR~variable, value.var="value", mean)
wqp <- wqp[complete.cases(wqp[,8:15]),] #4686 complete records

#### Skip to line 258 for compositional PCA ####
#transform data to approx normal dist and convert to mmol from meq
#log transform raw concentration measurements
param <- c("Calcium", "Chloride", "Potassium", "Magnesium", "Sodium", "Sulfate", "Silica")
charge <- c(2,1,1,2,1,2,1)
for(i in 1:length(param)){
  for(j in 1:nrow(wqp)){
    if(wqp[j,param[i]] != 0){
  wqp[j,param[i]] <- wqp[j,param[i]]/charge[i] #convert to mmol/L
  wqp[j,param[i]] <- log10(wqp[j,param[i]]) #log transform all except HCO3
  }
  }
}

#correlation matrix for wqp and wqbc
library(GGally)
ggpairs(wqp, columns= 8:15)
#scale the data - mean of 0 and scale to unit standard deviation
raw.scaled <- as.data.frame(scale(wqp[,8:15], #select all concentration columns
                      center=TRUE,scale=TRUE),stringsAsFactors = FALSE)
base::summary(raw.scaled)

#examine data ranges - create range plots
#Raw
ranges = vapply(raw.scaled, FUN=function(col) c(min(col), max(col)), numeric(2))
rownames(ranges) = c("vmin", "vmax") 
rframe = as.data.frame(t(ranges))  # make ymin/ymax the columns
rframe$varName = rownames(rframe)

setwd(fname)
pdf("CS.raw.varranges.pdf")
barbell_plot(rframe,"varName", "vmin", "vmax") +
  coord_flip() + ggtitle("x scaled variables: Colorado raw concentration")
dev.off()

#calculate MCD
library(robustbase)
mcd <- covMcd(raw.scaled, cor=T)
covPlot(raw.scaled, which="all", labels.id=wqp$NEAR)

raw.rob <- raw.scaled[mcd$best,] #subset to "best" points
#principal components analysis
#Raw
raw.pc <- as.matrix(raw.scaled)
princ <- prcomp(raw.pc, center = TRUE,scale. = TRUE) 
base::summary(princ)
#raw robust
raw.pc1 <- as.matrix(raw.rob)
princ1 <- prcomp(raw.pc1, center = TRUE,scale. = TRUE) 
base::summary(princ1)

#plot magnitude of PCs
#Raw
setwd(fname)
pdf("CS_PCs.pdf")
fviz_eig(princ, title = "Raw PCs")
dev.off()

#Compositional
setwd(fname)
pdf("CSnorm_PCs.pdf")
fviz_eig(princ1, title = "Colorado Surface Water PCs, Robust")
dev.off()


#look at variable loadings for first 5 pcs
#Raw
rot5 <- extractProjection(5,princ)
rotf = as.data.frame(rot5)
rotf$varName = rownames(rotf)
rotflong = gather(rotf, "PC", "loading", starts_with("PC"))

ttest <- lapply(abs(rotf[,1:5]),t.test)
rotflong <- split(rotflong, rotflong$PC)
for(i in 1:length(ttest)){
  for(j in 1:nrow(rotflong[[i]])){
    if(abs(rotflong[[i]]$loading[j]) > ttest[[i]]$conf.int[2]){
      rotflong[[i]]$color[j] <- "significant"
    } else {
      rotflong[[i]]$color[j] <- "not significant"
    }
  }
}
rotflong <- do.call(rbind,rotflong)

setwd(fname)
pdf("CS.varloadings.5pc.pdf")
dotplot_identity(rotflong, "varName", "loading", colorvar="color") + 
  facet_wrap(~PC,nrow=1) + coord_flip() + 
  scale_color_manual(values=c("significant"="#1b9e77", "not significant" = "black"))+
  ggtitle("x scaled variable loadings, first 5 principal components") 
dev.off()

#Robust
rot5 <- extractProjection(5,princ1)
rotf = as.data.frame(rot5)
rotf$varName = rownames(rotf)
rotflong = gather(rotf, "PC", "loading", starts_with("PC"))

#Run t-tests of variable loadings
ttest <- lapply(abs(rotf[,1:5]),t.test)
rotflong <- split(rotflong, rotflong$PC)
for(i in 1:length(ttest)){
  for(j in 1:nrow(rotflong[[i]])){
    if(abs(rotflong[[i]]$loading[j]) > ttest[[i]]$conf.int[2]){
      rotflong[[i]]$color[j] <- "significant"
    } else {
      rotflong[[i]]$color[j] <- "not significant"
    }
  }
}
rotflong <- do.call(rbind,rotflong)

setwd(fname)
pdf("CSnorm.varloadings.5pc.pdf")
dotplot_identity(rotflong, "varName", "loading", colorvar="color") + 
  facet_wrap(~PC,nrow=1) + coord_flip() + 
  scale_color_manual(values=c("significant" = "#1b9e77", "not significant" = "black"))+
  ggtitle("x scaled variable loadings, first 5 principal components") 
dev.off()

#### skip ####
#repeat sampling of PC, determine lower CI
set.seed(1)
mydf <- CS.pc
times <- 1000
ll <- vector(mode = "list", length = times)
for (i in seq_len(times)) {
  tempdf  <- mydf[sample(nrow(mydf), replace = TRUE), ]
  ll[[i]] <- abs(prcomp(tempdf)$rotation) ## NOTE: abs(...)
}

data.frame(apply(simplify2array(ll), 1:2, quantile, probs = 0.025)) #lower limit PC loadings

data.frame(apply(simplify2array(ll), 1:2, quantile, probs = 0.975)) #upper limit PC loadings
#########

#look at var projected with first 2 PC
#CS
library(ggfortify)
wqp$NEAR <- factor(wqp$NEAR, levels=c(2,3,4,5,6,7,8,9,10,12,13,14,15))
wqp$NEAR <- as.numeric(wqp$NEAR)
wqp$month <- format(as.Date(wqp$date, format="%Y-%m-%d"), "%m")
wqp$season <- NA
wqp$season <- ifelse(wqp$month == "12"|wqp$month == "01"|wqp$month == "02","W",
                     ifelse(wqp$month == "03"|wqp$month == "04"|wqp$month == "05","SP",
                            ifelse(wqp$month == "06"|wqp$month == "07"|wqp$month == "08", "SU", "F")))

setwd(fname)
pdf("CS.projecteddat.pdf")
autoplot(princ, data=wqp, colour='NEAR', size=4,
         loadings=TRUE, loadings.colour="black", loadings.label=TRUE,
         loadings.label.size=5, loadings.label.vjust=0, loadings.label.hjust=0, loadings.label.repel=TRUE,
         loadings.label.color="black")+
  #theme(legend.position="none")+
  #ggtitle("Colorado Surface Water, robust")+
  viridis::scale_color_viridis()+
  labs(color="Zone")+
  theme_bw()+
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12))
dev.off()

### Compositional Data analysis ####
#transform to compositional
wqp <- na.omit(reshape2::melt(wqpc, id.vars=c("date", "siteno", "lat", "long", "CBerror", "NEAR")))
wqp <- split(wqp, wqp$variable) #remove Phosphorus, Nitrate, Alkalinity, and Flow
wqp <- do.call(rbind,wqp)
wqp <- wqp[!(wqp$variable == "Phosphorus"|wqp$variable == "Alkalinity"|wqp$variable == "Flow"|
               wqp$variable == "Aluminum"|wqp$variable == "Bromide"|wqp$variable == "Iron"|wqp$variable == "Nitrate"),]
wqp$value <- as.numeric(wqp$value)
wqp <- reshape2::dcast(wqp, date+siteno+lat+long+NEAR~variable, value.var="value", mean)

wqp <- wqp[complete.cases(wqp[,8:15]),] #4863 complete records
param <- c("Bicarbonate","Calcium", "Chloride", "Potassium", "Magnesium", "Sodium", "Sulfate", "Silica")

#convert to compositional
wqp[,param] <- wqp[,param]/rowSums(wqp[,param])
rowSums(wqp[,param])

#look at variation matrix for ilr data
library(robCompositions)
var <- variation(wqp[,param])  #create variation matrix
print(var)
autoplot(var)
setwd(dname)
write.csv(var, "CS_comp_variationmatrix_robpivot.csv")

#look at pivot coordinates as well once understand better
#outlier detection
out <- outCoDa(wqp[,param], quantile = 0.975, method = "robust", alpha = 0.5, coda = TRUE)
wqcomp$Mahaldist <- out$mahalDist
wqcomp$outlier <- out$outlierIndex

plot <- ggplot(data=wqp, aes(x=long, y=Mahaldist, color=outlier))+
  geom_point()+
  theme_bw()
print(plot)

#robust compositional PCA
#wqcomp <- lcg
pccomp <-pcaCoDa(wqp[,param])
x <- base::summary(pccomp)
write.csv(x, "pcsummary.csv")

#create biplot
library(ggplot2)
PCoutput <- data.frame(pccomp[['princompOutputClr']]$scores)
PCload <- data.frame(pccomp[['princompOutputClr']]$loadings)
PCload$x <- 0
PCload$label <- rownames(PCload)
PCoutput$NEAR <- wqp$NEAR
PCoutput$siteno <- wqp$siteno
PCoutput$date <- wqp$date

#### Skip, this is for S1 ####
z <- PCoutput[which(PCoutput$Comp.2 > 2.5),]
z <- z[which(z$siteno == "USGS-08120700"),]
x <- wqp[which(wqp$siteno == "USGS-08121000"),]

ggplot(x, aes(x=date, y=Potassium))+
  geom_point()

wqp$NEAR <- factor(wqp$NEAR, levels=c(2,3,4,5,6,7,8,9,10,12))
wqp$trib <- factor(wqp$trib, levels=c(0,1))
wqp$tlabel <- ifelse(wqp$trib == 1, "Tributary", "Main Stem")
wqp$month <- format(as.Date(wqp$date, format="%Y-%m-%d"), "%m")
wqp$season <- NA
wqp$season <- ifelse(wqp$month == "12"|wqp$month == "01"|wqp$month == "02","W",
                     ifelse(wqp$month == "03"|wqp$month == "04"|wqp$month == "05","SP",
                            ifelse(wqp$month == "06"|wqp$month == "07"|wqp$month == "08", "SU", "F")))

setwd(fname)
jpeg("PCAszn1.jpeg",width=1000, height=800, quality=100)
plot <- ggplot()+
  geom_point(data=PCoutput, aes(x=Comp.1, y=Comp.2, color=wqp$season), size=4.5, alpha=0.5)+
  #scale_shape_manual(name="Position",values=c("Main Stem"=16, "Tributary"=3))+
  geom_segment(data=PCload, aes(x=x, y=x, xend=4*Comp.1, yend=4*Comp.2),
               arrow = arrow(length = unit(0.1,"cm")))+
  geom_text(data=PCload, aes(x=5.5*Comp.1, y=5.5*Comp.2, label=label), size=6)+
  labs(color="Zone",shape="Position", x="PC1 (76%)", y="PC2 (12%)")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=20))
print(plot)
dev.off()
############

#t test significance of loadings 
rotf1 = as.data.frame(pccomp[['princompOutputClr']]$loadings)
rotf1$varName = rownames(rotf1)
rotflong1 = reshape2::melt(rotf1[,c(1:5,8)], id.vars=c("varName"))
colnames(rotflong1) <- c("varName", "PC", "loading")

#Run t-tests of variable loadings
ttest1 <- lapply(abs(rotf1[,1:5]),t.test, conf.level=0.9)
rotflong1 <- split(rotflong1, rotflong1$PC)
for(i in 1:length(ttest1)){
  for(j in 1:nrow(rotflong1[[i]])){
    if(abs(rotflong1[[i]]$loading[j]) > ttest1[[i]]$conf.int[2]){
      rotflong1[[i]]$color[j] <- "significant"
    } else {
      rotflong1[[i]]$color[j] <- "not significant"
    }
  }
}
rotflong1 <- do.call(rbind,rotflong1)
z <- as.data.frame(t(rotflong1))

#save as csv
setwd(dname)
write.csv(rotflong1, 'rPCAloadingsig.csv')

#plot loadings
z <- pccomp[['princompOutputClr']]$loadings
dotplot_identity(rotflong1, "varName", "loading", colorvar="color") + 
  facet_wrap(~PC,nrow=1) + coord_flip() + 
  scale_color_manual(values=c("significant"="#1b9e77", "not significant" = "black"))+
  ggtitle("variable loadings, first 5 principal components") +
  theme_bw()

#biplot
setwd(fname)
jpeg("CR_PCA_06212022.jpeg",width=1000, height=800, quality=100)
plot <- ggplot()+
  geom_point(data=PCoutput, aes(x=Comp.1, y=Comp.2, color=wqp$NEAR), size=4.5, alpha=0.5)+
  geom_segment(data=PCload, aes(x=x, y=x, xend=4*Comp.1, yend=4*Comp.2),
               arrow = arrow(length = unit(0.1,"cm")))+
  geom_text(data=PCload, aes(x=5.5*Comp.1, y=5.5*Comp.2, label=label), size=6)+
  labs(color="Zone",shape="Position", x="PC1 (76%)", y="PC2 (12%)")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=20))
print(plot)
dev.off()

## Don't use for this dataset but save:: Code for box-cox transformation ##
#calculate compositional magnitude as separate df
param <- c("Bicarbonate","Calcium", "Chloride", "Potassium", "Magnesium", "Sodium", "Sulfate", "Silica")
wqcomp <- wqp[,1:5]
lambda <- c()
i <- 1
j <- 1
for(i in 1:length(param)){
  for(j in 1:nrow(wqp)){
    wqcomp[j,param[i]] <- wqp[j,param[i]]/rowSums(wqp[j,6:13])
  }
  x <- wqcomp[!(wqcomp[,param[i]] == 0),] #add this because there are 7 0 measurements for silica which cause an error
  x <- EnvStats::boxcox(x[,param[i]], optimize=TRUE)
  lambda[i] <- x$lambda
}
rowSums(wqcomp[,6:13]) #all sum to 1
#apply box cox transformation using lambdas
wqbc <- wqp[,1:5]
lambda <- round(lambda,1)
for(i in 1:length(param)){
  if(lambda[i] == 0.0){
    wqbc[,param[i]] <- log10(wqcomp[,param[i]])
  }else{
    wqbc[,param[i]] <- ((wqcomp[,param[i]]^lambda[i])-1)/lambda[i]
  }
}

