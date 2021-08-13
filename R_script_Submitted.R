setwd("D:\\")
###Create a Figure folder in your directory

###Add Data_S1.xlsx" into your directory 
###Export the first sheet with local data as TableS1.csv 
###Export the second sheet with migrant data as TableS2.csv



#!! Function to check whether a package is installed
is.installed <- function(mypkg){ is.element(mypkg, installed.packages()[,1]) } 

#!! install packages needed
########################## INSTALL PACKAGES ##################################################

if (!is.installed("raster"))              install.packages("raster",dependencies = TRUE)
if (!is.installed("rgdal"))               install.packages("rgdal",dependencies = TRUE)
if (!is.installed("maps"))                install.packages("maps",dependencies = TRUE)
if (!is.installed("rnaturalearth"))       install.packages("rnaturalearth",dependencies = TRUE)
if (!is.installed("rnaturalearthdata"))   install.packages("rnaturalearthdata",dependencies = TRUE)
if (!is.installed("assignR"))             install.packages("assignR",dependencies = TRUE)
if (!is.installed("cowplot"))       install.packages("cowplot",dependencies = TRUE)
if (!is.installed("rasterVis"))     install.packages("rasterVis",dependencies = TRUE)
if (!is.installed("colorRamps"))    install.packages("colorRamps",dependencies = TRUE)
if (!is.installed("fossil"))        install.packages("fossil",dependencies = TRUE)
if (!is.installed("mapmisc"))       install.packages("mapmisc",dependencies = TRUE)
if (!is.installed("gridExtra"))     install.packages("gridExtra",dependencies = TRUE)
if (!is.installed("ggspatial"))     install.packages("ggspatial",dependencies = TRUE)
if (!is.installed("ggplot2"))       install.packages("ggplot2",dependencies = TRUE)
if (!is.installed("mixtools"))      install.packages("mixtools",dependencies = TRUE)
if (!is.installed("maptools"))      install.packages("maptools",dependencies = TRUE)
if (!is.installed("mvnmle"))        install.packages("nvnmle",dependencies = TRUE)
if (!is.installed("sf"))            install.packages("sf",dependencies = TRUE)
if (!is.installed("sp"))            install.packages("sp",dependencies = TRUE)
if (!is.installed("readxl"))            install.packages("readxl",dependencies = TRUE)

##########################SET LIBRARIES##################################################
library("rgdal")#
library("raster")#
library("maps")#
library("maptools")#
library("rasterVis")#
library("mvnmle")#
library("mixtools")
library("fossil")#
library(cowplot)#
library("assignR")#
library(ggplot2)#
library(colorRamps)#
library(mapmisc)#
library("ggspatial")#
library("rnaturalearth")#
library("rnaturalearthdata")#
library(sf)
library(sp)
library(readxl)


#!! Choose what is relevant for your OS
os="unix"  # this includes osx on a Mac
os="win"

#!! Helper function to store (print) files in sub folders in different OS's
#!! so far only implemented in plots at the beginning of the script but could be applied alter too if so desired.
os_path <- function(folder_path,filename){
  if (os=="win") return(paste(folder_path,"\\",filename,sep="")) else
    return(paste(folder_path,"/",filename,sep="")) } 

#!! Create or use sub folders to keep things tidy!
Fig_Path <- "Figure"  
Out_Path <- "Output"
Ras_Path <- "Projected_rasters"

mydir  = getwd() 
dir.create(file.path(mydir,Fig_Path, fsep = .Platform$file.sep), showWarnings = FALSE)
dir.create(file.path(mydir,Ras_Path, fsep = .Platform$file.sep), showWarnings = FALSE)
dir.create(file.path(mydir,Out_Path, fsep = .Platform$file.sep), showWarnings = FALSE)


####################################################################################
## useful raster functions
####################################################################################
# subset(stack,which()) -- select layers from a stack
# density(raster) -- plots density of cell values
# quantile(rster,probs=) -- finds values at specificed quantiles
# stackApply(stack,fun) -- compute across cells for each layer in a stack
# calc(stack,fun) -- compute across layers for each cell
# cellStats(raster,fun) -- compute across cells for single raster

####################################################################################
# Define new functions
####################################################################################

## function to rescale precip for butterfly wing, generic rescaling from Brattström et al. 2010
#rescale <- function(x){1.1*x-41}
## generic function to rescale using lm output
#rescale <- function(x,b0,b1){b1*x+b0}

## function returns raster of posterior probabilities for univariate normal data
calcCellProb <- function(x,isoscape,std){
  m <- getValues(isoscape)
  s <- getValues(std) # use this is you have raster of variances
  m <- dnorm(x,mean=m,sd=s)
  cell.dens <- setValues(isoscape,m)
  return(cell.dens)
}

## function returns raster of posterior probabilities for bivariate normal data
## x is the unknown tissue of interest, will have two values, one for each isotope
## m is a 2-D vector, all the values in the raster for each isotope
## v is the same as m, but for variances
## r is a single number - the covariance. Can be vector if estimated as non-stationary
## ras is a raster that will serve as a template for the final product
calcCellProb2D <- function(x,m,v,r,ras) {
    pd <- 1/(2*pi*sqrt(v[,1])*sqrt(v[,2])*sqrt(1-r^2))*exp(-(1/(2*(1-r^2)))*
    ((x[1]-m[,1])^2/v[,1]+(x[2]-m[,2])^2/v[,2]-(2*r*(x[1]-m[,1])*
    (x[2]-m[,2]))/(sqrt(v[,1])*sqrt(v[,2]))))
    pdras <- setValues(ras,pd)
    return(pdras)
}

## function returns raster of posterior probability distribution
calcPostProb <- function(x){
    pp <- x/cellStats(x,sum)
    return(pp)
}

## function returns raster identifying upper quantiles across raster values
qtlRaster <- function(ras,q=0.95) {
    qras <- ras >= quantile(ras,probs=q)
    return(qras)
}

## function returns raster of log likelihoods for each cell
llRaster <- function(x,m,v,r,ras) {
    n <- nrow(x)
    k <- ncol(x)
    for (i in 1:nrow(m)) {
        mu <- m[i,]
        cv <- matrix(v[i,1],r,v[i,2],r,nrow=2)
        for (s in 1:n) {
            xs <- x[s,]
            kernSum  <- crossprod(x,solve(cv))%*%x
        }
        ll[i] <- - (n * k / 2) * (log(2*pi)) - (n/2) * log(abs(cv)) - kernSum/2
    }
    llras <- setValues(ras,ll)
    return(llras)
}

## function to get normalized cell probabilites
calcNormProb <- function(x){
    np <- x/cellStats(x,max)
    return(np)
}

## function to convert degrees (DD) to radians
deg2rad <- function(deg) return(deg*pi/180)

## function to get geodesic distance between two points using Spherical Law of Cosines
## xy specified by radian latitude/longitude 
calcDist <- function(long1, lat1, long2, lat2) {
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  R <- 6371 # earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # distance in km
}

####################################################################################
## Data input: read in tissue isotope data, GIS raster models
####################################################################################

###Input local hawck dataset

iso.data<-read.csv("TableS1.csv")

data(wrld_simpl)



### Calibrate a d2H isoscape for birds using the assignR package
d2H_local<-subset(iso.data,select=c("Lat","Long","d2Hrev"))
d2H_local$sd<-d2H_local$d2Hrev/d2H_local$d2Hrev*2
d2H_local$Site_ID<-iso.data$Site_ID
rename(d2H_local, c("Lat", "Long", "marker", "marker.sd", "Site_ID"))
d2H_local<-subset(d2H_local,d2H_local$Site_ID != "L18")  ###Remove one outlier
d2H_local<-subset(d2H_local,d2H_local$Site_ID != "L11") ###Remove one outlier
coordinates(d2H_local) <- ~Long+Lat
crs(d2H_local)<-"+proj=longlat +datum=WGS84"
par(mfrow=c(1,2))
r = calRaster(known = d2H_local, isoscape = d2h_world)
par(mfrow=c(1,3))
plot(r$isoscape.rescale)

x<-assignR::QA(d2H_local, d2h_world, bySite=FALSE, valiStation = 5, valiTime = 50, by=2, mask=NULL,
               setSeed = TRUE, name = NULL)


pdf(os_path(Fig_Path,"Fig2A.pdf"),width=12, height=4)
par(mar=c(2,2,2,2))
par(mfrow=c(1,4))
plot(x)
dev.off()



###Import strontium isoscape from Bataille et al. 2020 
### The global Sr isoscape is available at: https://drive.google.com/open?id=1g9rCGo3Kd3hz2o5JKkSbgNsGJclvsuQm

rf_sr<-raster("rf_plantsoilmammal1.tif")
sr.se<-raster("srse.tif")


# Reproject North America mask to same projection as strontium isoscape
America<-c(-168,-50,20,75)
rangemap<-crop(r$isoscape.rescale$mean/r$isoscape.rescale$mean, America, snap='near')
rf_sr<-raster("rf_plantsoilmammal1.tif")
sr.se<-raster("srse.tif")+0.0005

# Make a SpatialPolygon from the extent of r2
r2extent <- as(extent(rangemap), 'SpatialPolygons')
# Assign this SpatialPolygon the good projection
proj4string(r2extent) <- proj4string(rangemap)
# Transform the projection to that of r1
r2extr1proj <- spTransform(r2extent, CRS(proj4string(rf_sr)))

sr<-crop(rf_sr, r2extr1proj, snap='near')
sr.se<-crop(sr.se,sr,snap='near')

mask<-c(-1.2e7,-4e6,2.6e6,8e6)
rangemap<-sr/sr
rangemap<-crop(rangemap, mask, snap='near')
sr<-crop(sr,rangemap,snap='near')
sr.se<-crop(sr.se,rangemap,snap='near')
na.dH_Hobson2012<-projectRaster(r$isoscape.rescale, rangemap, method="ngb")

iso.data <- read.csv("TableS2.csv")

r_Sr<-brick(sr,sr.se)

Sr_local<-subset(iso.data,select=c("Lat","Long","Sr"))
Sr_local$sd<-iso.data$Sr/iso.data$Sr*0.00005
Sr_local$Site_ID<-iso.data$Site_ID
rename(Sr_local, c("Lat", "Long", "marker", "marker.sd", "Site_ID"))
Sr_local<-subset(Sr_local,Sr_local$Sr != 0.70725)###Remove one outlier Port Mackenzie
Sr_local<-subset(Sr_local,Sr_local$Site_ID != "L15")###Remove one outlier Kellogs
coordinates(Sr_local) <- ~Long+Lat
crs(Sr_local)<-"+proj=longlat +datum=WGS84"
Sr_local_proj<-spTransform(Sr_local, CRS(proj4string(rf_sr)))
par(mfrow=c(1,2))
rSr = calRaster(known = Sr_local_proj, isoscape = r_Sr)
par(mfrow=c(1,1))
plot(rSr$isoscape.rescale)

###test isoscape
x1<-QA(Sr_local_proj, r_Sr, bySite=FALSE, valiStation = 5, valiTime = 5, by=2,
       mask = NULL, setSeed = TRUE, name = NULL)

pdf(os_path(Fig_Path,"Fig2B.pdf"),width=12, height=4)
#par(mar=c(2,2,2,2))
par(mfrow=c(1,4))
plot(x1)
dev.off()


###Read breeding range###
Cooper1 <- read_sf("C:\\Users\\Clement\\Desktop\\Geolocation_hawcks_2019\\Breeding_range\\Accipiter_cooperii.shp")
SS1<-read_sf("C:\\Users\\Clement\\Desktop\\Geolocation_hawcks_2019\\Breeding_range\\Accipiter_striatus.shp")
Merlin1<-read_sf("C:\\Users\\Clement\\Desktop\\Geolocation_hawcks_2019\\Breeding_range\\Falco_columbarius.shp")

#Cooper1 <- calc(Cooper1, fun=function(x){ x[x = 0] <- NA; return(x)} )
#SS1 <- calc(SS1, fun=function(x){ x[x = 0] <- NA; return(x)} )
#Merlin1 <- calc(Merlin1, fun=function(x){ x[x = 0] <- NA; return(x)} )

#Cooper1_p<-projectRaster(Cooper1, rangemap, method="ngb")
#SS1_p<-projectRaster(SS1, rangemap, method="ngb")
#Merlin1_p<-projectRaster(Merlin1, rangemap, method="ngb")

Cooper1_p<-st_transform(Cooper1, crs(rangemap))
Cooper1_p<-st_crop(Cooper1_p, extent(rangemap))
SS1_p<-st_transform(SS1, crs(rangemap))
SS1_p<-st_crop(SS1_p, extent(rangemap))
Merlin1_p<-st_transform(Merlin1, crs(rangemap))
Merlin1_p<-st_crop(Merlin1_p, extent(rangemap))



###Fig 2A and B
wrld_simpl_p<-spTransform(wrld_simpl, crs(sr))
wrld_simpl_p<-gSimplify(wrld_simpl_p, tol = 0.00001)
wrld_simpl_p<-gBuffer(wrld_simpl_p, byid=TRUE, width=0)
rgeos::set_RGEOS_CheckValidity(2L)
wrld_simpl_p<-crop(wrld_simpl_p, extent(sr))
breakpoints<-c(0.703,0.705,0.707,0.709,0.711,0.713,0.715,0.720,0.730,0.750)
breakpoints2<-seq(-200,20,20)
pdf(os_path(Fig_Path,"Fig1.pdf"),width=11, height=5.5)
par(mfrow=c(1,2))
plot(rSr$isoscape.rescale$mean, col= matlab.like2(10), breaks=breakpoints, axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
points(iso_proj,pch=10,col="black",cex=1)
text(-1.4E7, 1E7, "A",cex = 2)
plot(na.dH_Hobson2012$mean, col=matlab.like(11), breaks=breakpoints2, axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
points(iso_proj1,pch=12,col="black",cex=1)
text(-1.4E7, 1E7, "B",cex = 2)
#scaleBar(crs(na.dH_Hobson2012$mean), "topright",cex=1, seg.len=2,box.color = NULL)
dev.off()

####################################################################################
## GEOGRAPHIC ASSIGMENT FOR LOCAL BIRDS
####################################################################################
iso.data <- read.csv("TableS2.csv") # tissue isotope data
iso.data$d2Hrev<-iso.data$d2H*0.8517+10.77
iso.data<-iso.data[order(iso.data$Species),]
iso_proj<- project(as.matrix(iso.data[,c("Long","Lat")]), "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
iso.data$X<-iso_proj[,1]
iso.data$Y<-iso_proj[,2]


####################################################################################
## univariate example
# 1. make assignment model
# 2. fetch probability at known origin
# 3. fetch distance from origin to max of posterior probability
# To assign unknown samples, rem lines related to known origin
####################################################################################

distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$d2Hrev[i],na.dH_Hobson2012$mean,na.dH_Hobson2012$sd) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  # rem the next 5 lines if assigning unknown origin samples
  #max.val <- xyFromCell(pp,which.max(pp)) # location(s) of max post prob
  #d <- vector() # declare vector for distances
  #for(j in 1:nrow(max.val)) # collect distances from origin to max (can be >1)
  #  d[j] <- calcDist(iso.data$Long[i],iso.data$Lat[i],max.val[j,1],max.val[j,2])
  #distance.list[[i]] <- d
}
#rm(i,pp,origins,max.val,j,d)
plot(origins)
breakpoints4<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
freq_h<-hist(origins,breakpoints4)


pdf(os_path(Fig_Path,"maps_d2H_local.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$Site_ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)
  #text(-1E7, 8E6, paste(iso.data$Sex[i]),cex = 1)
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

###########################SR#########################################################
## univariate example
# 1. make assignment model
# 2. fetch probability at known origin
# 3. fetch distance from origin to max of posterior probability
# To assign unknown samples, rem lines related to known origin
####################################################################################
distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$Sr[i],rSr$isoscape.rescale$mean,rSr$isoscape.rescale$sd) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  # rem the next 5 lines if assigning unknown origin samples
  max.val <- xyFromCell(pp,which.max(pp)) # location(s) of max post prob
  #d <- vector() # declare vector for distances
  #for(j in 1:nrow(max.val)) # collect distances from origin to max (can be >1)
  #  d[j] <- calcDist(iso.data$Long[i],iso.data$Lat[i],max.val[j,1],max.val[j,2])
  #distance.list[[i]] <- d
}
#rm(i,pp,origins,max.val,j,d)

plot(origins)
breakpoints4<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
freq_sr<-hist(origins,breakpoints4)


# Show maps of prob densities and distances
# this is hard coded for a bird dataset example
# (i.e. calls to text() are specific to the range considered in these data)
pdf(os_path(Fig_Path,"maps_Sr_local.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$Site_ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)
  #text(-1.2E7, 8E6, paste(iso.data$Sex[i]),cex = 1)
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)


#######################SRb= and d2h DUAL ASSIGNMENT#############################################################
## multivariate example
# 1. make assignment model
# 2. fetch probability at known origin
# To assign unknown samples, rem the lines related to the origin
####################################################################################

## function returns raster of posterior probabilities for bivariate normal data
## x is the unknown tissue of interest, will have two values, one for each isotope
## m is a 2-D vector, all the values in the raster for each isotope
## v is the same as m, but for variances
## r is a single number - the covariance. Can be vector if estimated as non-stationary
## ras is a raster that will serve as a template for the final product
calcCellProb2D <- function(x,m,v,r,ras) {
  pd <- 1/(2*pi*sqrt(v[,1])*sqrt(v[,2])*sqrt(1-r^2))*exp(-(1/(2*(1-r^2)))*
                                                           ((x[1]-m[,1])^2/v[,1]+(x[2]-m[,2])^2/v[,2]-(2*r*(x[1]-m[,1])*                                                                                                    (x[2]-m[,2]))/(sqrt(v[,1])*sqrt(v[,2]))))
  pdras <- setValues(ras,pd)
  return(pdras)
}

# get values from isoscapes for mean of model - convert from rasters to vectors
mu <- cbind(getValues(rSr$isoscape.rescale$mean),getValues(na.dH_Hobson2012$mean))

# compute Sr-O variance-covariance matrix for the two rasters and collect covariance
# only do this is you can't estimate the covariance in the tissues of interest
# use covariance between isotopes from tissue of interest if at all possible
rho <- mlest(mu)$sigmahat[1,2]

# combine sampling and model variances into 2-D vector
# var.c and var.n are (stationary) within-site variances (among individuals)
# ignore/remove if you don't have them
#vars <- cbind(getValues(na.precip.se)^2+var.c,getValues(na.precip.se)^2+var.n)
vars <- cbind(getValues(rSr$isoscape.rescale$sd)^2,getValues(na.dH_Hobson2012$sd)^2)
###Create raster mask
rc<-setValues(na.dH_Hobson2012$mean,NA) 


#####################################################################


# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=iso.data[,1])) {
  tm <- c(iso.data$Sr[i],iso.data$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}


freq_srh<-hist(origins,breakpoints4)

pdf(os_path(Fig_Path,"maps_Srd2H_local.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  #plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=2)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$Site_ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)
}
par(opar)
dev.off()
rm(i)

pdf(os_path(Fig_Path,"local_dual_performance.pdf",width=5, height=6))
par(mfrow=c(1,1))
plot(freq_srh$mids,(1-cumsum(freq_srh$counts)/sum(freq_srh$counts))*100,ty="l",col="darkgreen",lwd=2,log="y",xlab="Probability",ylab="Proportion of area",ylim=c(0.01,100))
lines(freq_sr$mids,(1-cumsum(freq_sr$counts)/sum(freq_sr$counts))*100,ty="l",col="red",lwd=2, log="y",add=TRUE)
lines(freq_h$mids,(1-cumsum(freq_h$counts)/sum(freq_h$counts))*100,ty="l",col="blue",lwd=2, log="y",add=TRUE)
dev.off()


####################################################################################
## MIGRANT BIRDS
####################################################################################
iso.data <- read.csv("TableS3.csv") # tissue isotope data
iso.data$d2Hrev<-iso.data$d2H*0.8517+10.77
iso.data<-iso.data[order(iso.data$Species),]
iso_proj<- project(as.matrix(iso.data[,c("Long","Lat")]), "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
iso.data$X<-iso_proj[,1]
iso.data$Y<-iso_proj[,2]



####################################################################################
## univariate example
# 1. make assignment model
# 2. fetch probability at known origin
# 3. fetch distance from origin to max of posterior probability
# To assign unknown samples, rem lines related to known origin
####################################################################################

distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$d2Hrev[i],na.dH_Hobson2012$mean,na.dH_Hobson2012$sd) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  # rem the next 5 lines if assigning unknown origin samples
  #max.val <- xyFromCell(pp,which.max(pp)) # location(s) of max post prob
  #d <- vector() # declare vector for distances
  #for(j in 1:nrow(max.val)) # collect distances from origin to max (can be >1)
  #  d[j] <- calcDist(iso.data$Long[i],iso.data$Lat[i],max.val[j,1],max.val[j,2])
  #distance.list[[i]] <- d
}
#rm(i,pp,origins,max.val,j,d)

breakpoints4<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
freq_h<-hist(origins,breakpoints4)


####################################################################################
## summary map plots
####################################################################################
distance.list <- vector("list", length(iso.data[,1]))
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
#spplot(origins,col.regions=rev(terrain.colors(20))) # this will make matrix
par(opar)

#beginCluster()
Sharped=origins$mean.1.1.1.1.2+origins$mean.2.1.1.1.2+origins$mean.1.2.1.1.2+origins$mean.2.2.1.1.2+origins$mean.1.1.2.1.2+origins$mean.2.1.2.1.2+origins$mean.1.2.2.1.2+origins$mean.2.2.2.1.2+origins$mean.1.1.1.2.2+origins$mean.2.1.1.2.2+origins$mean.1.2.1.2.2+origins$mean.2.2.1.2.2+origins$mean.1.1.2.2.2+origins$mean.2.1.2.2.2+origins$mean.1.2.2.2.2+origins$mean.2.2.2.2.2+origins$mean.1.1+origins$mean.2.1
Cooper=origins$mean.2.1.2.1.1+origins$mean.1.2.2.1.1+origins$mean.2.2.2.1.1+origins$mean.1.1.1.2.1+origins$mean.2.1.1.2.1+origins$mean.1.2.1.2.1+origins$mean.2.2.1.2.1+origins$mean.1.1.2.2.1+origins$mean.2.1.2.2.1+origins$mean.1.2.2.2.1+origins$mean.2.2.2.2.1
TMerlin=origins$mean+origins$mean.2+origins$mean.1+origins$mean.2.2
Bmerlin=origins$mean.1.1.1.1.1+origins$mean.2.1.1.1.1+origins$mean.1.2.1.1.1+origins$mean.2.2.1.1.1+origins$mean.1.1.2.1.1
Merlin=Bmerlin+TMerlin

Sharped=100*Sharped/18
Cooper=100*Cooper/11
TMerlin=100*TMerlin/4
Bmerlin=100*Bmerlin/5
Merlin=100*Merlin/9

#endCluster()

plot(st_geometry(nc), col = sf.colors(12, categorical = TRUE), border = 'grey', 
     axes = TRUE)

cols <- c("darkblue","darkred","gold2","darkorange")
cols1<-c("darkblue","gold2","darkred","darkorange")
pal.1 <- colorRampPalette(c("white","gray30"), bias=1)
breakpoints3<-c(0,10,20,30,40,50,60,70,80,90,100)

pdf("Figures\\d2H_byspecies_nofish.pdf",width=11, height=3)
par(mfrow=c(1,3))
par(mar=c(1,1,1,1))
plot(Sharped,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Sharped Shinned",cex = 1)
par(mar=c(1,1,1,1))
plot(Cooper,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cooper",cex = 1)
par(mar=c(1,1,1,1))
plot(Merlin,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols1,lty=c(3,6,1,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Merlin",cex = 1)
dev.off()

pdfos_path(Fig_Path,"maps_d2H_migrant.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)
  #text(-1E7, 8E6, paste(iso.data$Sex[i]),cex = 1)
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

pdf(os_path(Fig_Path,"maps_d2H_migrant_SS.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  #plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.5E6,paste(iso.data$ID[i]))
  text(-6E6, 8.5E6, paste(iso.data$Species[i]),cex = 1)
}
par(opar)
dev.off()
rm(i)

pdf(os_path(Fig_Path,"maps_d2H_migrant_Cooper.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  #plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.5E6,paste(iso.data$ID[i]))
  text(-6E6, 8.5E6, paste(iso.data$Species[i]),cex = 1)}
par(opar)
dev.off()

pdf(os_path(Fig_Path,"maps_d2H_migrant_Merlin.pdf"))
opar<-par()
par(mfrow=c(2,2),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)}
par(opar)
dev.off()




###########################SR#########################################################
## univariate example
# 1. make assignment model
# 2. fetch probability at known origin
# 3. fetch distance from origin to max of posterior probability
# To assign unknown samples, rem lines related to known origin
####################################################################################
distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$Sr[i],rSr$isoscape.rescale$mean,rSr$isoscape.rescale$sd) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  # rem the next 5 lines if assigning unknown origin samples
  max.val <- xyFromCell(pp,which.max(pp)) # location(s) of max post prob
  #d <- vector() # declare vector for distances
  #for(j in 1:nrow(max.val)) # collect distances from origin to max (can be >1)
  #  d[j] <- calcDist(iso.data$Long[i],iso.data$Lat[i],max.val[j,1],max.val[j,2])
  #distance.list[[i]] <- d
}
#rm(i,pp,origins,max.val,j,d)

freq_sr<-hist(origins,breakpoints4)

#beginCluster()
Sharped_Sr=origins$mean.1.1.1.1.2+origins$mean.2.1.1.1.2+origins$mean.1.2.1.1.2+origins$mean.2.2.1.1.2+origins$mean.1.1.2.1.2+origins$mean.2.1.2.1.2+origins$mean.1.2.2.1.2+origins$mean.2.2.2.1.2+origins$mean.1.1.1.2.2+origins$mean.2.1.1.2.2+origins$mean.1.2.1.2.2+origins$mean.2.2.1.2.2+origins$mean.1.1.2.2.2+origins$mean.2.1.2.2.2+origins$mean.1.2.2.2.2+origins$mean.2.2.2.2.2+origins$mean.1.1+origins$mean.2.1
Cooper_Sr=origins$mean.1.2.2.1.1+origins$mean.2.2.2.1.1+origins$mean.1.1.1.2.1+origins$mean.2.1.1.2.1+origins$mean.1.2.1.2.1+origins$mean.2.2.1.2.1+origins$mean.1.1.2.2.1+origins$mean.2.1.2.2.1+origins$mean.1.2.2.2.1+origins$mean.2.2.2.2.1
TMerlin_Sr=origins$mean+origins$mean.2+origins$mean.1+origins$mean.2.2
Bmerlin_Sr=origins$mean.1.1.1.1.1+origins$mean.2.1.1.1.1+origins$mean.1.2.1.1.1+origins$mean.2.2.1.1.1+origins$mean.1.1.2.1.1
Merlin_Sr=Bmerlin_Sr+TMerlin_Sr





Sharped_Sr=100*Sharped_Sr/18
Cooper_Sr=100*Cooper_Sr/10
TMerlin_Sr=100*TMerlin_Sr/4
Bmerlin_Sr=100*Bmerlin_Sr/5
Merlin_Sr=100*Merlin_Sr/9



#endCluster()


pdf(os_path(Fig_Path,"Sr_byspecies_nofish.pdf",width=11, height=3))
par(mfrow=c(1,3))
par(mar=c(1,1,1,1))
plot(Sharped_Sr,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Sharped Shinned",cex = 1)
par(mar=c(1,1,1,1))
plot(Cooper_Sr,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cooper",cex = 1)
par(mar=c(1,1,1,1))
plot(Merlin_Sr,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols1,lty=c(3,6,1,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Merlin",cex = 1)
dev.off()
####################################################################################
## summary map plots
####################################################################################
distance.list <- vector("list", length(iso.data[,1]))
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
#spplot(origins,col.regions=rev(terrain.colors(20))) # this will make matrix
par(opar)

#spplot(origins[[1]],col.regions=rev(terrain.colors(20))) # just one individual

# Show maps of prob densities and distances
# this is hard coded for a bird dataset example
# (i.e. calls to text() are specific to the range considered in these data)
pdf(os_path(Fig_Path,"maps_Sr_migrant.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)
  #text(-1.2E7, 8E6, paste(iso.data$Sex[i]),cex = 1)
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

pdf(os_path(Fig_Path,"maps_Sr_migrant_SS.pdf"))
opar<-par()
par(mfrow=c(2,2),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.5E6,paste(iso.data$ID[i]))
  text(-6E6, 8.5E6, paste(iso.data$Species[i]),cex = 1)
}
par(opar)
dev.off()
rm(i)

pdf(os_path(Fig_Path,"maps_Sr_migrant_Cooper.pdf"))
opar<-par()
par(mfrow=c(2,2),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.5E6,paste(iso.data$ID[i]))
  text(-6E6, 8.5E6, paste(iso.data$Species[i]),cex = 1)}
par(opar)
dev.off()

pdf(os_path(Fig_Path,"maps_Sr_migrant_Merlin.pdf"))
opar<-par()
par(mfrow=c(2,2),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.5E6,paste(iso.data$ID[i]))
  text(-6E6, 8.5E6, paste(iso.data$Species[i]),cex = 1)}
par(opar)
dev.off()

#######################SRb= and d2h DUAL ASSIGNMENT#############################################################
## multivariate example
# 1. make assignment model
# 2. fetch probability at known origin
# To assign unknown samples, rem the lines related to the origin
####################################################################################

## function returns raster of posterior probabilities for bivariate normal data
## x is the unknown tissue of interest, will have two values, one for each isotope
## m is a 2-D vector, all the values in the raster for each isotope
## v is the same as m, but for variances
## r is a single number - the covariance. Can be vector if estimated as non-stationary
## ras is a raster that will serve as a template for the final product
calcCellProb2D <- function(x,m,v,r,ras) {
  pd <- 1/(2*pi*sqrt(v[,1])*sqrt(v[,2])*sqrt(1-r^2))*exp(-(1/(2*(1-r^2)))*
                                                           ((x[1]-m[,1])^2/v[,1]+(x[2]-m[,2])^2/v[,2]-(2*r*(x[1]-m[,1])*                                                                                                    (x[2]-m[,2]))/(sqrt(v[,1])*sqrt(v[,2]))))
  pdras <- setValues(ras,pd)
  return(pdras)
}

# get values from isoscapes for mean of model - convert from rasters to vectors
mu <- cbind(getValues(rSr$isoscape.rescale$mean),getValues(na.dH_Hobson2012$mean))

# compute Sr-O variance-covariance matrix for the two rasters and collect covariance
# only do this is you can't estimate the covariance in the tissues of interest
# use covariance between isotopes from tissue of interest if at all possible
rho <- mlest(mu)$sigmahat[1,2]

# combine sampling and model variances into 2-D vector
# var.c and var.n are (stationary) within-site variances (among individuals)
# ignore/remove if you don't have them
#vars <- cbind(getValues(na.precip.se)^2+var.c,getValues(na.precip.se)^2+var.n)
vars <- cbind(getValues(rSr$isoscape.rescale$sd)^2,getValues(na.dH_Hobson2012$sd)^2)
###Create raster mask
rc<-setValues(na.dH_Hobson2012$mean,NA) 


#####################################################################




# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=iso.data[,1])) {
  tm <- c(iso.data$Sr[i],iso.data$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}


Sharped_Srdh=origins$mean.1.1.1.1.2+origins$mean.2.1.1.1.2+origins$mean.1.2.1.1.2+origins$mean.2.2.1.1.2+origins$mean.1.1.2.1.2+origins$mean.2.1.2.1.2+origins$mean.1.2.2.1.2+origins$mean.2.2.2.1.2+origins$mean.1.1.1.2.2+origins$mean.2.1.1.2.2+origins$mean.1.2.1.2.2+origins$mean.2.2.1.2.2+origins$mean.1.1.2.2.2+origins$mean.2.1.2.2.2+origins$mean.1.2.2.2.2+origins$mean.2.2.2.2.2+origins$mean.1.1+origins$mean.2.1
Cooper_Srdh=origins$mean.1.2.2.1.1+origins$mean.2.2.2.1.1+origins$mean.1.1.1.2.1+origins$mean.2.1.1.2.1+origins$mean.1.2.1.2.1+origins$mean.2.2.1.2.1+origins$mean.1.1.2.2.1+origins$mean.2.1.2.2.1+origins$mean.1.2.2.2.1+origins$mean.2.2.2.2.1
TMerlin_Srdh=origins$mean+origins$mean.2+origins$mean.1+origins$mean.2.2
Bmerlin_Srdh=origins$mean.1.1.1.1.1+origins$mean.2.1.1.1.1+origins$mean.1.2.1.1.1+origins$mean.2.2.1.1.1+origins$mean.1.1.2.1.1
Merlin_Srdh=Bmerlin_Srdh+TMerlin_Srdh


Sharped_Srdh=100*Sharped_Srdh/18
Cooper_Srdh=100*Cooper_Srdh/10
TMerlin_Srdh=100*TMerlin_Srdh/4
Bmerlin_Srdh=100*Bmerlin_Srdh/5
Merlin_Srdh=100*Merlin_Srdh/9

breakpoints3<-c(0,10,20,30,40,50,60,70,80,90,100)
#pal.1 <- colorRampPalette(c("white","gold","orange","orangered","red","red3"), bias=1)
pdf(os_path(Fig_Path,"Srdh_byspecies.pdf",width=11, height=3))
par(mfrow=c(1,3))
par(mar=c(1,1,1,1))
plot(Sharped_Srdh,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=0.5)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Sharped Shinned",cex = 1)
par(mar=c(1,1,1,1))
plot(Cooper_Srdh,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cooper",cex = 1)
par(mar=c(1,1,1,1))
plot(Merlin_Srdh,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols1,lty=c(3,6,1,4),lwd=0.5)
points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Merlin",cex = 1)
dev.off()




pdf(os_path(Fig_Path,"maps_Srd2H_migrant.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  #plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=2)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)
}
par(opar)
dev.off()
rm(i)

pdf(os_path(Fig_Path,"maps_Srd2H_migrant_Cooper.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  plot(wrld_simpl_p,add=T)
  #plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)}
par(opar)
dev.off()

pdf(os_path(Fig_Path,"maps_Srd2H_migrant_Merlin.pdf"))
opar<-par()
par(mfrow=c(5,4),mar=c(1,1,1,1))
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]], axes = FALSE)
  plot(wrld_simpl_p,add=T)
  #plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
  # plot known origin and individual label info
  points(iso.data$X[i],iso.data$Y[i],col="black",cex=0.7, pch=1)
  text(-1E7, 8.3E6,paste(iso.data$ID[i]))
  text(-7E6, 2.2E6, paste(iso.data$Species[i]),cex = 1)}
par(opar)
dev.off()

qtlRaster(origins, threshold=0.1, thresholdType = "prob", genplot = TRUE, outDir = NULL)


freq_srh<-hist(origins,breakpoints4)
par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(freq_srh$mids,freq_srh$counts/sum(freq_srh$counts),ty="l",col="darkgreen",log = "y")
lines(freq_sr$mids,freq_sr$counts/sum(freq_sr$counts),ty="l",col="red",log = "y", add=TRUE)
lines(freq_h$mids,freq_h$counts/sum(freq_h$counts),ty="l",col="blue",log = "y", add=TRUE)




pdf(os_path(Fig_Path,"migrant_dual_performance.pdf",width=5, height=6))
par(mfrow=c(1,1))
plot(freq_srh$mids,(1-cumsum(freq_srh$counts)/sum(freq_srh$counts))*100,ty="l",col="darkgreen",lwd=2,log = "y",xlab="Probability",ylab="Proportion of area",ylim=c(0.01,100))
lines(freq_sr$mids,(1-cumsum(freq_sr$counts)/sum(freq_sr$counts))*100,ty="l",col="red",lwd=2,log = "y", add=TRUE)
lines(freq_h$mids,(1-cumsum(freq_h$counts)/sum(freq_h$counts))*100,ty="l",col="blue",lwd=2,log = "y", add=TRUE)
dev.off()




####CLuster
library(NbClust)
cluster<-subset(iso.data, select=c(d2Hrev))
cluster<-cluster[complete.cases(cluster), ]
#how many clusters should we make?
NbClust(data = cluster, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 6, method = "kmeans") #According to the majority rule, the best number of clusters is  3 
#K-Means Cluster Analysis
fit <- kmeans(cluster, 3)

mean_cluster<-aggregate(cluster,by=list(fit$cluster),FUN=mean)#get cluster means 
iso.data <- data.frame(iso.data, fit$cluster)# append cluster assignment
iso.data1<-subset(iso.data, select=c(d2Hrev,Sr,ID,Species,fit.cluster))

Merlin<-iso.data[which (iso.data$Species== "Taiga Merlin"|iso.data$Species== "Black Merlin"),]
Merlin3<-Merlin[which (Merlin$fit.cluster==3),]
Merlin1<-Merlin[which (Merlin$fit.cluster==1),]

SS<-iso.data[which (iso.data$Species== "Sharp-shinned"),]
SS3<-SS[which (SS$fit.cluster==3),]
SS2<-SS[which (SS$fit.cluster==2),]

Cooper<-iso.data[which (iso.data$Species== "Cooper's"),]
Cooper3<-Cooper[which (Cooper$fit.cluster==3),]
Cooper2<-Cooper[which (Cooper$fit.cluster==2),]

# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=Merlin3[,1])) {
  tm <- c(Merlin3$Sr[i],Merlin3$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}


merlin_cluster3<- origins$mean.1.1+origins$mean.2.1+origins$mean.1.2+origins$mean.2.2+origins$mean.1+origins$mean.2+origins$mean
merlin_cluster3<-merlin_cluster3/7


# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=Merlin1[,1])) {
  tm <- c(Merlin1$Sr[i],Merlin1$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}

merlin_cluster1<- origins$mean.1+origins$mean
merlin_cluster1<-merlin_cluster1/2



breakpoints3<-c(0,10,20,30,40,50,60,70,80,90,100)
#pal.1 <- colorRampPalette(c("white","gold","orange","orangered","red","red3"), bias=1)
pdf(os_path(Fig_Path,"Srdh_Merlincluster.pdf",width=7.5, height=3))
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
plot(merlin_cluster1*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols1,lty=c(3,6,1,4),lwd=0.5)
points(Merlin1$X,Merlin1$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Merlin",cex = 1)
plot(merlin_cluster3*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols1,lty=c(3,6,1,4),lwd=0.5)
points(Merlin3$X,Merlin3$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Merlin",cex = 1)
dev.off()


# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=SS2[,1])) {
  tm <- c(SS2$Sr[i],SS2$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}


SS_cluster2<-origins$mean.1.1.1+origins$mean.2.1.1+origins$mean.1.2.1+origins$mean.2.2.1+origins$mean.1.1.2+origins$mean.2.1.2+origins$mean.1.2.2+origins$mean.2.2.2+origins$mean.1+origins$mean.2+origins$mean
SS_cluster2<-SS_cluster2/11


# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=SS3[,1])) {
  tm <- c(SS3$Sr[i],SS3$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}

SS_cluster3<- origins$mean.1.1+origins$mean.2.1+origins$mean.1.2+origins$mean.2.2+origins$mean.1+origins$mean.2+origins$mean
SS_cluster3<-SS_cluster3/7



breakpoints3<-c(0,10,20,30,40,50,60,70,80,90,100)
#pal.1 <- colorRampPalette(c("white","gold","orange","orangered","red","red3"), bias=1)
pdf(os_path(Fig_Path,"Srdh_SS.pdf",width=7.5, height=3))
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
plot(SS_cluster2*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,6,1,4),lwd=0.5)
points(SS2$X,SS2$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "SS",cex = 1)
plot(SS_cluster3*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,6,1,4),lwd=0.5)
points(SS3$X,SS3$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "SS",cex = 1)
dev.off()

# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=Cooper2[,1])) {
  tm <- c(Cooper2$Sr[i],Cooper2$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}


Cooper_cluster2<-origins$mean



# make the assignment models
origins <- stack()
origin.val <- NULL  
for (i in seq(along=Cooper3[,1])) {
  tm <- c(Cooper3$Sr[i],Cooper3$d2Hrev[i])
  rasta <- calcCellProb2D(tm,mu,vars,rho,rc)
  np <- rasta/cellStats(rasta,max)
  origins <- stack(origins, np)
  # rem next 2 lines for unknown origin samples - validation measure for known origins
  #xy <- data.frame(x=iso.data$lon[i],y=iso.data$lat[i])
  #origin.val[i] <- extract(np,xy) # find value at the origin
}

Cooper_cluster3<- origins$mean.2.1.1+origins$mean.1.2.1+origins$mean.2.2.1+origins$mean.1.1.2+origins$mean.2.1.2+origins$mean.1.2.2+origins$mean.2.2.2+origins$mean.1+origins$mean.2
Cooper_cluster3<-Cooper_cluster3/9



breakpoints3<-c(0,10,20,30,40,50,60,70,80,90,100)
#pal.1 <- colorRampPalette(c("white","gold","orange","orangered","red","red3"), bias=1)
pdf(os_path(Fig_Path,"Srdh_Cooper.pdf",width=7.5, height=3))
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
plot(Cooper_cluster2*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,6,1,4),lwd=0.5)
points(Cooper2$X,Cooper2$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cooper",cex = 1)
plot(Cooper_cluster3*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,6,1,4),lwd=0.5)
points(Cooper3$X,Cooper3$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cooper",cex = 1)
dev.off()

cluster1<-merlin_cluster1/2
cluster2<-(11*SS_cluster2+Cooper_cluster2)/12
cluster3<-(9*Cooper_cluster3+7*SS_cluster3+7*merlin_cluster3)/23

breakpoints3<-c(0,10,20,30,40,50,60,70,80,90,100)
#pal.1 <- colorRampPalette(c("white","gold","orange","orangered","red","red3"), bias=1)
pdf(os_path(Fig_Path,"Srdh_bycluster.pdf",width=11, height=3))
par(mfrow=c(1,3))
par(mar=c(1,1,1,1))
plot(cluster1*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
#plot(SS1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=0.5)
#points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cluster 1",cex = 1)
par(mar=c(1,1,1,1))
plot(cluster2*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
#plot(Cooper1_p["SEASONAL"],add=TRUE,col = "NA", border=cols,lty=c(3,1,6,4),lwd=1)
#points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cluster 2",cex = 1)
par(mar=c(1,1,1,1))
plot(cluster3*100,col= pal.1(10), breaks=breakpoints3,axes = FALSE)
plot(wrld_simpl_p,add=TRUE,ylim=c(2599837,7999837))
#plot(Merlin1_p["SEASONAL"],add=TRUE,col = "NA", border=cols1,lty=c(3,6,1,4),lwd=0.5)
#points(iso.data$X,iso.data$Y,col="red",cex=1, pch=1)
text(-5E6, 3E6, "Cluster 3",cex = 1)
dev.off()
