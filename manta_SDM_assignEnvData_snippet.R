##### MANTA SPECIES DISTRIBUTION MODEL #####
#     Developed by Nick Farmer, 
#     Tim Gowan, John Carlson, 
#     Lance Garrison
#
#     Last modified: 4/26/2019
############################################

#### SECTION 1.0 - Packages ####
#1.1 Packages
install.packages("devtools")
install.packages("foreign")
#devtools::install_github("tidyverse/ggplot2")
install.packages("ggplot2")
install.packages("rgdal")
install.packages("data.table")
install.packages("plyr")
install.packages("RColorBrewer")
install.packages("grid")
install.packages("gridExtra")
install.packages("maps")
install_github("olafmersmann/truncnorm")
install.packages("truncnorm")
devtools::install_github("DistanceDevelopment/mrds")
devtools::install_github("DistanceDevelopment/Distance")
install.packages("RColorBrewer")
install.packages("raster")
install.packages("AICcmodavg")
install.packages("rgeos")
install.packages("maptools")
install.packages("mapproj")
install.packages("sp")
install.packages("ggmap")
install.packages("sf")
install.packages("leaflet")
install.packages("ncdf4")
install.packages("maps")
install.packages("matlab")
#install.packages("grDevices")
#devtools::install_github("ropensci/rerddap")
install.packages("lubridate")
install.packages("tibble")
install.packages("marmap")
install.packages("car")
install.packages("pscl")
install.packages("mgcv")
install.packages("glmmADMB",
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
install.packages("imager")
install.packages("gamclass")
install.packages("rlang")
install.packages("dsm")
install.packages("captioner")
install.packages("tweedie")
install.packages("clickR")
install.packages("swfscMisc")
install.packages("remotes")
remotes::install_github("juoe/sdmflow")

#### SECTION 1.2 - Libraries ####
library(sjlabelled)
#library(devtools)
library(foreign)
library(ggplot2)
library(rgdal)        # for readOGR(...)
library(data.table)
library(plyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(maps)
library(truncnorm)
library(mrds)
library(Distance)
library(AICcmodavg)
library(rgeos)
library(maptools)
library(raster)
library(dplyr)
library(mapproj)
library(sp)
library(ggmap)
library(sf)
library(rgdal)
library(ncdf4)
library(matlab)
library(grDevices)
library(rerddap)
library(lubridate)
library(tibble)
library(arm)
library(lme4)
library(MASS)
library(Hmisc) 
library(plyr)  
# library(R2admb)
# library(glmmADMB)
library(usdm)
library(splines)
library(coda)
# library(scapeMCMC)
# library(coefplot2)
library(ncf)
library(ade4)
library(ecodist)
library(boot)
library(reshape)
library(marmap)
library(car)
library(pscl)
library(mgcv)
library(colorRamps)
library(adehabitatHR)
library(imager)
# library(spatial.tools)
library(classInt)
library(httr)
library(grec)
library(tidyselect)
library(sjPlot)
library(gamclass)
library(magick)
library(stringi)
library(mapdata)
library(rasterVis)
library(animation)
# library(mapmate)
library(purrr)
# library(wvtool)
library(oceanmap)
# library(geostatsp)
library(spatstat)
library(rlang)
library(dsm)
library(ROCR)
library(gam)
library(captioner)
library(tweedie)
library(RODBC)
library(knitr)
library(ggpubr)
library(tidyr)
library(grid)
library(gridExtra)
library(swfscMisc)
library(sf)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(tigris)
library(maps)
library(ggsn)
library(ggspatial)

#### SECTION 2 - Initialize workspace ####
set.seed(11123)

# define projections
wgs.84    <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# Albers Equal Area Conic http://spatialreference.org/ref/esri/north-america-albers-equal-area-conic/
albers <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#### SECTION 3 - GOMMAPPS SURVEYS ####

# LPG: YOU CAN IGNORE THIS SECTION, AS I COMBINE ALL SURVEYS BELOW
# ONLY THING TO NOTE HERE IS THAT THE OUTPUTS FROM GOMMAPPS WHEN FIT TO THE ENTIRE MODEL DOMAIN HAD HIGHER AUC THAN OUTPUTS FROM SEFSC SURVEYS COMBINED WHEN FIT TO THE ENTIRE MODEL DOMAIN
# SIMILARLY, THE NARWC DATA, WHICH I CONSIDERED LESS RESOLVED, ALSO HAD HIGHER AUC THAN OUTPUTS FROM THE SEFSC SURVEYS COMBINED
# AND THE NARWC DATA ALSO HAD HIGHER RESIDUAL DEVIANCE EXPLAINED

#### 3.1 GOMMAPPS - Data Imports ####
setwd("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort")    
gommapps_2018f<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="GoMMAPPS_Aerial_2018F_OnEffort_gridded_env",stringsAsFactors = F) # this is in equal area projection
gommapps_2018w<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="GoMMAPPS_Aerial_2018W_OnEffort_gridded_env",stringsAsFactors = F) # this is in equal area projection
gommapps_2017su<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="GoMMAPPS_Aerial_2017Su_OnEffort_grid_env",stringsAsFactors = F) # this is in equal area projection
nrda_2011f<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="NRDA_11FAL_EFFORT_gridded_env",stringsAsFactors = F) # this is in equal area projection
nrda_2011sp<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="NRDA_11SPR_EFFORT_gridded_env",stringsAsFactors = F) # this is in equal area projection
nrda_2011su<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="NRDA_11SUM_EFFORT_gridded_env",stringsAsFactors = F) # this is in equal area projection
nrda_2012w<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="NRDA_12WIN_EFFORT_gridded_env",stringsAsFactors = F) # this is in equal area projection

gommapps_2018f_manta_T1<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Fall 2018/GoMMAPPS_Aerial_2018F_MantaSightings_T1.csv",stringsAsFactors = F)
gommapps_2018f_manta_T1$Team<-1
gommapps_2018f_manta_T2<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Fall 2018/GoMMAPPS_Aerial_2018F_MantaSightings_T2.csv",stringsAsFactors = F)
gommapps_2018f_manta_T2$Team<-2
gommapps_2018w_manta_T1<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Winter 2018/GoMMAPPS_Aerial_2018W_MantaSightings_T1.csv",stringsAsFactors = F)
gommapps_2018w_manta_T1$Team<-1
gommapps_2018w_manta_T2<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Winter 2018/GoMMAPPS_Aerial_2018W_MantaSightings_T2.csv",stringsAsFactors = F)
gommapps_2018w_manta_T2$Team<-2
gommapps_2017s_manta_T1<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Summer 2017/GoMMAPPS_Aerial_2017Su_MantaSightings_T1.csv",stringsAsFactors = F)
gommapps_2017s_manta_T1$Team<-1
gommapps_2017s_manta_T2<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Summer 2017/GoMMAPPS_Aerial_2017Su_MantaSightings_T2.csv",stringsAsFactors = F)
gommapps_2017s_manta_T2$Team<-2

gommapps_2017s_surveyTrack<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Summer 2017/SurveyTrack_T1_FOR_2017s.csv",stringsAsFactors = F)
gommapps_2018f_surveyTrack<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Fall 2018/SurveyTrack_T1_FOR_2018F.csv",stringsAsFactors = F)
gommapps_2018w_surveyTrack<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/Winter 2018/SurveyTrack_T1_FOR_2018w.csv",stringsAsFactors = F)
surveytrack<-rbind(gommapps_2017s_surveyTrack,gommapps_2018f_surveyTrack,gommapps_2018w_surveyTrack)
surveytrack$datetime<-as.POSIXlt(paste(as.Date(surveytrack$Date,format="%Y-%m-%d"),as.character(surveytrack$Time,format="%H:%M:%S")))

gommapps<-rbind(gommapps_2017su,gommapps_2018f,gommapps_2018w)
summary(gommapps)
names(nrda_2011f); names(nrda_2011sp); names(nrda_2011su); names(nrda_2012w)
nrda<-rbind(nrda_2011f,nrda_2011sp,nrda_2011su,nrda_2012w)

gommapps$link<-paste(as.Date(gommapps$Date, "%Y-%m-%d"),"_",gommapps$TRANSECT,"_",gommapps$FID_ModelD, sep="")
gommapps$segmentLength_m<-gLength(gommapps,byid=T)

domain<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray",layer="ModelDomain3_10km_grid")
crs(domain)
crs(gommapps_2018f)

#create link by date, transect, and model grid
gommapps_2018f_manta_T1_spt <- SpatialPointsDataFrame(coords = gommapps_2018f_manta_T1[,c("EntryLon","EntryLat")], data = gommapps_2018f_manta_T1,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

gommapps_2018f_manta_T1_spt<-spTransform(gommapps_2018f_manta_T1_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(gommapps_2018f_manta_T1_spt,domain)
gommapps_2018f_manta_T1_spt@data$ModelD<-test$DomainFID

gommapps_2018f_manta_T1_spt@data$link<-paste(as.Date(gommapps_2018f_manta_T1_spt@data$Date, "%m/%d/%Y"),"_",gommapps_2018f_manta_T1_spt@data$Transect,"_",gommapps_2018f_manta_T1_spt@data$ModelD,sep="")

head(gommapps_2018f_manta_T1_spt@data)

gommapps_2018w_manta_T1_spt <- SpatialPointsDataFrame(coords = gommapps_2018w_manta_T1[,c("EntryLon","EntryLat")], data = gommapps_2018w_manta_T1,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

gommapps_2018w_manta_T1_spt<-spTransform(gommapps_2018w_manta_T1_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(gommapps_2018w_manta_T1_spt,domain)
gommapps_2018w_manta_T1_spt@data$ModelD<-test$DomainFID

gommapps_2018w_manta_T1_spt@data$link<-paste(as.Date(gommapps_2018w_manta_T1_spt@data$Date, "%m/%d/%Y"),"_",gommapps_2018w_manta_T1_spt@data$Transect,"_",gommapps_2018w_manta_T1_spt@data$ModelD,sep="")

gommapps_2017s_manta_T1_spt <- SpatialPointsDataFrame(coords = gommapps_2017s_manta_T1[,c("EntryLon","EntryLat")], data = gommapps_2017s_manta_T1,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

gommapps_2017s_manta_T1_spt<-spTransform(gommapps_2017s_manta_T1_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(gommapps_2017s_manta_T1_spt,domain)
gommapps_2017s_manta_T1_spt@data$ModelD<-test$DomainFID

gommapps_2017s_manta_T1_spt@data$link<-paste(as.Date(gommapps_2017s_manta_T1_spt@data$Date, "%m/%d/%Y"),"_",gommapps_2017s_manta_T1_spt@data$Transect,"_",gommapps_2017s_manta_T1_spt@data$ModelD,sep="")

gommapps_2018f_manta_T2_spt <- SpatialPointsDataFrame(coords = gommapps_2018f_manta_T2[,c("EntryLon","EntryLat")], data = gommapps_2018f_manta_T2,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

gommapps_2018f_manta_T2_spt<-spTransform(gommapps_2018f_manta_T2_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(gommapps_2018f_manta_T2_spt,domain)
gommapps_2018f_manta_T2_spt@data$ModelD<-test$DomainFID

gommapps_2018f_manta_T2_spt@data$link<-paste(as.Date(gommapps_2018f_manta_T2_spt@data$Date, "%m/%d/%Y"),"_",gommapps_2018f_manta_T2_spt@data$Transect,"_",gommapps_2018f_manta_T2_spt@data$ModelD,sep="")

head(gommapps_2018f_manta_T2_spt@data)

gommapps_2018w_manta_T2_spt <- SpatialPointsDataFrame(coords = gommapps_2018w_manta_T2[,c("EntryLon","EntryLat")], data = gommapps_2018w_manta_T2,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

gommapps_2018w_manta_T2_spt<-spTransform(gommapps_2018w_manta_T2_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(gommapps_2018w_manta_T2_spt,domain)
gommapps_2018w_manta_T2_spt@data$ModelD<-test$DomainFID

gommapps_2018w_manta_T2_spt@data$link<-paste(as.Date(gommapps_2018w_manta_T2_spt@data$Date, "%m/%d/%Y"),"_",gommapps_2018w_manta_T2_spt@data$Transect,"_",gommapps_2018w_manta_T2_spt@data$ModelD,sep="")

gommapps_2017s_manta_T2_spt <- SpatialPointsDataFrame(coords = gommapps_2017s_manta_T2[,c("EntryLon","EntryLat")], data = gommapps_2017s_manta_T2,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

gommapps_2017s_manta_T2_spt<-spTransform(gommapps_2017s_manta_T2_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(gommapps_2017s_manta_T2_spt,domain)
gommapps_2017s_manta_T2_spt@data$ModelD<-test$DomainFID

gommapps_2017s_manta_T2_spt@data$link<-paste(as.Date(gommapps_2017s_manta_T2_spt@data$Date, "%m/%d/%Y"),"_",gommapps_2017s_manta_T2_spt@data$Transect,"_",gommapps_2017s_manta_T2_spt@data$ModelD,sep="")

#### 3.2 GOMMAPPS - Prepare Data for Distance Function #### 
#code received from Gina Rappucci and Lance Garrison on 11/25/2019
#3_Aerial_Extract_Turtle_Sightings_afterFAL2012_25May18_ForT2PosInT1_07Jun18_LPG.r
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552872/

gommapps_mantas_T1<-rbind(gommapps_2018f_manta_T1_spt@data,gommapps_2018w_manta_T1_spt@data,gommapps_2017s_manta_T1_spt@data)
gommapps_mantas_T2<-rbind(gommapps_2018f_manta_T2_spt@data,gommapps_2018w_manta_T2_spt@data,gommapps_2017s_manta_T2_spt@data)

gommapps_mantas_T1$altitude<-183
gommapps_mantas_T2$altitude<-183
gommapps_mantas_T1$distance<-gommapps_mantas_T1$altitude*tan(gommapps_mantas_T1$SightingAngle*pi/180)
gommapps_mantas_T2$distance<-gommapps_mantas_T2$altitude*tan(gommapps_mantas_T2$SightingAngle*pi/180)
head(gommapps_mantas_T1)


#manually adjust the distance by subtracting 3.194, the minimum distance in both T1 and T2, which effectively moves the centerline out
gommapps_mantas_T1$distance<-gommapps_mantas_T1$distance-3.194
gommapps_mantas_T2$distance<-gommapps_mantas_T2$distance-3.194

#Histogram of sightings by observer teams
png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/GOMMAPPS_manta_sighting_distance_hist.png",sep=""),width=7,height=7,res=600,units="in")
par(mfrow=c(2,1))
hist(gommapps_mantas_T1$distance,xlab="Distance (m)",main="Team 1")
hist(gommapps_mantas_T2$distance,xlab="Distance (m)",main="Team 2")
dev.off()

#identify shared sightings
gommapps_mantas_T1$datetime <- as.POSIXlt(paste(as.Date(gommapps_mantas_T1$Date, format = "%m/%d/%Y"), trimws(as.character(format(gommapps_mantas_T1$Time, format = "%I:%M:%S %p")),which=c("right"))),format="%Y-%m-%d %I:%M:%S %p")
gommapps_mantas_T2$datetime <- as.POSIXlt(paste(as.Date(gommapps_mantas_T2$Date, format = "%m/%d/%Y"), trimws(as.character(format(gommapps_mantas_T2$Time, format = "%I:%M:%S %p")),which=c("right"))),format="%Y-%m-%d %I:%M:%S %p")
gommapps$datetime<-as.POSIXlt(gommapps$DATEBEG)

write.csv(gommapps_mantas_T1,"C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/gommapps_mantas_T1.csv")
write.csv(gommapps_mantas_T2,"C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/gommapps_mantas_T2.csv")
Team1<-gommapps_mantas_T1
Team2<-gommapps_mantas_T2

#### Calculate angle for Team1 (Forward Team)
#### and side of plane

Team1$Angle <- 0
Team1$Side <- "U"

n.sights.team1 <- dim(Team1) [1]

for (i in 1:n.sights.team1) {
  if (is.na(Team1$SightingAngle[i]) == FALSE) {
    if (Team1$SightingAngle[i] >= 0) {
      Team1$Angle[i] <- Team1$SightingAngle[i]
    }
    #Handle negative angles
    if (Team1$SightingAngle[i] < 0) {
      Team1$Angle[i] <- -1*Team1$SightingAngle[i]
    }
  }
  
  if (is.na(Team1$SightingAngle[i]) == TRUE) {
    Team1$Angle[i] <- 5 + (Team1$BubbleIncrement[i]-1)*10
  }
  
  if (Team1$OnEffort[i] == 0) {
    Team1$Angle[i] <- NA
  }
  if (Team1$OnEffort[i] == 0) {
    Team1$Side[i] <- "U"
  }
  
  if (Team1$OnEffort[i] == 1) {
    if (Team1$SightingPos[i] == 1) {
      #Handle negative angles   
      if (Team1$SightingAngle[i] >= 0) {
        Team1$Side[i] <- "L"
      }
      
      if (Team1$SightingAngle[i] < 0) {
        Team1$Side[i] <- "R"
      }
    }
    
    if (Team1$SightingPos[i] == 2) {
      #Handle negative angles   
      if (Team1$SightingAngle[i] >= 0) {
        Team1$Side[i] <- "R"
      }
      
      if (Team1$SightingAngle[i] < 0) {
        Team1$Side[i] <- "L"
      }
    }
  }
}


#### Get side and increment from Belly for T1 
#### This is to address T2 positions on T1 table when plane is flying with one team only

Team1$BellyInc <- as.integer(substr(as.character(Team1$BellyIncrement),1,1))
Team1$BellySide <- substr(as.character(Team1$BellyIncrement),2,2)

n.sights.team1 <- dim(Team1) [1]
Team1$ID <- c(1:n.sights.team1)

for (i in 1:n.sights.team1) {
  if (is.na(Team1$SightingAngle[i]) == FALSE) {
    if (Team1$SightingAngle[i] >= 0) {
      Team1$Angle[i] <- Team1$SightingAngle[i]
    }
    #Handle negative angles
    if (Team1$SightingAngle[i] < 0) {
      Team1$Angle[i] <- -1*Team1$SightingAngle[i]
    }
  }
  
  if (is.na(Team1$SightingAngle[i]) == TRUE) {
    
    if(Team1$SightingPosition[i] == 3) {#Belly Only
      Team1$Angle[i] <- 5 + (Team1$BellyInc[i] - 1) * 10 
    }
    
    if(Team1$SightingPosition[i] == 5) {#Belly and Bubble Only
      Team1$Angle[i] <- 5 + (Team1$BubbleInc[i] - 1) * 10 
    }
    
    if(Team1$SightingPosition[i] == 2) {#Bubble Only
      Team1$Angle[i] <- 5 + (Team1$BubbleInc[i] - 1) * 10 
    }
  }
  
  if(Team1$SightingPosition[i] == 5) {
    if (Team1$SightingAngle[i] >= 0) {
      Team1$Side[i] <- "R"
    }
    #Handle negative angles with sight pos 5 
    if (Team1$SightingAngle[i] < 0) {
      Team1$Side[i] <- "L"
    }
  }
  if(Team1$SightingPosition[i] == 2) {
    if(Team1$SightingAngle[i] >= 0) {
      Team1$Side[i] <- "R"
    }
    #Handle negative angles with sight pos 2 
    if (Team1$SightingAngle[i] < 0) {
      Team1$Side[i] <- "L"
    }
  }
  
  if(Team1$SightingPosition[i] == 3) {#Belly Only
    Team1$Side[i] <- Team1$BellySide[i]
  }
  if (Team1$OnEffort[i] == 0) {
    Team1$Angle[i] <- NA
  }
  if (Team1$OnEffort[i] == 0) {
    Team1$Side[i] <- "U"
  }
}


#### Calculate angle for Team2 (Aft Team)
#### and side of plane
#### Get side and increment from Belly
Team2$BellyInc <- as.integer(substr(as.character(Team2$BellyIncrement),1,1))
Team2$BellySide <- substr(as.character(Team2$BellyIncrement),2,2)

n.sights.team2 <- dim(Team2) [1]
Team2$ID <- c(1:n.sights.team2)

for (i in 1:n.sights.team2) {
  if (is.na(Team2$SightingAngle[i]) == FALSE) {
    if (Team2$SightingAngle[i] >= 0) {
      Team2$Angle[i] <- Team2$SightingAngle[i]
    }
    #Handle negative angles
    if (Team2$SightingAngle[i] < 0) {
      Team2$Angle[i] <- -1*Team2$SightingAngle[i]
    }
  }
  
  if (is.na(Team2$SightingAngle[i]) == TRUE) {
    
    if(Team2$SightingPosition[i] == 3) {#Belly Only
      Team2$Angle[i] <- 5 + (Team2$BellyInc[i] - 1) * 10 
    }
    
    if(Team2$SightingPosition[i] == 5) {#Belly and Bubble Only
      Team2$Angle[i] <- 5 + (Team2$BubbleInc[i] - 1) * 10 
    }
    
    if(Team2$SightingPosition[i] == 2) {#Bubble Only
      Team2$Angle[i] <- 5 + (Team2$BubbleInc[i] - 1) * 10 
    }
  }
  
  if(Team2$SightingPosition[i] == 5) {
    if (Team2$SightingAngle[i] >= 0) {
      Team2$Side[i] <- "R"
    }
    #Handle negative angles with sight pos 5 
    if (Team2$SightingAngle[i] < 0) {
      Team2$Side[i] <- "L"
    }
  }
  if(Team2$SightingPosition[i] == 2) {
    if(Team2$SightingAngle[i] >= 0) {
      Team2$Side[i] <- "R"
    }
    #Handle negative angles with sight pos 2 
    if (Team2$SightingAngle[i] < 0) {
      Team2$Side[i] <- "L"
    }
  }
  
  if(Team2$SightingPosition[i] == 3) {#Belly Only
    Team2$Side[i] <- Team2$BellySide[i]
  }
  if (Team2$OnEffort[i] == 0) {
    Team2$Angle[i] <- NA
  }
  if (Team2$OnEffort[i] == 0) {
    Team2$Side[i] <- "U"
  }
}


#### Use the Team1 sightings as the reference
#### step through each one and get Team2 matching sightings
#### match terms are...time < 15 seconds, same side of plane (when known)
#### anglediff <= 15 degrees, and Number of Animals is equal.

n.sight <- dim(Team1) [1]

Team1$MasterIndex<-paste(Team1$datetime,Team1$OtherIndex,sep="_")
Team2$MasterIndex<-paste(Team2$datetime,Team2$OtherIndex,sep="_")

team1.index <- Team1$MasterIndex #instead of OtherIndex because I've combined all surveys
team2.index <- rep(NA, n.sight) #holds matching index

#need to flag whether or not a sighting has been matched already
Team2$matched <- 0

for (i in 1:n.sight) {
  if (Team1$OnEffort[i] == -1) { 
    
    sight.time <- as.POSIXlt(Team1$datetime[i])
    sight.angle <- Team1$Angle[i]
    
    #include groupsize as part of match criteria
    sight.number <- Team1$Number[i]
    
    #deal with angle side
    sight.side <- Team1$Side[i]
    
    time.diff <- abs(difftime(sight.time, Team2$datetime, units = "sec"))
    order <- order(time.diff) 
    time.diff.o <- time.diff[order]
    Team2.o <- Team2[order,]
    Team2.candidates <- Team2.o[time.diff.o < 15,] 
    
    #these are the possible matches, ordered by time.diff
    
    
    #filter out those that have already been matched to something
    Team2.candidates <- Team2.candidates[Team2.candidates$matched == 0,] 
    n.cand <- dim(Team2.candidates) [1]
    
    n.match <- 0
    
    if (n.cand > 0) {
      angle.match <- rep(0, n.cand)
      
      #try to deal with difference in side.  Make all (known) left angles negative for both teams
      #team1 angle - the thing we are trying to match
      if (sight.side == "L") {
        sight.angle.temp <- sight.angle * -1} else {sight.angle.temp <- sight.angle}
      
      #create Angle.Temp in candidates
      Team2.candidates$Angle.Temp <- Team2.candidates$Angle
      Team2.candidates$Angle.Temp[Team2.candidates$Side == "L"] <- Team2.candidates$Angle[Team2.candidates$Side == "L"] * -1
      
      #Handle Unknown side for Team2. Just put them on the same side as the team1 side
      if (sight.side == "L") {
        Team2.candidates$Angle.Temp[Team2.candidates$Side == "U"] <- Team2.candidates$Angle[Team2.candidates$Side == "U"] * -1
      } else {
        Team2.candidates$Angle.Temp[Team2.candidates$Side == "U"] <- Team2.candidates$Angle[Team2.candidates$Side == "U"]
      }
      
      #now use the angle difference.  This will allow match between angles close to the line, but on different sides.  For example 3 L and 3 R will match
      #because angle diff < 15, but 10L and 10R will not.
      
      angle.diff <- abs(sight.angle.temp - Team2.candidates$Angle.Temp)
      angle.match[angle.diff < 15] <- 1
      
      #Turn angle.match back to zero if the number is different.  
      number.diff <- sight.number - Team2.candidates$Number
      angle.match[number.diff != 0] <- 0
      
      Team2.matches <- Team2.candidates[angle.match == 1,]
      n.match <- dim(Team2.matches) [1]
    }
    
    if (n.match > 0) {
      team2.index[i] <- Team2.matches$MasterIndex[1] #will be the first one, since closest in time
      #set matched to 1 so it isn't used again
      Team2$matched[Team2$MasterIndex == Team2.matches$MasterIndex[1]] <- 1
    }
    
    if (n.match == 0) {
      team2.index[i] <- NA
    }
    
  } #end of oneffort
} #end loop

linked.indices <- data.frame(cbind(team1.index, team2.index))
manta.links <- linked.indices[!is.na(linked.indices$team2.index),]

#create fields to hold datetimes with dummy values
my.lt <- strptime("1970-09-12 03:00:02", format="%Y-%m-%d %H:%M:%S")
manta.links$team1.datetime <- my.lt
manta.links$team2.datetime <- my.lt

#and side,angle, and species
side <- "U"
angle <- 0
species <- "Manta"


manta.links$team1.side <- side
manta.links$team1.angle <- angle
manta.links$team1.species <- species
manta.links$team1.number <- 1


manta.links$team2.side <- side
manta.links$team2.angle <- angle
manta.links$team2.species <- species
manta.links$team2.number <- 1

manta.links$side <- "L"
manta.links$species <- "Manta"
manta.links$angle <- 0
manta.links$number <- 1
manta.links$mis.id <- 0

n.link <- dim(manta.links) [1]

for (i in 1:n.link) {
  
  manta.links$team1.datetime[i] <- Team1$datetime[Team1$MasterIndex == manta.links$team1.index[i]]
  manta.links$team2.datetime[i] <- Team2$datetime[Team2$MasterIndex == manta.links$team2.index[i]]
  
  manta.links$team1.side[i] <- Team1$Side[Team1$MasterIndex == manta.links$team1.index[i]]
  manta.links$team1.angle[i] <- Team1$Angle[Team1$MasterIndex == manta.links$team1.index[i]]
  manta.links$team1.species[i] <- as.character(Team1$Description[Team1$MasterIndex == manta.links$team1.index[i]])
  manta.links$team1.number[i] <-Team1$Number[Team1$MasterIndex == manta.links$team1.index[i]]
  
  manta.links$team2.side[i] <- Team2$Side[Team2$MasterIndex == manta.links$team2.index[i]]
  manta.links$team2.angle[i] <- Team2$Angle[Team2$MasterIndex == manta.links$team2.index[i]]
  manta.links$team2.species[i] <- as.character(Team2$Description[Team2$MasterIndex == manta.links$team2.index[i]])
  manta.links$team2.number[i] <- Team2$Number[Team2$MasterIndex == manta.links$team2.index[i]]
  
  manta.links$team1.latitude[i] <- Team1$EntryLat[Team1$MasterIndex == manta.links$team1.index[i]]
  manta.links$team1.longitude[i]<- Team1$EntryLon[Team1$MasterIndex == manta.links$team1.index[i]]
  
  manta.links$team2.latitude[i] <- Team2$EntryLat[Team2$MasterIndex == manta.links$team2.index[i]]
  manta.links$team2.longitude[i] <- Team2$EntryLon[Team2$MasterIndex == manta.links$team2.index[i]]
  
  manta.links$team1.trackindex[i] <- Team1$TrackIndex[Team1$MasterIndex == manta.links$team1.index[i]]
  manta.links$team2.trackindex[i] <- Team2$TrackIndex[Team2$MasterIndex == manta.links$team2.index[i]]
  
  # manta.links$t1.behavior[i] <- as.character(Team1$mantaBehavior[Team1$OtherIndex == manta.links$team1.index[i]])
  # manta.links$t2.behavior[i] <- as.character(Team2$mantaBehavior[Team2$OtherIndex == manta.links$team2.index[i]])
  
  
  
  #create "final" ouput for the link table by combining fields across teams
  #use forward team as reference for side if they are the same side
  
  if (manta.links$team1.side[i] == manta.links$team2.side[i]) {
    manta.links$side[i] <- manta.links$team1.side[i]
    #average the two angles together and round
    manta.links$angle[i] <-  round((manta.links$team1.angle[i] +  manta.links$team2.angle[i])/2,0)
  }
  
  #if they are not the same side - then you have to deal with the angles...set left angles to negative and average..
  if (manta.links$team1.side[i] != manta.links$team2.side[i]) {
    
    if (manta.links$team1.side[i] == "L") {
      temp.team1.angle <- manta.links$team1.angle[i] * -1
    } else {
      temp.team1.angle <- manta.links$team1.angle[i]
    }
    
    if (manta.links$team2.side[i] == "L") {
      temp.team2.angle <- manta.links$team2.angle[i] * -1
    } else {#deal with unknowns
      if (manta.links$team2.side[i] == "R") {
        temp.team2.angle <- manta.links$team2.angle[i]
      } else { #unknowns
        temp.team2.angle <- manta.links$team2.angle[i]
        if (manta.links$team1.side[i] == "L") {temp.team2.angle <- temp.team2.angle *-1}
      }
    }
    
    #average the two angles together and round
    manta.links$angle[i] <-  round((temp.team1.angle +  temp.team2.angle)/2,0)
    
    if (manta.links$angle[i] < 0) {
      manta.links$side[i] <- "L"
    } else {
      manta.links$side[i] <- "R"
    }
    #make the angle positive at the end
    manta.links$angle[i] <- abs(manta.links$angle[i])
    
  }
  
  
  #average the two group sizes
  manta.links$number[i] <- round((manta.links$team1.number[i] +  manta.links$team2.number[i])/2,0)
  # 
  # #use species rules...
  # team1.species <- manta.links$team1.species[i]
  # team2.species <- manta.links$team2.species[i]
  # 
  # #case 1 - they are equal
  # if (team1.species == team2.species) {manta.links$species[i] <- team1.species}
  # 
  # #case 2 - one team has hardshell, the other has an id - set it to the more specific id
  # if (team1.species == "Manta" & team2.species != "Manta") {manta.links$species[i] <- team2.species}
  # if (team2.species == "Manta" & team1.species != "Manta") {manta.links$species[i] <- team1.species}
  # 
  # #case 3 - both teams have an id (neither is hardshell) but it is mismatched
  # if (team1.species != "Manta" & team2.species != "Manta") {
  #   if (team1.species != team2.species) {
  #     manta.links$species[i] <- "Manta"
  #     manta.links$mis.id[i] <- 1}
  # } 
  
}

nrow(linked.indices)
nrow(manta.links)

head(manta.links)
t1sub<-Team1[,c("MasterIndex","ModelD","link")]
t1sub$team1.index<-t1sub$MasterIndex
manta.links<-merge(manta.links,t1sub,by=c("team1.index"))

#manta.links holds only linked records.
#Get the unlinked records for Team1
Team1.nomatch <- Team1[!(Team1$MasterIndex %in% manta.links$team1.index),]
Team1.nomatch$team1.index <- Team1.nomatch$MasterIndex
Team1.nomatch$team1.datetime <- Team1.nomatch$datetime
Team1.nomatch$team1.side <- Team1.nomatch$Side 
Team1.nomatch$team1.angle <- Team1.nomatch$Angle
Team1.nomatch$team1.species <- as.character(Team1.nomatch$Description)
Team1.nomatch$team1.number <- Team1.nomatch$Number
Team1.nomatch$team1.trackindex <- Team1.nomatch$TrackIndex

Team1.nomatch$team1.latitude <- Team1.nomatch$EntryLat
Team1.nomatch$team1.longitude <- Team1.nomatch$EntryLon
#Team1.nomatch$t1.behavior <- as.character(Team1.nomatch$mantaBehavior)

#fill in side, angle, number, species, and mis.id = -1
Team1.nomatch$side <- Team1.nomatch$team1.side
Team1.nomatch$angle <- Team1.nomatch$team1.angle
Team1.nomatch$species <- Team1.nomatch$team1.species
Team1.nomatch$number <- Team1.nomatch$team1.number
Team1.nomatch$mis.id <- -1

keeps <- c("team1.index", "team1.datetime", "team1.side", "team1.angle", "team1.species", "team1.number", "team1.trackindex",
           "team1.latitude", "team1.longitude", "side", "angle", "species", "number", "mis.id","ModelD","link","MasterIndex")
Team1.nomatch <- Team1.nomatch[, names(Team1.nomatch) %in% keeps] 

#fill in empty values for team2 sightings
Team1.nomatch$team2.index <- 0
Team1.nomatch$team2.datetime <- NA
Team1.nomatch$team2.side <- "N"
Team1.nomatch$team2.angle <- NA
Team1.nomatch$team2.species <- "none"
Team1.nomatch$team2.number <- 0
Team1.nomatch$team2.latitude <- NA
Team1.nomatch$team2.longitude <- NA
Team1.nomatch$team2.trackindex <- 0
#Team1.nomatch$t2.behavior <- "none"

rm(all.mantas) #get rid of any previous versions of this data
all.mantas <- rbind(manta.links, Team1.nomatch)

#Get the unlinked records for Team2
Team2.nomatch <- Team2[!(Team2$OtherIndex %in% manta.links$team2.index),]
Team2.nomatch$team2.index <- Team2.nomatch$MasterIndex
Team2.nomatch$team2.datetime <- Team2.nomatch$datetime
Team2.nomatch$team2.side <- Team2.nomatch$Side 
Team2.nomatch$team2.angle <- Team2.nomatch$Angle
Team2.nomatch$team2.species <- as.character(Team2.nomatch$Description)
Team2.nomatch$team2.number <- Team2.nomatch$Number
Team2.nomatch$team2.trackindex <- Team2.nomatch$TrackIndex

Team2.nomatch$team2.latitude <- Team2.nomatch$EntryLat
Team2.nomatch$team2.longitude <- Team2.nomatch$EntryLon
#Team2.nomatch$t2.behavior <- as.character(Team2.nomatch$mantaBehavior)

#fill in side, angle, number, species, and mis.id = -1
Team2.nomatch$side <- Team2.nomatch$team2.side
Team2.nomatch$angle <- Team2.nomatch$team2.angle
Team2.nomatch$species <- Team2.nomatch$team2.species
Team2.nomatch$number <- Team2.nomatch$team2.number
Team2.nomatch$mis.id <- -1

keeps <- c("team2.index", "team2.datetime", "team2.side", "team2.angle", "team2.species", "team2.number", "team2.trackindex",
           "team2.latitude", "team2.longitude", "side", "angle", "species", "number", "mis.id", "ModelD", "link","MasterIndex")

Team2.nomatch <- Team2.nomatch[, names(Team2.nomatch) %in% keeps] 

#fill in empty values for team1 sightings
Team2.nomatch$team1.index <- 0
Team2.nomatch$team1.datetime <- Sys.time()
Team2.nomatch$team1.side <- "N"
Team2.nomatch$team1.angle <- NA
Team2.nomatch$team1.species <- "none"
Team2.nomatch$team1.number <- 0
Team2.nomatch$team1.latitude <- NA
Team2.nomatch$team1.longitude <- NA
Team2.nomatch$team1.trackindex <- 0
#Team2.nomatch$t1.behavior <- "none"

nrow(Team2.nomatch)
#use team2.datetime to look up closest team1.trackindex.
#then can use trackindex to look up environmental variables for every sighting
n.team2 <- dim(Team2.nomatch) [1]
for (i in 1:n.team2) {
  date.time <- Team2.nomatch$team2.datetime[i]
  survey.track.recs <- surveytrack[(surveytrack$datetime >= date.time - 20) & (surveytrack$datetime <= date.time + 20),]
  survey.track.recs$time.diff <- abs(date.time - survey.track.recs$datetime)
  survey.track.index <- survey.track.recs$Index[survey.track.recs$time.diff == min(survey.track.recs$time.diff)] [1] #in case there is a tie
  
  Team2.nomatch$team1.trackindex[i] <- survey.track.recs$Index[survey.track.recs$Index == survey.track.index]
  Team2.nomatch$team1.datetime[i] <- survey.track.recs$datetime[survey.track.recs$Index == survey.track.index]
  Team2.nomatch$team1.latitude[i] <- survey.track.recs$Latitude[survey.track.recs$Index == survey.track.index]
  Team2.nomatch$team1.longitude[i] <- survey.track.recs$Longitude[survey.track.recs$Index == survey.track.index]
}

all.mantas <- rbind(all.mantas, Team2.nomatch)

#add unique id
unique.id <- c(1:(dim(all.mantas) [1]))
all.mantas$unique.id <- unique.id

#add avail.aft - eq to 0 if side = L and angle > 40
all.mantas$avail.aft <- 0
all.mantas$avail.aft[all.mantas$angle <= 40] <- 1  #anything with angle < 40
all.mantas$avail.aft[all.mantas$team2.angle >= 0] <- 1 #anything actually seen by team 2
all.mantas$avail.aft[all.mantas$angle > 40 & all.mantas$side == "R"] <- 1

#INCLUDE CUSTOM CHANGES IN AVAILABILITY - EG ONE TEAM NOT ACTIVE - HERE
# e.g.: all.mantas$avail.aft[all.mantas$team1.index == 598] <- 0
#create history based on angles
all.mantas$history <- "00"
all.mantas$history[(!is.na(all.mantas$team1.angle) & is.na(all.mantas$team2.angle))] <- "10"
all.mantas$history[(!is.na(all.mantas$team1.angle) & !is.na(all.mantas$team2.angle))] <- "11"
all.mantas$history[(is.na(all.mantas$team1.angle) & !is.na(all.mantas$team2.angle))] <- "01"

#location and time variables based on team1 survey track - since this links to the effort files.
all.mantas$lat <- all.mantas$team1.latitude
all.mantas$lon <- all.mantas$team1.longitude
all.mantas$surveytrack.index <- all.mantas$team1.trackindex
all.mantas$init.datetime <- all.mantas$team1.datetime

#create effort status variables
all.mantas$t1.effort <- 0
all.mantas$t2.effort <- 0
all.mantas$t1.effort[!is.na(all.mantas$team1.angle)] <- 1
all.mantas$t2.effort[!is.na(all.mantas$team2.angle)] <- 1

#calculate distance based on average of angles
all.mantas$dist <- 183 * tan(all.mantas$angle * pi/180)

#merge in environmental variables using surveytrack.index
mantacopy<-all.mantas
mantacopy$team1.datetime<-as.POSIXct(mantacopy$team1.datetime)
mantacopy$team2.datetime<-as.POSIXct(mantacopy$team2.datetime)
mantacopy$init.datetime<-as.POSIXct(mantacopy$init.datetime)
mantacopy$Date<-as.Date(mantacopy$team1.datetime)
surveytrack2<-surveytrack
surveytrack2$datetime<-as.POSIXct(surveytrack2$datetime)
surveytrack2$Date<-as.Date(surveytrack2$datetime)

# head(mantacopy)
# all.mantas.env <- merge(x=all.mantas, y=surveytrack, by.x = "surveytrack.index", by.y = "Index",all.x=T,all.y=F)
all.mantas.env <- left_join(mantacopy, unique(surveytrack2),by=c("surveytrack.index"="Index","Date"="Date"))
all.mantas.env <- all.mantas.env[order(all.mantas.env$unique.id),]

#clean up
drops <- c("team1.datetime", "team2.datetime", "team1.side", "team1.angle", "team1.number", "team2.side", "team2.angle", "team2.number", 
           "team1.latitude", "team1.longitude", "team2.latitude", "team2.longitude", "team1.trackindex", "team2.trackindex", "RecType", 
           "Survey", "Date", "Time", "Latitude", "Longitude", "Speed", "Heading", "WaterTemp", "LeftObs", "RightObs", "BellyObs", "Recorder", 
           "OnEffort", "OnTransect", "datetime")

all.mantas.f <- all.mantas.env[,!(names(all.mantas.env) %in% drops)]
head(all.mantas.f)
#rename columns
names(all.mantas.f) [1] <- "t1.index"
names(all.mantas.f) [2] <- "t2.index"
names(all.mantas.f) [3] <- "t1.species"
names(all.mantas.f) [4] <- "t2.species"
names(all.mantas.f) [6] <- "CommonName"
names(all.mantas.f) [8] <- "groupsize"
names(all.mantas.f) [37] <- "track"
names(all.mantas.f)

write.csv(all.mantas.f, "C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GoMMAPPS/all_mantas.csv", row.names = FALSE)


#### 3.3 GOMMAPPS - Run distance analysis (stratified) ####
#based on "Aerial_Stratified_TwoTeam_Distance_Analysis.r" from LPG
effort.init <- gommapps
effort.init$effort<-as.numeric(effort.init$EFFORTTY)
summary(effort.init)
#only on effort parts
effort.init <- effort.init[effort.init$effort == 1,]
effort.init$stratid<-"default"
#aggregate by stratum and transect
effort <- aggregate(effort.init$segmentLength_m, by=list(Region.Label=effort.init$stratid,
                                                         Sample.Label=effort.init$TRANID), FUN=sum) #effort.init$ FID_ModelD
effort$Effort <- effort$x
effort$x <- NULL
effort <- effort[effort$Effort >= 1,] #drop short segments...

summary(all.mantas.f)

sightings<-all.mantas.f
sightings$Region.Label <- "default" #experiment with stratification by latitude later?
sightings$Sample.Label <- sightings$ModelD
sightings$distance <- sightings$dist
sightings$size <- sightings$groupsize
sightings$object <- sightings$unique.id

#drop sightings coded as off effort for both teams
sightings.f <- sightings[!(sightings$history== "00"),]

#eliminate any sightings where aft team could not have seen the manta
sightings.f <- sightings.f[(sightings.f$avail.aft==1),]

#Create the two teams
sightings.t1 <- sightings.f
sightings.t1$observer <- 1
sightings.t1$detected <- 0
sightings.t1$detected[sightings.t1$t1.effort == 1] <- 1

sightings.t2 <- sightings.f
sightings.t2$observer <- 2
sightings.t2$detected <- 0
sightings.t2$detected[sightings.t2$t2.effort == 1] <- 1

sightings.distance <- rbind(sightings.t1, sightings.t2)
sightings.distance <- sightings.distance[order(sightings.distance$object, sightings.distance$observer),]

hist(sightings.distance$distance, breaks = 15, col = "gray", main = "Sighting Distances", xlab = "PSD(m)", ylab = "Number of Sightings")

par(mfrow = c(2,1))
hist(sightings.distance$distance[sightings.distance$observer == 1 & sightings.distance$detected == 1], breaks = 15, col = "gray", main = "Team 1", xlab = "PSD(m)", ylab = "Number of Sightings",xlim=c(0,325))
hist(sightings.distance$distance[sightings.distance$observer == 2 & sightings.distance$detected == 1], breaks = 15, col = "gray", main = "Team 2", xlab = "PSD(m)", ylab = "Number of Sightings",xlim=c(0,325))

#right truncate at 300
RT.distance <- 300

cds.null.hr <-ddf(method='ds',dsmodel=~mcds(key='hr', formula = ~ 1),
                  data=sightings.distance, 
                  meta.data=list(binned=F, width = RT.distance))

cds.null.hn <-ddf(method='ds',dsmodel=~mcds(key='hn', formula = ~ 1),
                  data=sightings.distance, 
                  meta.data=list(binned=F, width = RT.distance))

plot(cds.null.hr)
plot(cds.null.hn)
summary(cds.null.hr)
ddf.gof(cds.null.hr)

summary(cds.null.hn)
ddf.gof(cds.null.hn)

strata<-data.frame("Region.Label"="default","Area"=area(domain)) ###NAF check here that area computed in square meters

#Compute density and abundance estimates and variances based on Horvitz-Thompson-like estimator.
dht.null <- dht(cds.null.hr, region.table=strata, sample.table=effort, se = TRUE)#, options = list(convert.units=0.001))
dht.null

#Generate density, n, and cv

dht.spe <- dht(cds.null.hr, strata, effort, sightings.distance, se = TRUE)#, options=list(convert.units=0.001))
n.est <- dht.spe$individuals$N$Estimate
d.est <- dht.spe$individuals$D$Estimate
cv.est <- dht.spe$individuals$N$cv

#Abundance estimates by species
estimates.out <- data.frame(D = d.est[1], N = n.est[1], cv = cv.est[1])
estimates.out

#MODEL SELECTION - STEP 1:  HR vs HN on null dsmodel and distance only mr model 
#- select based on min AIC/GOF tests

mrds.null.hr <-ddf(method='io',dsmodel=~mcds(key='hr', formula = ~ 1),
                   mrmodel=~glm(link='logit',formula=~distance), data=sightings.distance, 
                   meta.data=list(binned=F, width = RT.distance))

mrds.null.hn <-ddf(method='io',dsmodel=~mcds(key='hn', formula = ~ 1),
                   mrmodel=~glm(link='logit',formula=~distance), data=sightings.distance, 
                   meta.data=list(binned=F, width = RT.distance))


plot(mrds.null.hn)
plot(mrds.null.hr)

ddf.gof(mrds.null.hn)
ddf.gof(mrds.null.hr)
summary(mrds.null.hn)
summary(mrds.null.hr)

ds.key <- "hr"

#look at covariates
sightings.distance$g<-ifelse(sightings.distance$side=="R",sightings.distance$RightGlareCover,sightings.distance$LeftGlareCover)
sightings.distance$gi<-ifelse(sightings.distance$side=="R",sightings.distance$RightGlareIntens,sightings.distance$LeftGlareIntens)

#unique sightings
sight.unique <- sightings.distance[sightings.distance$observer == 1,]
e.names <- c("SeaState","SunPen","Weather","g","gi","Turbidity","CloudCover","Conditions","Haze","Fog","loggs")
num.env <- length(e.names)
sight.unique$loggs <- log(sight.unique$size)
head(sight.unique)

windows()
par(mfrow = c(3,4))
plot(sight.unique$SeaState, sight.unique$distance, ylab = "PSD", xlab = "Sea State")
plot(sight.unique$SunPen, sight.unique$distance, ylab = "PSD", xlab = "Sun Penetration")
plot(sight.unique$Weather, sight.unique$distance, ylab = "PSD", xlab = "Weather")
plot(sight.unique$g, sight.unique$distance, ylab = "PSD", xlab = "Glare Cover")
plot(sight.unique$gi, sight.unique$distance, ylab = "PSD", xlab = "Glare Intensity")
plot(sight.unique$Turbidity, sight.unique$distance, ylab = "PSD", xlab = "Turbidity")
plot(sight.unique$CloudCover, sight.unique$distance, ylab = "PSD", xlab = "Cloud Cover")
plot(sight.unique$Conditions, sight.unique$distance, ylab = "PSD", xlab = "Conditions")
plot(sight.unique$Haze, sight.unique$distance, ylab = "PSD", xlab = "Haze")
plot(sight.unique$Fog, sight.unique$distance, ylab = "PSD", xlab = "Fog")
plot(sight.unique$loggs, sight.unique$distance, ylab = "PSD", xlab = "lgGS")

#recode variables
#Sea state has pattern, sun penetration not obvious pattern, weather no contrast, glare cover insufficient data
#Turbidity no pattern, Cloud Cover undecided, Conditions no clear pattern, Haze no clear pattern, Fog no clear pattern
sight<-sightings.distance
sight$SeaState2 <- sight$SeaState
sight$SeaState2[sight$SeaState >= 3] <- 3 #combine ss 3-5

#test correlation with log(Size) - include log(size) in DS model if p < 0.1
size.cor <- cor.test(sight.unique$distance, log(sight.unique$size)) 
size.cor
# p-value = 0.3563, so drop

n.var <- 9
##create formulas for ds
DF <- data.frame(Class=1:n.var, SeaState=1:n.var, SunPen=1:n.var, Haze=1:n.var, 
                 g=1:n.var, gi=1:n.var, Turbidity=1:n.var, CloudCover=1:n.var,
                 Conditions=1:n.var)
Cols <- names(DF)
Cols <- Cols[! Cols %in% "Class"]
n <- length(Cols)
id<-unlist(lapply(1:n,function(i)combn(1:n,i,simplify=F)),recursive=F)
Formulas<-sapply(id,function(i)paste("~",paste(Cols[i],collapse="+")))

#Add the null models
Formulas[length(Formulas) + 1] <- "~1"

##Loop through all formulas for the ds and save the AIC
ds.AIC <-vector("numeric",length(Formulas))


for (i in 1:length(Formulas)){
  try({
    mcds.test <-ddf(method="ds",dsmodel=~mcds(key=ds.key, formula = Formulas[i]), data= sight, 
                    meta.data=list(binned=F, width = RT.distance))
    ds.AIC[i]<-mcds.test$criterion	
  }, silent = FALSE)
}

ds.AIC.out <- data.frame(model = Formulas, AIC = ds.AIC)
ds.AIC.out<-ds.AIC.out[order(ds.AIC.out$AIC),]
ds.AIC.out<-ds.AIC.out[(ds.AIC.out$AIC > 0),]
head(ds.AIC.out)
# model      AIC
# 234              ~ SeaState+Haze+g+gi+Turbidity+CloudCover 1149.399
# 80                             ~ Haze+Turbidity+CloudCover 1151.648
# 252 ~ SeaState+SunPen+g+gi+Turbidity+CloudCover+Conditions 1156.619
# 158                            ~ g+gi+Turbidity+CloudCover 1156.878
# 133                            ~ SunPen+Haze+gi+CloudCover 1159.873
# 159                            ~ g+gi+Turbidity+Conditions 1161.225

model<-models[[2]]
temp<-c(model$ds$aux$ddfobj$type,
        model$ddf$ds$aux$ddfobj$scale$formula,
        model$ddf$criterion,
        ddf.gof(model, qq=FALSE)$dsgof$CvM$p,
        summary(model)$average.p,
        summary(model)$average.p.se
)
##find the best AIC and the matching formula
#https://workshops.distancesampling.org/duke-spatial-2015/practicals/1-detection-functions-solutions.html
make_table <- function(models){
  
  # this function extracts the model data for a single model (row)
  extract_model_data <- function(model){
    c(model$ds$aux$ddfobj$type,
      model$ds$aux$ddfobj$scale$formula,
      model$criterion,
      ddf.gof(model, qq=FALSE)$dsgof$CvM$p,
      summary(model)$average.p,
      summary(model)$average.p.se
    )
  }
  
  temp<-list()
  for(i in 1:length(models)) {
    if(!is.na(models[i])) {
      model<-models[[i]]
      temp[i]<-as.data.frame(extract_model_data(model),stringsAsFactors=FALSE)
    }
  }
  
  res<-as.data.frame(do.call("rbind", temp),stringsAsFactors=FALSE)
  
  # 
  # 
  # # applying that to all the models then putting it into a data.frame
  # res <- as.data.frame(t(as.data.frame(lapply(na.omit(models), extract_model_data))),
  #                      stringsAsFactors=FALSE)
  
  # making sure the correct columns are numeric
  res[,3] <- as.numeric(res[,3])
  res[,4] <- as.numeric(res[,4])
  res[,5] <- as.numeric(res[,5])
  res[,6] <- as.numeric(res[,6])
  
  # giving the columns names
  colnames(res) <- c("Key function", "Formula", "AIC", "Cramer-von Mises $p$-value",
                     "$\\hat{P_a}$", "se($\\hat{P_a}$)")
  
  # creating a new column for the AIC difference to the best model
  res[["$\\Delta$AIC"]] <- res$AIC - min(res$AIC, na.rm=TRUE)
  # ordering the model by AIC score
  res <- res[order(res$AIC),]
  
  # returning the data.frame
  return(res)
}

mcds.formula <- as.character(ds.AIC.out$model[ds.AIC.out$AIC == min(ds.AIC.out$AIC)])
#"~ SeaState+Haze+g+gi+Turbidity+CloudCover"

m.formula<-c(8,1)
for(i in 1:8) {
  m.formula[i]<-as.character(ds.AIC.out[i,"model"])
}

models<-list()
for(i in 1:8) {
  models[[i]] = tryCatch({
    ddf(method="ds",dsmodel=~mcds(key=ds.key, formula = m.formula[i]), data= sight, 
        meta.data=list(binned=F, width = RT.distance))
  }, warning = function(w) {
    NA  }, error = function(e) {
      NA  })
}

model_table <- make_table(models)
kable(model_table, digits=3)
model_table

# Key function                                                           Formula      AIC Cramer-von Mises $p$-value $\\hat{P_a}$ se($\\hat{P_a}$)
# 1           hr                                    ~Haze + Turbidity + CloudCover 1151.649               0.0663490707    0.6332310       0.14571894
# 2           hr ~SeaState + SunPen + g + gi + Turbidity + CloudCover + Conditions 1156.619               0.0363421596    0.6174853       0.06425632
# 4           hr                                  ~g + gi + Turbidity + Conditions 1161.225               0.0105741157    0.6579777       0.06271257
# 5           hr                         ~SeaState + Haze + Turbidity + CloudCover 1162.488               0.0926436752    0.6251197       0.02348093
# 6           hr                                     ~SunPen + g + gi + CloudCover 1164.079               0.0210283554    0.6460645       0.02098960
# 3           hr                                  ~SunPen + Haze + gi + CloudCover 1244.013               0.0001208716    1.0000000       0.00000000
# $\\Delta$AIC
# 1     0.000000
# 2     4.970491
# 4     9.576316
# 5    10.838988
# 6    12.430588
# 3    92.364481

test<-ddf(method="ds",dsmodel=~mcds(key=ds.key, formula = model_table[1,"Formula"]), data= sight, 
          meta.data=list(binned=F, width = RT.distance))
ddf.gof(test)

#restrict responses to models with CVM p-value > 0.05, indicating goodness of fit
#https://workshops.distancesampling.org/online-course/exercisepdfs/Ch7/E7-1-ducknests-sol.pdf 
valid_dsfits<-model_table[model_table$`Cramer-von Mises $p$-value`>0.05,]
mcds.formula<-valid_dsfits[1,"Formula"]
mcds.best<-ddf(method="ds",dsmodel=~mcds(key=ds.key, formula = mcds.formula), data= sight, 
               meta.data=list(binned=F, width = RT.distance))
summary(mcds.best)
plot(mcds.best,main=mcds.formula)
ddf.gof(mcds.best)

#generate formulas for mr model. Include distance in all models, exclude size
n.var <- 9
##create formulas for ds
DF <- data.frame(Class=1:n.var, SeaState=1:n.var, SunPen=1:n.var, Haze=1:n.var, 
                 g=1:n.var, gi=1:n.var, Turbidity=1:n.var, CloudCover=1:n.var,
                 Conditions=1:n.var)
Cols <- names(DF)
Cols <- Cols[! Cols %in% "Class"]
n <- length(Cols)
id<-unlist(lapply(1:n,function(i)combn(1:n,i,simplify=F)),recursive=F)
Formulas<-sapply(id,function(i)paste("~distance +",paste(Cols[i],collapse="+")))
Formulas2<-sapply(id,function(i)paste("~distance * observer +",paste(Cols[i],collapse="+")))

#Add the null models
Formulas[length(Formulas) + 1] <- "~1"
Formulas[length(Formulas) + 1] <- "~distance * observer"

mrds.Formulas <- c(Formulas, Formulas2)

##Loop through all formulas for the ds and save the AIC
mrds.AIC <-vector("numeric",length(mrds.Formulas))

for (i in 1:length(mrds.Formulas)){
  try({
    mrds.test <-ddf(method="io",
                    dsmodel=~mcds(key=ds.key, formula = mcds.formula),
                    mrmodel=~glm(link='logit', formula= mrds.Formulas[i]),
                    data= sight, 
                    meta.data=list(binned=F, width = RT.distance))
    mrds.AIC[i]<-mrds.test$criterion
  }, silent = TRUE)		
  print(i)
}

mrds.AIC.out <- data.frame(model = mrds.Formulas, AIC = mrds.AIC)
mrds.AIC.out<-mrds.AIC.out[order(mrds.AIC.out$AIC),]
mrds.AIC.out<-mrds.AIC.out[which(mrds.AIC.out$AIC > 0),]
head(mrds.AIC.out)

mrds.formula <- as.character(mrds.AIC.out$model[mrds.AIC.out$AIC == min(mrds.AIC[mrds.AIC > 0])])

#Set best model
mrds.best <-ddf(method='io',dsmodel=~mcds(key= ds.key, formula = mcds.formula),
                mrmodel=~glm(link='logit',formula= mrds.formula), data=sight, 
                meta.data=list(binned=F, width = RT.distance))

summary(mrds.best)
ddf.gof(mrds.best,main="no interaction")

mrds.best.nullds <-ddf(method='io',dsmodel=~mcds(key= ds.key, formula = ~1),
                       mrmodel=~glm(link='logit',formula= ~distance * observer), data=sight, 
                       meta.data=list(binned=F, width = RT.distance))


summary(mrds.best.nullds)
ddf.gof(mrds.best.nullds)

mrds.best.int <-ddf(method='io',dsmodel=~mcds(key= ds.key, formula = mcds.formula),
                    mrmodel=~glm(link='logit',formula= ~distance * observer), data=sight, 
                    meta.data=list(binned=F, width = RT.distance))


summary(mrds.best.int)
ddf.gof(mrds.best.int,main="interaction")

windows()
par(mfrow = c(3,2))
plot(mrds.best)

windows()
par(mfrow = c(3,2))
plot(mrds.best.nullds)

windows()
par(mfrow = c(3,2))
plot(mrds.best.int)

mrds.best$criterion
mrds.best.nullds$criterion
mrds.best.int$criterion

dht <- dht(mrds.best, strata, effort, se = TRUE)#, options = list(convert.units=0.001))
dht
n.est <- dht$individuals$N$Estimate
d.est <- dht$individuals$D$Estimate
cv.est <- dht$individuals$N$cv

#Abundance estimates by species
estimates.out <- data.frame(cbind(D = d.est, N = n.est, cv = cv.est))
estimates.out

#### 3.4 GOMMAPPS - Determine Effective Strip Width ####
###next steps: http://distancesampling.org/R/vignettes/mexico-analysis.html

#3.4.1. filter for unique manta sightings (across both teams)
positives<-sight[sight$detected==1,]
uniquepos<-positives[!duplicated(positives$object),]

# manta_spdf<-SpatialPointsDataFrame(coords = cbind(uniquepos$lon,uniquepos$lat), data = uniquepos, proj4string = crs(wgs.84))
# manta_spdf<-spTransform(manta_spdf,crs(gommapps))
# windows()
# plot(gommapps)
# plot(manta_spdf,pch=3,cex=3,col="red",add=T)
# 
# #points need to become small polygons so they have area and can "intersect" the trackline
# manta_spdf.poly<-rgeos::gBuffer(manta_spdf,byid=T,width=0.001)

#3.4.2 link unique manta observations to the tracklines
surveytrack$link2<-paste(surveytrack$datetime,"_",surveytrack$Index,sep="")
surveytrack$surveytrack.index<-surveytrack$Index
surveytrack$SeaState2 <- surveytrack$SeaState
surveytrack$SeaState2[surveytrack$SeaState >= 3] <- 3 #combine ss 3-5
# uniquepos$link2<-paste(uniquepos$init.datetime,"_",uniquepos$surveytrack.index,sep="")
# track.predict<-merge(uniquepos,surveytrack,by=c("link2","surveytrack.index","GPSAltitude","SurveyAltitude","SeaState", "CloudCover",
#                                                 "LeftGlareCover","RightGlareCover","LeftGlareIntens","RightGlareIntens", 
#                                                 "Weather","Conditions","Haze","Fog","SunPen","Turbidity","SeaState2","AuditingCode"),
#                      all.y=T)
# nrow(track.predict) #more rows because some segments had multiple (distinct) mantas observed
# nrow(surveytrack)
# 
# # test<-over(gommapps_2017s_manta_T2_spt,domain)
# # gommapps_2017s_manta_T2_spt@data$ModelD<-test$DomainFID
# 
# 
# #3.4.3 fit predicted detection probabilities to the tracklines (try to retain SE in detection probability to allow bootstrapping later)
# track.predict$distance<-ifelse(is.na(track.predict$distance),RT.distance,track.predict$distance)
# track.predict$detected<-ifelse(is.na(track.predict$detected),1,track.predict$detected)
# track.predict$size<-ifelse(is.na(track.predict$size),0,track.predict$size)
# track.predict$observer<-ifelse(is.na(track.predict$observer),1,track.predict$observer)
# track.predict$object<-ifelse(is.na(track.predict$object),0,track.predict$object)
# track.predict$Region.Label<-"default"
# track.predict$Sample.Label<-0
# 
# #compute ESW for the right side of the plane
# track.predict$g<-track.predict$RightGlareCover
# track.predict$gi<-track.predict$RightGlareIntens
# track.predict$CommonName<-"Manta Ray"
# head(track.predict)
# test<-predict(mrds.best,sight,compute=T,esw=T)
# head(test$fitted)
# head(sight)
# eswr<-predict(mrds.best,newdata=track.predict,compute=T,esw=T)
# head(eswr$fitted)
# length(eswr$fitted)

surveytrack$observer<-1
surveytrack2<-surveytrack
surveytrack2$observer<-2

surveytrack.double<-rbind(surveytrack,surveytrack2)
sight2<-sight
sight2$link2<-paste(sight2$init.datetime,sight2$surveytrack.index,sight2$observer,sep="_")
surveytrack.double$link2<-paste(surveytrack.double$datetime,surveytrack.double$Index,surveytrack.double$observer,sep="_")
track.predict2<-merge(sight2,surveytrack.double,by=c("link2","surveytrack.index","GPSAltitude","SurveyAltitude","SeaState", "CloudCover",
                                                     "LeftGlareCover","RightGlareCover","LeftGlareIntens","RightGlareIntens", "observer",
                                                     "Weather","Conditions","Haze","Fog","SunPen","Turbidity","SeaState2","AuditingCode"),
                      all.y=T)

nrow(surveytrack.double)
nrow(track.predict2) #note there are extra rows here because sometimes Team 1 observed an individual and Team 2 observed a different individual
track.predict2$distance<-ifelse(is.na(track.predict2$distance),RT.distance,track.predict2$distance)
#track.predict2$distance<-ifelse(track.predict2$distance==0,RT.distance,track.predict2$distance) #NAF check this step
track.predict2$detected<-ifelse(is.na(track.predict2$detected),1,track.predict2$detected)
track.predict2$size<-ifelse(is.na(track.predict2$size),0,track.predict2$size)
track.predict2$observer<-ifelse(is.na(track.predict2$observer),1,track.predict2$observer)
track.predict2$object<-ifelse(is.na(track.predict2$object),0,track.predict2$object)
track.predict2$Region.Label<-"default"
track.predict2$Sample.Label<-0

#3.4.3 - compute ESW
#compute ESW for the right side of the plane
track.predict2$g<-track.predict2$RightGlareCover
track.predict2$gi<-track.predict2$RightGlareIntens
track.predict2$CommonName<-"Manta Ray"
detp_r<-predict(mrds.best,newdata=track.predict2,compute=F,esw=F)
#compute ESW for the left side of the plane
track.predict2$g<-track.predict2$LeftGlareCover
track.predict2$gi<-track.predict2$LeftGlareIntens
detp_l<-predict(mrds.best,newdata=track.predict2,compute=F,esw=F)

gommapps.detp<-data.frame(detp_r$fitted,detp_l$fitted)
gommapps.detp$esw<-(gommapps.detp$detp_r.fitted+gommapps.detp$detp_r.fitted)*RT.distance #NAF added because esw was just producing detection probs
summary(gommapps.detp)
track.predict2<-cbind(track.predict2,gommapps.detp)
head(track.predict2)
write.csv(track.predict2,"C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort/track_predict2.csv")

#3.4.4 assign environmental parameters to tracklines
#performed this step in "assign_satellite_data_sefsc_survey.R"

#3.4.5 compute segment length for each trackline segment
track_predict2_env<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort", layer="track_predict2_env") # this is in equal area projection
#tail(track_predict2_env@data) #error in reading in the data
gommapps.env<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort/track_predict2_env3.csv",stringsAsFactors = F)
#Remove records for transit with buggy times in the "SurveyTrack_T1_FOR_2018W" data (index 1:11)
gommapps.env$Date<-as.Date(gommapps.env$Date,format="%Y-%m-%d")
d.gommapps<-gommapps.env[-(which(gommapps.env$Date==as.Date("2018-01-18") & gommapps.env$surveytrack.index <= 11 )),]

#3.4.6 compute strip area (LPG suggests length of the segment * 2 * truncation distance; NAF used length * esw)
d.gommapps$datetime<-as.POSIXct(d.gommapps$datetime)
#temp<-track_predict2_env@data[order(track_predict2_env@data[,"Date"],track_predict2_env@data[,"surveytrack.index"],track_predict2_env@data[,"datetime"]),]

d.gommapps$striplength<-NA
for(i in 1:(nrow(d.gommapps)-2)) { #data steps through 2 observations per timestep
  t1<-d.gommapps[i,"datetime"]
  t2<-d.gommapps[(i+2),"datetime"]
  if(difftime(t2,t1,units="secs")<120) { #valid next point
    pt1<-SpatialPoints(cbind((d.gommapps[i,"Longitude"]),d.gommapps[i,"Latitude"])) #identify the current point
    pt2<-SpatialPoints(cbind((d.gommapps[(i+2),"Longitude"]),d.gommapps[(i+2),"Latitude"])) #identify the next point
    d.gommapps[i,"striplength"]<-pointDistance(pt1,pt2,lonlat=T) #distance in meters between points
  }
}
#temp<-d.gommapps[which(d.gommapps$striplength>2000),]; View(temp) #verified no absurd speeds
d.gommapps$striparea<-d.gommapps$esw*d.gommapps$striplength #in meters squared NAF verify units look correct
summary(d.gommapps)

#3.4.7 aggregate data within larger grids
d.gommapps_spt <- SpatialPointsDataFrame(coords = d.gommapps[,c("Longitude","Latitude")], data = d.gommapps,
                                         proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

d.gommapps_spt<-spTransform(d.gommapps_spt, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test<-over(d.gommapps_spt,domain)
d.gommapps_spt@data$DomainFID<-test$DomainFID
summary(d.gommapps_spt@data)

DomainFID=summarize(d.gommapps_spt@data$striparea,d.gommapps_spt@data$DomainFID,sum,na.rm=T)[1]
striparea=summarize(d.gommapps_spt@data$striparea,d.gommapps_spt@data$DomainFID,sum,na.rm=T)
N=summarize(d.gommapps_spt@data$size,d.gommapps_spt@data$DomainFID,sum,na.rm=T)
Depth_m=summarize(d.gommapps_spt@data$Depth_m,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
Front=summarize(d.gommapps_spt@data$Front,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
Front_Z=summarize(d.gommapps_spt@data$Front_Z,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
SST=summarize(d.gommapps_spt@data$SST,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
Chla=summarize(d.gommapps_spt@data$ChlA,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
waterv=summarize(d.gommapps_spt@data$waterv,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
oscarv=summarize(d.gommapps_spt@data$oscarv,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
pp=summarize(d.gommapps_spt@data$pp,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
k490=summarize(d.gommapps_spt@data$k490,d.gommapps_spt@data$DomainFID,mean,na.rm=T)
wave=summarize(d.gommapps_spt@data$wave,d.gommapps_spt@data$DomainFID,mean,na.rm=T)

gommapps_grid<-Reduce(function(x,y) merge(x = x, y = y, by = "d.gommapps_spt@data$DomainFID"), 
                      list(striparea,N,Depth_m,Front,Front_Z,SST,Chla,waterv,oscarv,pp,k490,wave))

colnames(gommapps_grid)<-c("DomainFID","striparea","N","Depth_m","Front","Front_Z","SST","Chla","waterv","oscarv","pp","k490","wave")
gommapps_grid$Nhat<-gommapps_grid$N/gommapps_grid$striparea
gommapps_grid<-gommapps_grid[which(gommapps_grid$Depth_m<0), ] #only ocean obs.
gommapps_grid<-remove_all_labels(gommapps_grid)
gommapps_grid$stripareaKM<-gommapps_grid$striparea/(1000*1000)
summary(gommapps_grid)

#### 3.5 GOMMAPPS - fit GAM ####
#generate formulas for gam model
cor(gommapps_grid$waterv,gommapps_grid$oscarv, use = "pairwise.complete.obs") #-0.05
cor(gommapps_grid$k490,gommapps_grid$pp, use = "pairwise.complete.obs") #0.73
cor(gommapps_grid$k490,gommapps_grid$Chla, use = "pairwise.complete.obs") #0.86
plot(gommapps_grid$k490~gommapps_grid$pp)
#use k490 as most complete proxy for Chla and pp
plot(gommapps_grid$waterv~gommapps_grid$oscarv)
#inclusion of wave, waterv and oscarv results in lots of NAs which causes issues for model comparisons; eliminate at this time
#eliminate wave assuming the sea conditions portion of the distance-weighted sampling handles that a bit

n.var <- 10
# DF <- data.frame(Class=1:n.var, Front_Z=1:n.var, SST=1:n.var, ChlA=1:n.var, waterv=1:n.var, 
#                  oscarv=1:n.var, pp=1:n.var, k490=1:n.var, wave=1:n.var,
#                  Depth_m=1:n.var)
DF <- data.frame(Class=1:n.var, Front=1:n.var, Front_Z=1:n.var, SST=1:n.var, k490=1:n.var,
                 Depth_m=1:n.var)
Cols <- names(DF)
Cols <- Cols[! Cols %in% "Class"]
n <- length(Cols)
id<-unlist(lapply(1:n,function(i)combn(1:n,i,simplify=F)),recursive=F)
#NAF note k=4 resulted in artificial "peakiness" at right end of SST and k490 splines; limiting all splines to 3 knots
Formulas<-sapply(id,function(i)paste("N~offset(log(striparea))+",paste('ti(',Cols[i],',k=3,bs="ts")',collapse="+",sep="")))
#+ti(Front_Z,SST,k=c(4,4),bs=c("ts","ts"))
Formulas2<-sapply(id,function(i)paste("N~offset(log(striparea))+ti(Front_Z,SST,k=c(3,3),bs=c('ts','ts'))+",paste('ti(',Cols[i],',k=3,bs="ts")',collapse="+",sep="")))
Formulas3<-sapply(id,function(i)paste("N~offset(log(striparea))+ti(Front_Z,SST,k=c(3,3),bs=c('ts','ts'))+ti(Front_Z,SST,k490,k=c(3,3,3),bs=c('ts','ts','ts'))+",paste('ti(',Cols[i],',k=3,bs="ts")',collapse="+",sep="")))
Formulas4<-c(Formulas,Formulas2,Formulas3)
#Add the null models
Formulas4[length(Formulas4) + 1] <- "N~offset(log(striparea))+1"

##Loop through all formulas for the ds and save the AIC
gam.AIC <-vector("numeric",length(Formulas4))

for (i in 1:length(Formulas4)){
  try({
    formula=as.formula(Formulas4[i])
    gam.test <-mgcv::gam(formula,family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML")
    gam.AIC[i]<-gam.test$aic
  }, silent = TRUE)		
  print(i)
}

gam.AIC.out <- data.frame(model = Formulas4, AIC = gam.AIC)
gam.AIC.out<-gam.AIC.out[which(gam.AIC.out$AIC > 0),] #!is.na(gam.AIC.out) && 
gam.AIC.out<-gam.AIC.out[order(gam.AIC.out$AIC),]
View(gam.AIC.out)
write.csv(gam.AIC.out,"C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gam.AIC.out_limited.csv")

gam.formula <- as.character(gam.AIC.out$model[gam.AIC.out$AIC == min(gam.AIC[gam.AIC > 0])])

gam.best <-mgcv::gam(as.formula(gam.formula),family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML")
summary(gam.best)
windows()
par(mfrow=c(2,2))
plot(gam.best)
# 
# gam.best2 <-mgcv::gam(N~offset(log(striparea))+ti(Front_Z, k = 4, bs = "ts") + ti(SST, k = 4, bs = "ts"),
#                       family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML")
# summary(gam.best2)
# gam.best2$aic
# plot(gam.best2)

#### 3.6 GOMMAPPS - Create Predictive Grid ####
area_msq<-rgeos::gArea(domain,byid=TRUE)
domain$striparea<-area_msq #area of cell in square meters
temp<-domain
temp@data$SST<-NA
temp@data$Front_Z<-NA
temp@data$pp<-NA
temp@data$Depth_m<-NA
temp@data$Front<-NA
temp@data$ChlA<-NA
temp@data$waterv<-NA #HYCOM surface water velocity (north)
temp@data$oscarv<-NA #OSCAR 5-d composite surface water velocity (north ~ meridional)
temp@data$k490<-NA #VIIRSN 8-d Diffuse attenuation coefficient at 490 nm, KD2 algorithm, m^-1
temp@data$wave<-NA #WaveWatch III (WW3) Global Wave Model significant wave height
temp@data$Nhat<-NA
temp@data$Nhat_SE<-NA


#assign monthly average data to grid
setwd("D:/Satellite/SEFSC")

#extent
min_Lat<-min(temp@data$Lat_midpt)-1
max_Lat<-max(temp@data$Lat_midpt)+1
min_Long<-min(temp@data$Long_midpt)-1
max_Long<-max(temp@data$Long_midpt)+1

#get depth
domain_bathy<-getNOAA.bathy(min_Long,max_Long,min_Lat,max_Lat,resolution=4)

env_grid_mo<-function(month) {
  temp@data$month<-month
  d1<-as.Date(paste("2017",month,"16",sep="-"))
  
  #time-specific environmental covariates
  day_start<-paste(as.character(d1))
  
  destFile <- paste("./MURsst_m_",day_start,".nc",sep="")
  if((!file.exists(destFile)) | (file.info(destFile)$size<1000)) {
    myURL <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(',min_Lat,'):1:(',max_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
    junk <- GET(myURL, write_disk(destFile,overwrite=TRUE), progress())
  }
  current_sst<-raster(destFile)
  current_front<-brick(detectFronts(current_sst[[1]],method = "BelkinOReilly2009",intermediate = FALSE))
  max_front<-cellStats(current_front,stat="max")
  
  destFile2 <- paste("./MH1chla_m_",day_start,".nc",sep="")
  if((!file.exists(destFile2)) | (file.info(destFile2)$size<1000)) {
    myURL2 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday.nc?chlorophyll[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
    junk2 <- GET(myURL2, write_disk(destFile2,overwrite=TRUE), progress())
  }
  current_chla<-brick(destFile2)
  
  if(d1 > as.Date("2016-04-17") & d1 < as.Date("2018-11-21")) {
    destFile3 <- paste("./HYCOM_v_d_",day_start,".nc",sep="")
    if((!file.exists(destFile3)) | (file.info(destFile3)$size<1000)) {
      myURL3 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/nrlHycomGLBu008e912D_LonPM180.nc?water_v[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(0.0):1:(0.0)][(',min_Lat,'):1:(',max_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
      junk3 <- GET(myURL3, write_disk(destFile3,overwrite=TRUE), progress())
    }
    current_waterv<-brick(destFile3)
  } 
  
  if(d1 < as.Date("2018-06-12")) {
    destFile4 <- paste("./OSCAR_v_d_",day_start,".nc",sep="")
    if((!file.exists(destFile4)) | (file.info(destFile4)$size<1000)) {
      myURL4 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplOscar.nc?v[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(15.0):1:(15.0)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long+360,'):1:(',max_Long+360,')]',sep="")
      junk4 <- GET(myURL4, write_disk(destFile4,overwrite=TRUE), progress())
    }
    current_oscarv<-brick(destFile4)
  } 
  
  destFile5 <- paste("./MODIS_PP_m_",day_start,".nc",sep="")
  if((!file.exists(destFile5)) | (file.info(destFile5)$size<1000)) {
    myURL5 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1ppmday.nc?productivity[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(0.0):1:(0.0)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
    junk5 <- GET(myURL5, write_disk(destFile5,overwrite=TRUE), progress())
  }
  current_pp<-brick(destFile5)
  
  if(d1>as.Date("2012-01-05")) {
    destFile6 <- paste("./VIIRSN_k490_m_",day_start,".nc",sep="")
    if((!file.exists(destFile6)) | (file.info(destFile6)$size<1000)) {
      myURL6 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018k490mday.nc?k490[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
      junk6 <- GET(myURL6, write_disk(destFile6,overwrite=TRUE), progress())
    }
    current_k490<-brick(destFile6)
  } 
  
  if(d1 >= as.Date("2013-01-01")) {
    destFile7 <- paste("./WAVE_d_",day_start,".nc",sep="")
    if((!file.exists(destFile7)) | (file.info(destFile7)$size<1000)) {
      myURL7 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NWW3_Global_Best.nc?Thgt[(',day_start,'T09:00:00Z):1:(',day_start,'T09:00:00Z)][(0.0):1:(0.0)][(',min_Lat,'):1:(',max_Lat,')][(',min_Long+360,'):1:(',max_Long+360,')]',sep="")
      junk7 <- GET(myURL7, write_disk(destFile7,overwrite=TRUE), progress())
    }
    current_wave<-brick(destFile7)
  }
  
  #assign data to points
  for(i in 1:nrow(temp@data)) {
    #time-specific environmental covariates
    
    pt<-SpatialPoints(cbind((temp@data[i,"Long_midpt"]),temp@data[i,"Lat_midpt"])) #identify the current point (note this source has negative long option, no need to add 360)
    
    #depth
    if(is.na(temp@data[i,"Depth_m"])) {
      temp@data[i,"Depth_m"]<-get.depth(domain_bathy,x=temp@data[i,"Long_midpt"],y=temp@data[i,"Lat_midpt"],locator=FALSE)$depth
    }
    temp@data[i,"SST"]<-raster::extract(current_sst,pt)
    temp@data[i,"Front"]<-raster::extract(current_front,pt) #match the frontal gradient at the point to the correct day within the brick
    temp@data[i,"Front_Z"]<-temp@data[i,"Front"]/max_front #match the frontal gradient at the point to the correct day within the brick
    temp@data[i,"ChlA"]<-raster::extract(current_chla,pt,method='bilinear') #to fill gaps
    
    if(d1 > as.Date("2016-04-17") & d1 < as.Date("2018-11-21")) {
      temp@data[i,"waterv"]<-raster::extract(current_waterv,pt)
    } else {temp@data[i,"waterv"]<-NA}
    
    if(d1 < as.Date("2018-06-12")) {
      temp@data[i,"oscarv"]<-raster::extract(current_oscarv,SpatialPoints(cbind((temp@data[i,"Long_midpt"]+360),temp@data[i,"Lat_midpt"])),method="bilinear") #to fill gaps
    } else {temp@data[i,"oscarv"]<-NA}
    
    temp@data[i,"pp"]<-raster::extract(current_pp,pt,method="bilinear") #to fill gaps
    
    if(d1>as.Date("2012-01-05")) {
      temp@data[i,"k490"]<-raster::extract(current_k490,pt,method="bilinear") #to fill gaps
    } else {temp@data[i,"k490"]<-NA}
    
    if(d1 < as.Date("2013-01-01")) {
      temp@data[i,"wave"]<-NA
    } else {
      temp@data[i,"wave"]<-raster::extract(current_wave,SpatialPoints(cbind((temp@data[i,"Long_midpt"]+360),temp@data[i,"Lat_midpt"])))
    }
    
    #incremental write due to IT backup restarts
    if(i %% 5000 == 0) {
      print(paste("month=",month,", i=",i,", ",i/nrow(temp@data)*100,"% complete"))
      write.csv(temp@data,paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort/pred_grid_",month,".csv",sep=""))
    }
  }
  
  write.csv(temp@data,paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/SEFSC_Projected_Effort/pred_grid_",month,".csv",sep=""))
  print(paste("month= ",month," complete!"))
  return(temp@data)
  
}#end function env_grid_mo

library(parallel)
cl <- makeCluster((detectCores()-1))
clusterEvalQ(cl, { 
  library(sjlabelled)
  library(devtools)
  library(foreign)
  library(ggplot2)
  library(rgdal)        # for readOGR(...)
  library(data.table)
  library(plyr)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(maps)
  library(truncnorm)
  library(mrds)
  library(Distance)
  library(AICcmodavg)
  library(rgeos)
  library(maptools)
  library(raster)
  library(dplyr)
  library(mapproj)
  library(sp)
  library(ggmap)
  library(sf)
  library(rgdal)
  library(ncdf4)
  library(maps)
  library(matlab)
  library(grDevices)
  library(rerddap)
  library(lubridate)
  library(tibble)
  library(devtools)
  library(arm)
  library(lme4)
  library(maps)
  library(MASS)
  library(Hmisc) 
  library(plyr)  
  library(AICcmodavg)
  library(foreign)
  library(lubridate)
  library(R2admb)
  library(usdm)
  library(splines)
  library(coda)
  library(ncf)
  library(ade4)
  library(ecodist)
  library(boot)
  library(reshape)
  library(marmap)
  library(car)
  library(pscl)
  library(mgcv)
  library(colorRamps)
  library(adehabitatHR)
  library(imager)
  library(spatial.tools)
  library(classInt)
  library(rerddap)
  library(ncdf4)
  library(httr)
  library(grec)
  library(tidyselect)
  library(sjPlot)
  library(gamclass)
  library(magick)
  library(rerddap)
  library(ncdf4)
  library(stringi)
  library(mapdata)
  library(ggplot2)
  library(raster)
  library(httr)
  library(lubridate)
  library(rasterVis)
  library(RColorBrewer)
  library(maps)
  library(animation)
  library(dplyr)
  library(purrr)
  library(colorRamps)
  library(ncdf4)
  library(raster)
  library(wvtool)
  library(oceanmap)
  library(imager)
  library(geostatsp)
  library(maptools)
  library(spatstat)
  library(grec)
  library(rlang)
  library(dsm)
  library(ROCR)
  library(gam)
  library(Distance)
  library(captioner)
  library(tweedie)
  library(RODBC)
  library(foreign)
  library(knitr)
  library(mgcv)
})
clusterExport(cl,"temp")
clusterEvalQ(cl,temp)
clusterExport(cl,"domain_bathy")
clusterEvalQ(cl,domain_bathy)
clusterExport(cl,"min_Lat")
clusterEvalQ(cl,min_Lat)
clusterExport(cl,"max_Lat")
clusterEvalQ(cl,max_Lat)
clusterExport(cl,"min_Long")
clusterEvalQ(cl,min_Long)
clusterExport(cl,"max_Long")
clusterEvalQ(cl,max_Long)

months<-seq(1,12)
pred_grid<-parLapply(cl,months,env_grid_mo)
stopCluster(cl)
head(pred_grid[[1]])
head(pred_grid[[12]])

#### 3.7 GOMMAPPS - Fit to Predictive Grid ####
#gam.best
#get se.fit
gommapps.predict_2017_se<-vector("list",12)
#gommapps.predict_2017_se2<-vector("list",12)
for(i in 1:12) {
  gommapps.predict_2017_se[[i]] <- predict(gam.best, pred_grid[[i]],type="response", se.fit=T) #check type="response" with LPG; default (link) resulted in negative predictions
  #gommapps.predict_2017_se2[[i]] <- predict(gam.best2, pred_grid[[i]],type="response", se.fit=T) #check type="response" with LPG; default (link) resulted in negative predictions
}
summary(gommapps.predict_2017_se[[1]]$fit)
summary(gommapps.predict_2017_se[[7]]$fit)
hist(gommapps.predict_2017_se[[1]]$fit)

rm(pred_grid_po,pred_grid_pt,pred_grid_df)
pred_grid_po<-vector("list",12)

for(i in 1:12) {
  pred_grid_po[[i]]<-domain
  pred_grid_po[[i]]@data$SST<-pred_grid[[i]]$SST
  pred_grid_po[[i]]@data$Front_Z<-pred_grid[[i]]$Front_Z
  pred_grid_po[[i]]@data$pp<-pred_grid[[i]]$pp
  pred_grid_po[[i]]@data$Depth_m<-pred_grid[[i]]$Depth_m
  pred_grid_po[[i]]@data$AvailabilityBias<-ifelse(pred_grid[[i]]$Depth_m>=-2,1,ifelse(pred_grid[[i]]$Depth_m>=-15,0.25,0.07)) #3 tags on east coast=24.7% [15% (Gillie), 25% (Aidon), 34% (Oceane)]
  pred_grid_po[[i]]@data$Front<-pred_grid[[i]]$Front
  pred_grid_po[[i]]@data$ChlA<-pred_grid[[i]]$ChlA
  pred_grid_po[[i]]@data$waterv<-pred_grid[[i]]$waterv #HYCOM surface water velocity (north)
  pred_grid_po[[i]]@data$oscarv<-pred_grid[[i]]$oscarv #OSCAR 5-d composite surface water velocity (north ~ meridional)
  pred_grid_po[[i]]@data$k490<-pred_grid[[i]]$k490 #VIIRSN 8-d Diffuse attenuation coefficient at 490 nm, KD2 algorithm, m^-1
  pred_grid_po[[i]]@data$wave<-pred_grid[[i]]$wave #WaveWatch III (WW3) Global Wave Model significant wave height
  pred_grid_po[[i]]@data$Nhat<-gommapps.predict_2017_se[[i]]$fit #N/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_SE<-gommapps.predict_2017_se[[i]]$se.fit #SE/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_exp<-(gommapps.predict_2017_se[[i]]$fit/pred_grid_po[[i]]@data$AvailabilityBias) #N perception & availability bias  #O:\Species Conservation\Giant Manta Ray\Research\FGBNMS Manta Tagging Aug 2019\Aug 7\Luana Tag 03\Luana Tag Data\181630-Series_20190912.xlsx
  pred_grid_po[[i]]@data$Nhat_SE_exp<-(gommapps.predict_2017_se[[i]]$se.fit/pred_grid_po[[i]]@data$AvailabilityBias) #SE perception & availability bias 
  # pred_grid_po[[i]]@data$Nhat2<-gommapps.predict_2017_se2[[i]]$fit #N/meters squared perception bias
  # pred_grid_po[[i]]@data$Nhat2_SE<-gommapps.predict_2017_se2[[i]]$se.fit #SE/meters squared perception bias
  # pred_grid_po[[i]]@data$Nhat2_exp<-(gommapps.predict_2017_se2[[i]]$fit/pred_grid_po[[i]]@data$AvailabilityBias) #N perception & availability bias  #O:\Species Conservation\Giant Manta Ray\Research\FGBNMS Manta Tagging Aug 2019\Aug 7\Luana Tag 03\Luana Tag Data\181630-Series_20190912.xlsx
  # pred_grid_po[[i]]@data$Nhat2_SE_exp<-(gommapps.predict_2017_se2[[i]]$se.fit/pred_grid_po[[i]]@data$AvailabilityBias) #SE perception & availability bias 
  pred_grid_po[[i]]$id<-rownames(as.data.frame(pred_grid_po[[i]]))
}

pred_grid_sf<-vector("list",12)
for(i in 1:12) {
  pred_grid_sf[[i]]<-st_as_sf(pred_grid_po[[i]])
}

quantile(pred_grid_sf[[12]]$Nhat_exp,na.rm=T)

#define max values for plotting to reduce effect of outliers on plot
for(i in 1:12) {
  pred_grid_sf[[i]]$Nhat_exp_plot<-ifelse(pred_grid_sf[[i]]$Nhat_exp>100,100,pred_grid_sf[[i]]$Nhat_exp)
  #pred_grid_sf[[i]]$Nhat2_exp_plot<-ifelse(pred_grid_sf[[i]]$Nhat2_exp>100,100,pred_grid_sf[[i]]$Nhat2_exp)
}

p<-vector("list",12)
#p2<-vector("list",12)
for(i in 1:12) {
  grob <- grobTree(textGrob(i, x=0.1,  y=0.8, hjust=0,
                            gp=gpar(col="white", fontsize=13, fontface="italic")))
  p[[i]]<-ggplot(pred_grid_sf[[i]]) + 
    geom_sf(aes(fill=Nhat_exp_plot),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))+
    annotation_custom(grob)+
    theme_void()+
    theme(legend.position = "none")
  # p2[[i]]<-ggplot(pred_grid_sf[[i]]) + 
  #   geom_sf(aes(fill=Nhat2_exp_plot),colour = NA) + 
  #   scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))+
  #   annotation_custom(grob)+
  #   theme_void()+
  #   theme(legend.position = "none")
}
p_temp<-ggplot(pred_grid_sf[[1]]) + 
  geom_sf(aes(fill=Nhat_exp),colour = NA) + 
  scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))

library(ggpubr)
leg<-get_legend(p_temp)
pleg<-as_ggplot(leg)

# png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest.png",sep=""),width=9,height=6.5,res=600,units="in")
# p1grid<-grid.arrange(arrangeGrob(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
#                 p[[10]],p[[11]],p[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
# p1grid
# dev.off()
# png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest2.png",sep=""),width=9,height=6.5,res=600,units="in")
# p2grid<-grid.arrange(arrangeGrob(p2[[1]],p2[[2]],p2[[3]],p2[[4]],p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],
#                             p2[[10]],p2[[11]],p2[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
# p2grid
# dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_formula.png",sep=""),width=9,height=6.5,res=600,units="in")
p1grid<-grid.arrange(arrangeGrob(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
                                 p[[10]],p[[11]],p[[12]],bottom=formula(gam.best),nrow=4),pleg,nrow=1,ncol=2,widths=c(9,1))
p1grid
grid.text(formula(gam.best), x = unit(0.9, "npc"),
          y = unit(0.5, "npc"),gp = gpar(fontsize=10, fontfamily="Times New Roman"),rot=90)
dev.off()

# png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest2_formula.png",sep=""),width=8.5,height=6.5,res=600,units="in")
# p2grid<-grid.arrange(arrangeGrob(p2[[1]],p2[[2]],p2[[3]],p2[[4]],p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],
#                                  p2[[10]],p2[[11]],p2[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(9,1))
# p2grid
# grid.text(formula(gam.best2), x = unit(0.9, "npc"),
#           y = unit(0.5, "npc"),gp = gpar(fontsize=10, fontfamily="Times New Roman"),rot=90)
# dev.off()

#### 3.8 GOMMAPPS - Fit GAM to full suite of parameters for comparison ####
DF <- data.frame(Class=1:n.var, Front=1:n.var, Front_Z=1:n.var, SST=1:n.var, ChlA=1:n.var, waterv=1:n.var,
                 pp=1:n.var, k490=1:n.var, wave=1:n.var, #eliminated oscar-v
                 Depth_m=1:n.var)
# DF <- data.frame(Class=1:n.var, Front_Z=1:n.var, SST=1:n.var, k490=1:n.var,
#                  Depth_m=1:n.var)
Cols <- names(DF)
Cols <- Cols[! Cols %in% "Class"]
n <- length(Cols)
id<-unlist(lapply(1:n,function(i)combn(1:n,i,simplify=F)),recursive=F)
Formulas<-sapply(id,function(i)paste("N~offset(log(striparea))+",paste('ti(',Cols[i],',k=3,bs="ts")',collapse="+",sep="")))
Formulas2<-sapply(id,function(i)paste("N~offset(log(striparea))+ti(Front_Z,SST,k=c(3,3),bs=c('ts','ts'))+",paste('ti(',Cols[i],',k=3,bs="ts")',collapse="+",sep="")))
Formulas3<-sapply(id,function(i)paste("N~offset(log(striparea))+ti(Front_Z,SST,k=c(3,3),bs=c('ts','ts'))+ti(Front_Z,SST,k490,k=c(3,3,3),bs=c('ts','ts','ts'))+",paste('ti(',Cols[i],',k=3,bs="ts")',collapse="+",sep="")))
Formulas4<-c(Formulas,Formulas2,Formulas3)
#Add the null models
Formulas4[length(Formulas4) + 1] <- "N~offset(log(striparea))+1"

##Loop through all formulas for the ds and save the AIC
gam.AIC.full <-vector("numeric",length(Formulas4))

for (i in 1:length(Formulas4)){
  try({
    formula=as.formula(Formulas4[i])
    gam.test <-mgcv::gam(formula,family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML")
    gam.AIC.full[i]<-gam.test$aic
  }, silent = TRUE)		
  print(i)
}

gam.AIC.out.full <- data.frame(model = Formulas4, AIC = gam.AIC.full)
gam.AIC.out.full<-gam.AIC.out.full[which(gam.AIC.out.full$AIC > 0),] #!is.na(gam.AIC.out) && 
gam.AIC.out.full<-gam.AIC.out.full[order(gam.AIC.out.full$AIC),]
View(gam.AIC.out.full)
write.csv(gam.AIC.out.full,"C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gam_AIC_out_full.csv")

gam.formula.full <- as.character(gam.AIC.out.full$model[gam.AIC.out.full$AIC == min(gam.AIC.full[gam.AIC.full > 0])])

gam.best.full <-mgcv::gam(as.formula(gam.formula.full),family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML")
summary(gam.best.full)
windows()
par(mfrow=c(3,3))
plot(gam.best.full)

#manually eliminate insigificant terms and highly correlated terms
#N ~ offset(log(striparea)) + ti(Front_Z, SST, k = c(3, 3), bs = c("ts","ts")) + ti(Front_Z, k = 3, bs = "ts") + ti(SST,k = 3, bs = "ts") + ti(waterv, k = 3, bs = "ts") + 
#ti(pp, k = 3, bs = "ts") + ti(k490, k = 3, bs = "ts") + ti(wave, k = 3, bs = "ts") + ti(Depth_m, k = 3, bs = "ts")

gam.formula.full2<-'N~offset(log(striparea))+ ti(Front_Z,k=3,bs="ts")+ti(SST,k=3,bs="ts")+ti(waterv,k=3,bs="ts")+ti(k490,k=3,bs="ts")'
gam.best.full2 <-mgcv::gam(as.formula(gam.formula.full2),family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML",select=T)
summary(gam.best.full2)
windows()
par(mfrow=c(2,2))
plot(gam.best.full2)

gam.formula.full3<-'N~offset(log(striparea))+ ti(Front_Z,k=3,bs="ts")+ti(SST,k=3,bs="ts")+ti(waterv,k=3,bs="ts")+ti(pp,k=3,bs="ts")'
gam.best.full3 <-mgcv::gam(as.formula(gam.formula.full3),family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML",select=T)
summary(gam.best.full3)
windows()
par(mfrow=c(2,2))
plot(gam.best.full3)

gam.formula.full4<-'N~offset(log(striparea))+ ti(Front,k=3,bs="ts")+ti(SST,k=3,bs="ts")+ti(waterv,k=3,bs="ts")+ti(k490,k=3,bs="ts")'
gam.best.full4 <-mgcv::gam(as.formula(gam.formula.full4),family=tw(),data=gommapps_grid,na.action=na.exclude,method="REML",select=T)
summary(gam.best.full4)
windows()
par(mfrow=c(2,2))
plot(gam.best.full4)


#### 3.9 GOMMAPPS - Fit GAM.full to Predictive Grid ####
gommapps.predict_2017_se.full<-vector("list",12)
gommapps.predict_2017_se.full2<-vector("list",12)
gommapps.predict_2017_se.full3<-vector("list",12)
gommapps.predict_2017_se.full4<-vector("list",12)
for(i in 1:12) {
  gommapps.predict_2017_se.full[[i]] <- predict(gam.best.full, pred_grid[[i]],type="response", se.fit=T) #check type="response" with LPG; default (link) resulted in negative predictions
  gommapps.predict_2017_se.full2[[i]] <- predict(gam.best.full2, pred_grid[[i]],type="response", se.fit=T) #check type="response" with LPG; default (link) resulted in negative predictions
  gommapps.predict_2017_se.full3[[i]] <- predict(gam.best.full3, pred_grid[[i]],type="response", se.fit=T) #check type="response" with LPG; default (link) resulted in negative predictions
  gommapps.predict_2017_se.full4[[i]] <- predict(gam.best.full4, pred_grid[[i]],type="response", se.fit=T) #check type="response" with LPG; default (link) resulted in negative predictions
}

for(i in 1:12) {
  pred_grid_po[[i]]@data$Nhat_full<-gommapps.predict_2017_se.full[[i]]$fit #N/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_SE_full<-gommapps.predict_2017_se.full[[i]]$se.fit #SE/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_exp_full<-(gommapps.predict_2017_se.full[[i]]$fit/pred_grid_po[[i]]@data$AvailabilityBias) #N perception & availability bias  #O:\Species Conservation\Giant Manta Ray\Research\FGBNMS Manta Tagging Aug 2019\Aug 7\Luana Tag 03\Luana Tag Data\181630-Series_20190912.xlsx
  pred_grid_po[[i]]@data$Nhat_SE_exp_full<-(gommapps.predict_2017_se.full[[i]]$se.fit/pred_grid_po[[i]]@data$AvailabilityBias) #SE perception & availability bias 
  pred_grid_po[[i]]@data$Nhat_full2<-gommapps.predict_2017_se.full2[[i]]$fit #N/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_SE_full2<-gommapps.predict_2017_se.full2[[i]]$se.fit #SE/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_exp_full2<-(gommapps.predict_2017_se.full2[[i]]$fit/pred_grid_po[[i]]@data$AvailabilityBias) #N perception & availability bias  #O:\Species Conservation\Giant Manta Ray\Research\FGBNMS Manta Tagging Aug 2019\Aug 7\Luana Tag 03\Luana Tag Data\181630-Series_20190912.xlsx
  pred_grid_po[[i]]@data$Nhat_SE_exp_full2<-(gommapps.predict_2017_se.full2[[i]]$se.fit/pred_grid_po[[i]]@data$AvailabilityBias) #SE perception & availability bias 
  pred_grid_po[[i]]@data$Nhat_full3<-gommapps.predict_2017_se.full3[[i]]$fit #N/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_SE_full3<-gommapps.predict_2017_se.full3[[i]]$se.fit #SE/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_exp_full3<-(gommapps.predict_2017_se.full3[[i]]$fit/pred_grid_po[[i]]@data$AvailabilityBias) #N perception & availability bias  #O:\Species Conservation\Giant Manta Ray\Research\FGBNMS Manta Tagging Aug 2019\Aug 7\Luana Tag 03\Luana Tag Data\181630-Series_20190912.xlsx
  pred_grid_po[[i]]@data$Nhat_SE_exp_full3<-(gommapps.predict_2017_se.full3[[i]]$se.fit/pred_grid_po[[i]]@data$AvailabilityBias) #SE perception & availability bias 
  pred_grid_po[[i]]@data$Nhat_full4<-gommapps.predict_2017_se.full4[[i]]$fit #N/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_SE_full4<-gommapps.predict_2017_se.full4[[i]]$se.fit #SE/meters squared perception bias
  pred_grid_po[[i]]@data$Nhat_exp_full4<-(gommapps.predict_2017_se.full4[[i]]$fit/pred_grid_po[[i]]@data$AvailabilityBias) #N perception & availability bias  #O:\Species Conservation\Giant Manta Ray\Research\FGBNMS Manta Tagging Aug 2019\Aug 7\Luana Tag 03\Luana Tag Data\181630-Series_20190912.xlsx
  pred_grid_po[[i]]@data$Nhat_SE_exp_full4<-(gommapps.predict_2017_se.full4[[i]]$se.fit/pred_grid_po[[i]]@data$AvailabilityBias) #SE perception & availability bias 
}

#define max values for plotting to reduce effect of outliers on plot
for(i in 1:12) {
  pred_grid_sf[[i]]$Nhat_full<-pred_grid_po[[i]]@data$Nhat_full
  pred_grid_sf[[i]]$Nhat_SE_full<-pred_grid_po[[i]]@data$Nhat_SE_full
  pred_grid_sf[[i]]$Nhat_exp_full<-pred_grid_po[[i]]@data$Nhat_exp_full
  pred_grid_sf[[i]]$Nhat_SE_exp_full<-pred_grid_po[[i]]@data$Nhat_SE_exp_full
  pred_grid_sf[[i]]$Nhat_exp_plot_full<-ifelse(pred_grid_sf[[i]]$Nhat_exp_full>100,100,pred_grid_sf[[i]]$Nhat_exp_full)
  pred_grid_sf[[i]]$Nhat_full2<-pred_grid_po[[i]]@data$Nhat_full2
  pred_grid_sf[[i]]$Nhat_SE_full2<-pred_grid_po[[i]]@data$Nhat_SE_full2
  pred_grid_sf[[i]]$Nhat_exp_full2<-pred_grid_po[[i]]@data$Nhat_exp_full2
  pred_grid_sf[[i]]$Nhat_SE_exp_full2<-pred_grid_po[[i]]@data$Nhat_SE_exp_full2
  pred_grid_sf[[i]]$Nhat_exp_plot_full2<-ifelse(pred_grid_sf[[i]]$Nhat_exp_full2>100,100,pred_grid_sf[[i]]$Nhat_exp_full2)
  pred_grid_sf[[i]]$Nhat_full3<-pred_grid_po[[i]]@data$Nhat_full3
  pred_grid_sf[[i]]$Nhat_SE_full3<-pred_grid_po[[i]]@data$Nhat_SE_full3
  pred_grid_sf[[i]]$Nhat_exp_full3<-pred_grid_po[[i]]@data$Nhat_exp_full3
  pred_grid_sf[[i]]$Nhat_SE_exp_full3<-pred_grid_po[[i]]@data$Nhat_SE_exp_full3
  pred_grid_sf[[i]]$Nhat_exp_plot_full3<-ifelse(pred_grid_sf[[i]]$Nhat_exp_full3>100,100,pred_grid_sf[[i]]$Nhat_exp_full3)
  pred_grid_sf[[i]]$Nhat_full4<-pred_grid_po[[i]]@data$Nhat_full4
  pred_grid_sf[[i]]$Nhat_SE_full4<-pred_grid_po[[i]]@data$Nhat_SE_full4
  pred_grid_sf[[i]]$Nhat_exp_full4<-pred_grid_po[[i]]@data$Nhat_exp_full4
  pred_grid_sf[[i]]$Nhat_SE_exp_full4<-pred_grid_po[[i]]@data$Nhat_SE_exp_full4
  pred_grid_sf[[i]]$Nhat_exp_plot_full4<-ifelse(pred_grid_sf[[i]]$Nhat_exp_full4>100,100,pred_grid_sf[[i]]$Nhat_exp_full4)
}

fgbnms<-st_read("C:/Users/nick.farmer/Documents/GIS/Base/Gulf/Management Areas/FGBNMS Boundary/fgbnms_pl_ed.shp")
all_manta_pts<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/Anecdotal Observations/ALL Manta_pts_all_20200323.csv",stringsAsFactors = F)
all_manta_pts$Date<-as.Date(all_manta_pts$Date,format="%m/%d/%Y")
all_manta_pts$Lat<-as.numeric(all_manta_pts$Lat)
all_manta_pts$Long<-as.numeric(all_manta_pts$Long)
all_manta_pts$Month<-as.numeric(all_manta_pts$Month)
all_manta_pts$Year<-as.numeric(all_manta_pts$Year)
summary(all_manta_pts)
all_manta_pts.sf<-st_as_sf(all_manta_pts[which(!is.na(all_manta_pts$Lat)),],coords=c("Long","Lat"),crs=wgs.84)

# windows()
# pt_mo<-all_manta_pts.sf[which(all_manta_pts.sf$Month==6 & all_manta_pts.sf$Year==2017),]
# ggplot(pred_grid_sf[[6]]) + labs(title=formula(gam.best.full),subtitle=paste("AIC = ",gam.best.full$aic)) +
#   geom_sf(aes(fill=Nhat_exp_full),colour = NA) + 
#   scale_fill_distiller(palette ="Spectral", name="N")+
#   #annotation_custom(grob)+
#   #theme_void()+
#   theme(legend.position = c(0.1,0.75))+
#   geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
#   geom_sf(data=pt_mo,shape=3,size=2,color="black")+
#   coord_sf(xlim = c(min(pred_grid_sf[[6]]$x), max(pred_grid_sf[[6]]$x)), ylim = c(min(pred_grid_sf[[6]]$y), max(pred_grid_sf[[6]]$y)), expand = FALSE)

p.full<-vector("list",12)
p.full2<-vector("list",12)
p.full3<-vector("list",12)
p.full4<-vector("list",12)
for(i in 1:12) {
  grob <- grobTree(textGrob(i, x=0.1,  y=0.8, hjust=0,
                            gp=gpar(col="white", fontsize=13, fontface="italic")))
  p.full[[i]]<-ggplot(pred_grid_sf[[i]]) + 
    geom_sf(aes(fill=Nhat_exp_plot_full),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))+
    annotation_custom(grob)+
    theme_void()+
    theme(legend.position = "none")+
    geom_sf(data=fgbnms,color="black",fill=NA,size=1)
  p.full2[[i]]<-ggplot(pred_grid_sf[[i]]) + 
    geom_sf(aes(fill=Nhat_exp_plot_full2),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))+
    annotation_custom(grob)+
    theme_void()+
    theme(legend.position = "none")+
    geom_sf(data=fgbnms,color="black",fill=NA,size=1)
  p.full3[[i]]<-ggplot(pred_grid_sf[[i]]) + 
    geom_sf(aes(fill=Nhat_exp_plot_full3),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))+
    annotation_custom(grob)+
    theme_void()+
    theme(legend.position = "none")+
    geom_sf(data=fgbnms,color="black",fill=NA,size=1)
  p.full4[[i]]<-ggplot(pred_grid_sf[[i]]) + 
    geom_sf(aes(fill=Nhat_exp_plot_full4),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))+
    annotation_custom(grob)+
    theme_void()+
    theme(legend.position = "none")+
    geom_sf(data=fgbnms,color="black",fill=NA,size=1)
}
p_temp<-ggplot(pred_grid_sf[[1]]) + 
  geom_sf(aes(fill=Nhat_exp_plot_full),colour = NA) + 
  scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))

library(ggpubr)
leg<-get_legend(p_temp)
pleg<-as_ggplot(leg)

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_full.png",sep=""),width=9,height=6.5,res=600,units="in")
p1grid.full<-grid.arrange(arrangeGrob(p.full[[1]],p.full[[2]],p.full[[3]],p.full[[4]],p.full[[5]],p.full[[6]],p.full[[7]],p.full[[8]],p.full[[9]],
                                      p.full[[10]],p.full[[11]],p.full[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
p1grid.full
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_formula_full.png",sep=""),width=9,height=6.5,res=600,units="in")
p1grid.full<-grid.arrange(arrangeGrob(p.full[[1]],p.full[[2]],p.full[[3]],p.full[[4]],p.full[[5]],p.full[[6]],p.full[[7]],p.full[[8]],p.full[[9]],
                                      p.full[[10]],p.full[[11]],p.full[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
p1grid.full
grid.text(formula(gam.best.full), x = unit(0.9, "npc"),
          y = unit(0.5, "npc"),gp = gpar(fontsize=6, fontfamily="Times New Roman"),rot=90)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_full2.png",sep=""),width=9,height=6.5,res=600,units="in")
pgrid.full2<-grid.arrange(arrangeGrob(p.full2[[1]],p.full2[[2]],p.full2[[3]],p.full2[[4]],p.full2[[5]],p.full2[[6]],p.full2[[7]],p.full2[[8]],p.full2[[9]],
                                      p.full2[[10]],p.full2[[11]],p.full2[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
pgrid.full2
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_formula_full2.png",sep=""),width=9,height=6.5,res=600,units="in")
pgrid.full2<-grid.arrange(arrangeGrob(p.full2[[1]],p.full2[[2]],p.full2[[3]],p.full2[[4]],p.full2[[5]],p.full2[[6]],p.full2[[7]],p.full2[[8]],p.full2[[9]],
                                      p.full2[[10]],p.full2[[11]],p.full2[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
pgrid.full2
grid.text(formula(gam.best.full2), x = unit(0.9, "npc"),
          y = unit(0.5, "npc"),gp = gpar(fontsize=6, fontfamily="Times New Roman"),rot=90)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_full3.png",sep=""),width=9,height=6.5,res=600,units="in")
pgrid.full3<-grid.arrange(arrangeGrob(p.full3[[1]],p.full3[[2]],p.full3[[3]],p.full3[[4]],p.full3[[5]],p.full3[[6]],p.full3[[7]],p.full3[[8]],p.full3[[9]],
                                      p.full3[[10]],p.full3[[11]],p.full3[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
pgrid.full3
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_formula_full3.png",sep=""),width=9,height=6.5,res=600,units="in")
pgrid.full3<-grid.arrange(arrangeGrob(p.full3[[1]],p.full3[[2]],p.full3[[3]],p.full3[[4]],p.full3[[5]],p.full3[[6]],p.full3[[7]],p.full3[[8]],p.full3[[9]],
                                      p.full3[[10]],p.full3[[11]],p.full3[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
pgrid.full3
grid.text(formula(gam.best.full3), x = unit(0.9, "npc"),
          y = unit(0.5, "npc"),gp = gpar(fontsize=6, fontfamily="Times New Roman"),rot=90)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_full3.png",sep=""),width=9,height=6.5,res=600,units="in")
pgrid.full4<-grid.arrange(arrangeGrob(p.full4[[1]],p.full4[[2]],p.full4[[3]],p.full4[[4]],p.full4[[5]],p.full4[[6]],p.full4[[7]],p.full4[[8]],p.full4[[9]],
                                      p.full4[[10]],p.full4[[11]],p.full4[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
pgrid.full4
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gambest_formula_full4.png",sep=""),width=9,height=6.5,res=600,units="in")
pgrid.full4<-grid.arrange(arrangeGrob(p.full4[[1]],p.full4[[2]],p.full4[[3]],p.full4[[4]],p.full4[[5]],p.full4[[6]],p.full4[[7]],p.full4[[8]],p.full4[[9]],
                                      p.full4[[10]],p.full4[[11]],p.full4[[12]],nrow=4,ncol=3),pleg,nrow=1,ncol=2,widths=c(10,1))
pgrid.full4
grid.text(formula(gam.best.full4), x = unit(0.9, "npc"),
          y = unit(0.5, "npc"),gp = gpar(fontsize=6, fontfamily="Times New Roman"),rot=90)
dev.off()

#plot vs. observations
p<-vector("list",12)
#p2<-vector("list",12)
p.full<-vector("list",12)
p.full2<-vector("list",12)
p.full3<-vector("list",12)
p.full4<-vector("list",12)
p_temp<-ggplot(pred_grid_sf[[1]]) + 
  geom_sf(aes(fill=Nhat_exp_plot_full),colour = NA) + 
  scale_fill_distiller(palette ="Spectral", name="N",limits=c(0,100))

library(ggpubr)
leg<-get_legend(p_temp)
pleg<-as_ggplot(leg)

for(i in 1:12) {
  grob <- grobTree(textGrob(i, x=0.1,  y=0.8, hjust=0,
                            gp=gpar(col="white", fontsize=13, fontface="italic")))
  pt_mo<-all_manta_pts.sf[which(all_manta_pts.sf$Month==i & all_manta_pts.sf$Year==2017),]
  p[[i]]<-ggplot(pred_grid_sf[[i]]) + labs(title=formula(gam.best),subtitle=paste("AIC = ",gam.best$aic)) +
    geom_sf(aes(fill=Nhat_exp_plot),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N")+
    annotation_custom(grob)+
    theme(legend.position = "none",plot.title = element_text(size=6),plot.subtitle = element_text(size=6))+
    geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
    geom_sf(data=pt_mo,shape=3,size=2,color="black")+
    coord_sf(xlim = c(min(pred_grid_sf[[i]]$x), max(pred_grid_sf[[i]]$x)), ylim = c(min(pred_grid_sf[[i]]$y), max(pred_grid_sf[[i]]$y)), expand = FALSE)
  # p2[[i]]<-ggplot(pred_grid_sf[[i]]) + labs(title=formula(gam.best2),subtitle=paste("AIC = ",gam.best2$aic)) +
  #   geom_sf(aes(fill=Nhat2_exp_plot),colour = NA) + 
  #   scale_fill_distiller(palette ="Spectral", name="N")+
  #   annotation_custom(grob)+
  #   theme(legend.position = "none",plot.title = element_text(size=6),plot.subtitle = element_text(size=6))+
  #   geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
  #   geom_sf(data=pt_mo,shape=3,size=2,color="black")+
  #   coord_sf(xlim = c(min(pred_grid_sf[[i]]$x), max(pred_grid_sf[[i]]$x)), ylim = c(min(pred_grid_sf[[i]]$y), max(pred_grid_sf[[i]]$y)), expand = FALSE)
  p.full[[i]]<-ggplot(pred_grid_sf[[i]]) + labs(title=formula(gam.best.full),subtitle=paste("AIC = ",gam.best.full$aic)) +
    geom_sf(aes(fill=Nhat_exp_plot_full),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N")+
    annotation_custom(grob)+
    theme(legend.position = "none",plot.title = element_text(size=6),plot.subtitle = element_text(size=6))+
    geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
    geom_sf(data=pt_mo,shape=3,size=2,color="black")+
    coord_sf(xlim = c(min(pred_grid_sf[[i]]$x), max(pred_grid_sf[[i]]$x)), ylim = c(min(pred_grid_sf[[i]]$y), max(pred_grid_sf[[i]]$y)), expand = FALSE)
  p.full2[[i]]<-ggplot(pred_grid_sf[[i]]) + labs(title=formula(gam.best.full2),subtitle=paste("AIC = ",gam.best.full2$aic)) +
    geom_sf(aes(fill=Nhat_exp_plot_full2),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N")+
    annotation_custom(grob)+
    theme(legend.position = "none",plot.title = element_text(size=6),plot.subtitle = element_text(size=6))+
    geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
    geom_sf(data=pt_mo,shape=3,size=2,color="black")+
    coord_sf(xlim = c(min(pred_grid_sf[[i]]$x), max(pred_grid_sf[[i]]$x)), ylim = c(min(pred_grid_sf[[i]]$y), max(pred_grid_sf[[i]]$y)), expand = FALSE)
  p.full3[[i]]<-ggplot(pred_grid_sf[[i]]) + labs(title=formula(gam.best.full3),subtitle=paste("AIC = ",gam.best.full3$aic)) +
    geom_sf(aes(fill=Nhat_exp_plot_full3),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N")+
    annotation_custom(grob)+
    theme(legend.position = "none",plot.title = element_text(size=6),plot.subtitle = element_text(size=6))+
    geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
    geom_sf(data=pt_mo,shape=3,size=2,color="black")+
    coord_sf(xlim = c(min(pred_grid_sf[[i]]$x), max(pred_grid_sf[[i]]$x)), ylim = c(min(pred_grid_sf[[i]]$y), max(pred_grid_sf[[i]]$y)), expand = FALSE)
  p.full4[[i]]<-ggplot(pred_grid_sf[[i]]) + labs(title=formula(gam.best.full3),subtitle=paste("AIC = ",gam.best.full3$aic)) +
    geom_sf(aes(fill=Nhat_exp_plot_full4),colour = NA) + 
    scale_fill_distiller(palette ="Spectral", name="N")+
    annotation_custom(grob)+
    theme(legend.position = "none",plot.title = element_text(size=6),plot.subtitle = element_text(size=6))+
    geom_sf(data=fgbnms,color="white",fill=NA,size=1)+
    geom_sf(data=pt_mo,shape=3,size=2,color="black")+
    coord_sf(xlim = c(min(pred_grid_sf[[i]]$x), max(pred_grid_sf[[i]]$x)), ylim = c(min(pred_grid_sf[[i]]$y), max(pred_grid_sf[[i]]$y)), expand = FALSE)
  
  png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_",i,".png",sep=""),width=9,height=6.5,res=600,units="in")
  #pgrid.mo<-grid.arrange(arrangeGrob(p[[i]],p2[[i]],p.full[[i]],p.full2[[i]],nrow=2,top=i),pleg,nrow=1,ncol=2,widths=c(10,1))
  pgrid.mo<-grid.arrange(arrangeGrob(p[[i]],p.full[[i]],p.full2[[i]],p.full3[[i]],nrow=2,top=i),pleg,nrow=1,ncol=2,widths=c(10,1))
  pgrid.mo
  dev.off()
  
  png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_full4_",i,".png",sep=""),width=9,height=6.5,res=600,units="in")
  #pgrid.mo<-grid.arrange(arrangeGrob(p[[i]],p2[[i]],p.full[[i]],p.full2[[i]],nrow=2,top=i),pleg,nrow=1,ncol=2,widths=c(10,1))
  pgrid.mo4<-grid.arrange(arrangeGrob(p.full4[[i]],nrow=1,top=i),pleg,nrow=1,ncol=2,widths=c(10,1))
  pgrid.mo4
  dev.off()
}

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_gam.best.full4.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
gam.check(gam.best.full4)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_gam.best.full3.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
gam.check(gam.best.full3)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_gam.best.full2.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
gam.check(gam.best.full2)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_gam.best.full.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
gam.check(gam.best.full)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_gam.best.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
gam.check(gam.best)
dev.off()

# png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_gamfits_gam.best2.png",sep=""),width=9,height=6.5,res=600,units="in")
# par(mfrow=c(2,2))
# gam.check(gam.best2)
# dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_visgam_gam.best_front.png",sep=""),width=9,height=6.5,res=600,units="in")
vis.gam(gam.best,view=c("SST","Front_Z"),type="response",color="heat",ticktype="detailed",theta=-35)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_visgam_gam.best.png",sep=""),width=9,height=6.5,res=600,units="in")
vis.gam(gam.best,view=c("SST","k490"),type="response",color="heat",ticktype="detailed",theta=-55)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_visgam_gam.best.full2.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
vis.gam(gam.best.full2,view=c("SST","k490"),type="response",color="heat",ticktype="detailed",theta=-35)
vis.gam(gam.best.full2,view=c("Front_Z","k490"),type="response",color="heat",ticktype="detailed",theta=-35)
vis.gam(gam.best.full2,view=c("waterv","k490"),type="response",color="heat",ticktype="detailed",theta=-35)
vis.gam(gam.best.full2,view=c("SST","waterv"),type="response",color="heat",ticktype="detailed",theta=-35)
dev.off()

png(paste("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/SEFSC Aerial Surveys/GOMMAPPS_visgam_gam.best.full3.png",sep=""),width=9,height=6.5,res=600,units="in")
par(mfrow=c(2,2))
vis.gam(gam.best.full3,view=c("SST","pp"),type="response",color="heat",ticktype="detailed",theta=-35)
vis.gam(gam.best.full3,view=c("Front_Z","pp"),type="response",color="heat",ticktype="detailed",theta=-35)
vis.gam(gam.best.full3,view=c("waterv","pp"),type="response",color="heat",ticktype="detailed",theta=-35)
vis.gam(gam.best.full3,view=c("SST","waterv"),type="response",color="heat",ticktype="detailed",theta=-35)
dev.off()

#### 3.10 GOMMAPPS - Assign satellite data to independent manta pts ####
#extent
min_Lat<-min(domain@data$Lat_midpt)
max_Lat<-max(domain@data$Lat_midpt)
min_Long<-min(domain@data$Long_midpt)
max_Long<-max(domain@data$Long_midpt)

#get depth
domain_bathy<-getNOAA.bathy(min_Long,max_Long,min_Lat,max_Lat,resolution=4)

setwd("C:/Users/nick.farmer/Desktop/temp")

start_i<-1

#assign environmental data to all independent observations
for(i in start_i:nrow(all_manta_pts)) {
  #time-specific environmental covariates
  d1<-all_manta_pts[i,"Date"]
  if(!(is.na(d1)) & !(is.na(all_manta_pts[i,"Lat"]))) {
    day_start<-paste(as.character(d1))
    
    print(paste("i=",i,", ",i/nrow(all_manta_pts)*100,"% complete"))
    
    pt<-SpatialPoints(cbind((all_manta_pts[i,"Long"]),all_manta_pts[i,"Lat"])) #identify the current point (note this source has negative long option, no need to add 360)
    
    if(((!is.na(all_manta_pts[(i-1),"Date"]) && (all_manta_pts[i,"Date"] != all_manta_pts[(i-1),"Date"])) | (is.na(all_manta_pts[(i-1),"Date"])))) { #don't need new data if it's not a new date
      
      if(d1>as.Date("2002-06-01")) {
        destFile <- paste("./MURsst_d_",day_start,".nc",sep="")
        if((!file.exists(destFile)) | (file.info(destFile)$size<1000)) {
          myURL <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.nc?analysed_sst[(',day_start,'T09:00:00Z):1:(',day_start,'T09:00:00Z)][(',min_Lat,'):1:(',max_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
          junk <- GET(myURL, write_disk(destFile,overwrite=TRUE), progress())
        }
        current_sst<-raster(destFile)
        current_front<-brick(detectFronts(current_sst[[1]],method = "BelkinOReilly2009",intermediate = FALSE))
      } else {
        current_sst<-NA
        current_front<-NA
      }
      
      if(d1>as.Date("2003-01-05")) {
        destFile2 <- paste("./MH1chla_d_",day_start,".nc",sep="")
        if((!file.exists(destFile2)) | (file.info(destFile2)$size<1000)) {
          myURL2 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla8day.nc?chlorophyll[(',day_start,'T12:00:00Z):1:(',day_start,'T12:00:00Z)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
          junk2 <- GET(myURL2, write_disk(destFile2,overwrite=TRUE), progress())
        }
        current_chla<-brick(destFile2)
      } else {
        current_chla<-NA
      }
      
      if(d1 > as.Date("2016-04-17") & d1 < as.Date("2018-11-21")) {
        destFile3 <- paste("./HYCOM_v_d_",day_start,".nc",sep="")
        if((!file.exists(destFile3)) | (file.info(destFile3)$size<1000)) {
          myURL3 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/nrlHycomGLBu008e912D_LonPM180.nc?water_v[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(0.0):1:(0.0)][(',min_Lat,'):1:(',max_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
          junk3 <- GET(myURL3, write_disk(destFile3,overwrite=TRUE), progress())
        }
        current_waterv<-brick(destFile3)
      } 
      
      #coordinate space for OSCARV is all screwed up from 20-420 degrees--eliminate    
      # if(d1 > as.Date("1992-10-21") & d1 < as.Date("2018-06-12")) {
      #   destFile4 <- paste("./OSCAR_v_d_",day_start,".nc",sep="")
      #   if((!file.exists(destFile4)) | (file.info(destFile4)$size<1000)) {
      #     myURL4 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplOscar.nc?v[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(15.0):1:(15.0)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long+360,'):1:(',max_Long+360,')]',sep="")
      #     junk4 <- GET(myURL4, write_disk(destFile4,overwrite=TRUE), progress())
      #   }
      #   current_oscarv<-brick(destFile4)
      # } 
      
      if(d1>as.Date("2003-01-05")) {
        destFile5 <- paste("./MODIS_PP_d_",day_start,".nc",sep="")
        if((!file.exists(destFile5)) | (file.info(destFile5)$size<1000)) {
          myURL5 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1pp8day.nc?productivity[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(0.0):1:(0.0)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
          junk5 <- GET(myURL5, write_disk(destFile5,overwrite=TRUE), progress())
        }
        current_pp<-brick(destFile5)
      } else {
        current_pp<-NA
      }
      
      if(d1>as.Date("2012-01-05")) {
        destFile6 <- paste("./VIIRSN_k490_d_",day_start,".nc",sep="")
        if((!file.exists(destFile6)) | (file.info(destFile6)$size<1000)) {
          myURL6 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018k4908day.nc?k490[(',day_start,'T00:00:00Z):1:(',day_start,'T00:00:00Z)][(',max_Lat,'):1:(',min_Lat,')][(',min_Long,'):1:(',max_Long,')]',sep="")
          junk6 <- GET(myURL6, write_disk(destFile6,overwrite=TRUE), progress())
        }
        current_k490<-brick(destFile6)
      } 
      
      if(d1 >= as.Date("2013-01-01")) {
        destFile7 <- paste("./WAVE_d_",day_start,".nc",sep="")
        if((!file.exists(destFile7)) | (file.info(destFile7)$size<1000)) {
          myURL7 <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NWW3_Global_Best.nc?Thgt[(',day_start,'T09:00:00Z):1:(',day_start,'T09:00:00Z)][(0.0):1:(0.0)][(',min_Lat,'):1:(',max_Lat,')][(',min_Long+360,'):1:(',max_Long+360,')]',sep="")
          junk7 <- GET(myURL7, write_disk(destFile7,overwrite=TRUE), progress())
        }
        current_wave<-brick(destFile7)
      }
    }#end if
    
    #assign data to point
    
    #depth
    if(all_manta_pts[i,"Long"]>=min_Long & all_manta_pts[i,"Long"]<=max_Long & all_manta_pts[i,"Lat"]>=min_Lat & all_manta_pts[i,"Lat"]<=max_Lat) {
      all_manta_pts[i,"Depth_m"]<-get.depth(domain_bathy,x=all_manta_pts[i,"Long"],y=all_manta_pts[i,"Lat"],locator=FALSE)$depth
      if(d1>as.Date("2002-06-01")) {
        all_manta_pts[i,"SST"]<-raster::extract(current_sst,pt)
        all_manta_pts[i,"Front"]<-raster::extract(current_front,pt) #match the frontal gradient at the point to the correct day within the brick
        all_manta_pts[i,"Front_Z"]<-all_manta_pts[i,"Front"]/cellStats(current_front,stat="max") #match the frontal gradient at the point to the correct day within the brick
      } else {
        all_manta_pts[i,"SST"]<-NA
        all_manta_pts[i,"Front"]<-NA
        all_manta_pts[i,"Front_Z"]<-NA
      }
      if(d1>as.Date("2003-01-05")) {
        all_manta_pts[i,"ChlA"]<-raster::extract(current_chla,pt,method='bilinear') #to fill gaps
      } else {
        all_manta_pts[i,"ChlA"]<-NA
      }
      
      if(d1 > as.Date("2016-04-17") & d1 < as.Date("2018-11-21")) {
        all_manta_pts[i,"waterv"]<-raster::extract(current_waterv,pt)
      } else {all_manta_pts[i,"waterv"]<-NA}
      
      # if(d1 > as.Date("1992-10-21") & d1 < as.Date("2018-06-12")) {
      #   all_manta_pts[i,"oscarv"]<-raster::extract(current_oscarv,SpatialPoints(cbind((all_manta_pts[i,"Long"]+360),all_manta_pts[i,"Lat"])),method="bilinear") #to fill gaps
      # } else {all_manta_pts[i,"oscarv"]<-NA}
      
      if(d1>as.Date("2003-01-05")) {
        all_manta_pts[i,"pp"]<-raster::extract(current_pp,pt,method="bilinear") #to fill gaps
      } else {
        all_manta_pts[i,"pp"]<-NA
      }
      
      if(d1>as.Date("2012-01-05")) {
        all_manta_pts[i,"k490"]<-raster::extract(current_k490,pt,method="bilinear") #to fill gaps
      } else {all_manta_pts[i,"k490"]<-NA}
      
      if(d1 < as.Date("2013-01-01")) {
        all_manta_pts[i,"wave"]<-NA
      } else {
        all_manta_pts[i,"wave"]<-raster::extract(current_wave,SpatialPoints(cbind((all_manta_pts[i,"Long"]+360),all_manta_pts[i,"Lat"])))
      }
    }#within the data range
  }#is.na(d1)
}#end loop

write.csv(all_manta_pts,"C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/all_manta_pts_env.csv")

all_manta_pts_env<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/all_manta_pts_env.csv",stringsAsFactors = F)
summary(all_manta_pts_env)

#### 7.5 NYSERDA - Add Environmental Data to Grids ####
d.nyserda <- read.csv("./nyserda_track_env.csv",stringsAsFactors = F)

d.nyserda_spt_wgs2 <- SpatialPointsDataFrame(coords = d.nyserda[,c("Longitude","Latitude")], data = d.nyserda,
                                             proj4string = CRS("+init=epsg:4326"))

d.nyserda_spt_alb2<-spTransform(d.nyserda_spt_wgs2, CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

d.nyserda_spt_alb2@data$DfromShore<-NA
for(i in 1:nrow(d.nyserda_spt_alb2)) {
  d.nyserda_spt_alb2@data[i,"DfromShore"]<-gDistance(d.nyserda_spt_alb2[i,],coast.proj)  #units are in meters.
}
summary(d.nyserda_spt_alb2@data$DfromShore)

crs(d.nyserda_spt_alb2)
crs(domain)
proj4string(d.nyserda_spt_alb2)<-crs(domain)
proj4string(domain)<-crs(domain)
test<-over(d.nyserda_spt_alb2,domain)
d.nyserda_spt_alb2@data$DomainFID<-test$DomainFID
d.nyserda_spt_alb2@data$Survey_date<-paste(d.nyserda_spt_alb2@data$SurveyID,as.Date(d.nyserda_spt_alb2@data$Datetime),sep="_")
d.nyserda_spt_alb2@data$Survey_date_DomainFID<-paste(d.nyserda_spt_alb2@data$SurveyID,d.nyserda_spt_alb2@data$Survey,as.Date(d.nyserda_spt_alb2@data$Datetime),d.nyserda_spt_alb2@data$DomainFID,sep="_")
d.nyserda_spt_alb2@data$striparea<-d.nyserda_spt_alb2@data$effort
#set CRM depth and slope
writeOGR(d.nyserda_spt_alb2,dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA",layer="d.nyserda_spt_alb2",driver="ESRI Shapefile")

d.nyserda_spt_alb2<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA",layer="d.nyserda_spt_alb2")
d_all<-d.nyserda_spt_alb2@data
depth<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA",layer="nyserda_crmDepth")
depth@data<-dplyr::rename(depth@data,depth_crm=RASTERVALU)
test<-cbind(d_all,"depth_crm"=depth@data$depth_crm)

slope<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA",layer="nyserda_crmSlope")
slope@data<-dplyr::rename(slope@data,slope_crm=RASTERVALU)
test2<-cbind(test,"slope_crm"=slope@data$slope_crm)
test2$depth_crm<-ifelse(test2$depth_crm==-9999,NA,test2$depth_crm)
test2$slope_crm<-ifelse(test2$slope_crm==-9999,NA,test2$slope_crm)
summary(test2)

depth2<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA",layer="nyserda_gebDepth")
depth2@data<-dplyr::rename(depth2@data,depth_geb=RASTERVALU)
test3<-cbind(test2,"depth_geb"=depth2@data$depth_geb)
summary(test3)

slope2<-readOGR(dsn="C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA",layer="nyserda_gebSlope")
slope2@data<-dplyr::rename(slope2@data,slope_geb=RASTERVALU)
test4<-cbind(test3,"slope_geb"=slope2@data$slope_geb)
test4$depth_geb<-ifelse(test4$depth_geb==-9999,NA,test4$depth_geb)
test4$slope_geb<-ifelse(test4$slope_geb==-9999,NA,test4$slope_geb)
summary(test4)

test4$Depth_m<-ifelse(is.na(test4$depth_crm),test4$depth_geb,test4$depth_crm)
test4$Slope_deg<-ifelse(is.na(test4$slope_crm),test4$slope_geb,test4$slope_crm)

write.csv(test4,paste("./nyserda_env_20201216.csv",sep=""))

nyserda_env<-read.csv("C:/Users/nick.farmer/Documents/GIS/Data/Protected Species/Manta Ray/NYSERDA/nyserda_env_20201216.csv",stringsAsFactors=F)


Survey_date_DomainFID=summarize(nyserda_env$striparea,nyserda_env$Survey_date_DomainFID,sum,na.rm=T)[1]
striparea=summarize(nyserda_env$striparea,nyserda_env$Survey_date_DomainFID,sum,na.rm=T)
N=summarize(nyserda_env$N,nyserda_env$Survey_date_DomainFID,sum,na.rm=T)
Latitude=summarize(nyserda_env$Latitude,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
Longitude=summarize(nyserda_env$Longitude,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
Depth_m=summarize(nyserda_env$Depth_m,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
SST=summarize(nyserda_env$SST,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
Front=summarize(nyserda_env$Front,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
Front_Z=summarize(nyserda_env$Front_Z,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
pp=summarize(nyserda_env$pp,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
ChlA=summarize(nyserda_env$ChlA,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
waterv=summarize(nyserda_env$waterv,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
k490=summarize(nyserda_env$k490,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
wave=summarize(nyserda_env$wave,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)
Slope_deg=summarize(nyserda_env$Slope_deg,nyserda_env$Survey_date_DomainFID,max,na.rm=T)
DfromShore=summarize(nyserda_env$Slope_deg,nyserda_env$Survey_date_DomainFID,mean,na.rm=T)

nyserda_grid<-Reduce(function(x,y) merge(x = x, y = y, by = "nyserda_env$Survey_date_DomainFID"), 
                     list(Latitude,Longitude,striparea,N,Depth_m,SST,Front,Front_Z,pp,ChlA,waterv,k490,wave,Slope_deg,DfromShore))

colnames(nyserda_grid)<-c("Survey_date_DomainFID","Latitude","Longitude","striparea","N","Depth_m","SST","Front","Front_Z","pp","ChlA","waterv","k490","wave","Slope_deg","DfromShore")
nyserda_grid$Nhat<-nyserda_grid$N/nyserda_grid$striparea
nyserda_grid<-tidyr::separate(nyserda_grid,"Survey_date_DomainFID",c("SurveyTeam","YrSeason","Datetime","DomainFID"),sep="_",remove=F)
nyserda_grid$Survey<-paste(nyserda_grid$SurveyTeam,nyserda_grid$YrSeason,sep="_")
#remove records with invalid spatial data
nyserda_grid<-nyserda_grid[which(nyserda_grid$Latitude != 0 & nyserda_grid$Longitude != 0),]
nyserda_grid<-nyserda_grid[which(nyserda_grid$Depth_m<0), ] #only ocean obs.
nyserda_grid<-remove_all_labels(nyserda_grid)
nyserda_grid$stripareaKM<-nyserda_grid$striparea/(1000*1000)
summary(nyserda_grid)

write.csv(nyserda_grid,"./nyserda_grid_20201216.csv")