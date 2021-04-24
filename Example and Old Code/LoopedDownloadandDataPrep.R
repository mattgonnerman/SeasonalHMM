require(dplyr)
require(move)
require(momentuHMM)
require(lubridate)
library(rgdal)

#set working directory
setwd("E:/Maine Drive/Analysis/momentuHMM") #at home
setwd("C:/Users/Matt Gonnerman/Google Drive/Analysis/momentuHMM") #at work

### Download through 1st nest ###

#List of birds that survived through April 14 2018
#unused birds = 258, 260, 259
April2018 <- c(280, 10376, 286, 361, 355, 354, 251)
#Format as nest ids in database
nestids_2018 <- paste(April2018, 2018, 1, sep = "-") 

#list of birds that made it through April 14 2019
April2019 <- c(410, 409, 833, 426, 375, 393, 444, 447, 398, 436, 435, 420, 354)
#Format as nest ids in database
nestids_2019 <- paste(April2019, 2019, 1, sep = "-") 


#Load Nest Monitoring Database
nestmonitor.raw <- read.csv("Nest Monitoring - Nest Status Checks.csv") 
nestdetails.raw <- read.csv("Nest Monitoring - Nest Info.csv")

### Get Start and End Dates for movebank download based on nest database
#Start is Nov 1 for year prior to nest
#End is last day bird was on nest
lastonDates <- nestmonitor.raw %>%
  dplyr::select(-Comment) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(as.numeric(substr(NestID, nchar(as.character(NestID)), nchar(as.character(NestID)))) == 1) %>%
  filter(Status == 0) %>%
  group_by(NestID) %>%
  summarize(LastOnNest = max(Date)) %>%
  mutate(MovebankStartDate = paste(year(LastOnNest)-1,"1101000000000",sep ="" )) %>%
  mutate(MovebankEndDate = paste(gsub("\\D","", LastOnNest), "000000000", sep = ""))
lastonDates <- as.data.frame(lastonDates)


### CHOOSE YOUR DATASET
# Set which years you want to use with this code
birdlist <- c(April2018, April2019)
nestlist <- c(nestids_2018, nestids_2019)

### Need a single dataframe for the crawlWrap function
# loop through the birds, download from movebank and combine into database for crawl
raw.move.df <- data.frame(location_long = numeric(),
                          location_lat = numeric(),
                          Timestamp = as.POSIXct(character()),
                          ID = character())
for(i in 1:length(birdlist)){
  animalname <- as.character(birdlist[i])
  nestid <- nestlist[i]
  if(nestid %in% unique(lastonDates$NestID)){
      timestart <- lastonDates[lastonDates$NestID==nestid,3]
      timeend <- lastonDates[lastonDates$NestID==nestid,4]
  }else{ #not all birds nested, in which case just use full track from Nov11-July31
    timestart <- paste(as.numeric(substr(nestid, nchar(nestid)-5, nchar(nestid)-2))-1, "1101000000000", sep = "")
    timeend <- paste(substr(nestid, nchar(nestid)-5, nchar(nestid)-2), "0801000000000", sep = "")
  }

  
  yearnest <- substr(timeend, 0, 4)

  indbird_mbdown <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", 
                                    login = login,
                                    animal = animalname,
                                    timestamp_start = timestart,
                                    timestamp_end = timeend)
  indbird.df <-as.data.frame(coordinates(indbird_mbdown)) %>%
    mutate(Timestamp = timestamps(indbird_mbdown)) %>%
    mutate(Timestamp = round_date(Timestamp, unit = "hour")) %>%
    mutate(ID = paste(animalname, yearnest, sep = "_"))
  
  raw.move.df <- rbind(raw.move.df, indbird.df)
}

### Need to have coordinates in UTM not LL
# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(raw.move.df[,1:2], proj4string=CRS(projection(turkeygps2018)))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=19 ellps=WGS84"))

# add UTM locations to data frame
raw.move.df$x <- attr(utmcoord,"coords")[,1]
raw.move.df$y <- attr(utmcoord,"coords")[,2]


### Duplicate roost locations for all night time hours
# includes random error based on gps location uncertainty (~12.7m)
roosts <- raw.move.df %>%
  rename(lon = location_long, lat = location_lat) %>%
  mutate(date = as.Date(Timestamp)) %>%
  filter(hour(Timestamp) < 6 & hour(Timestamp) > 3)

require(suncalc)
suntimes <- getSunlightTimes(data = roosts, tz = "GMT", keep = c("sunrise", "sunset"))
roosts$sunrise <- suntimes$sunrise
roosts$sunset <- suntimes$sunset

for(j in 1:nrow(roosts)){
  roosttime <- roosts$Timestamp[j]
  hr_before_sun <- floor(as.numeric(roosts$sunrise[j]-roosts$Timestamp[j]))-1
  hr_after_sun <- floor(24 - as.numeric(roosts$sunset[j]-roosts$Timestamp[j]))-1
  
  working.df <- roosts[rep(j, hr_before_sun + hr_after_sun),] %>%
    dplyr::select(location_long = lon, location_lat = lat, Timestamp, ID, x, y)
  
  for(i in 1:hr_after_sun){
    working.df$Timestamp[i] <- working.df$Timestamp[i] - (i*60*60) #change time
    #working.df$x[i] <- working.df$x[i] + rnorm(1,0,6.5) #change x UTM
    #working.df$y[i] <- working.df$y[i] + rnorm(1,0,6.5) #change y UTM
  }
  for(i in 1:hr_before_sun){
    working.df$Timestamp[i+hr_after_sun] <- working.df$Timestamp[i+hr_after_sun] + (i*60*60) #change time
    #working.df$x[i+hr_after_sun] <- working.df$x[i+hr_after_sun] + rnorm(1,0,6.5) #change x UTM
    #working.df$y[i+hr_after_sun] <- working.df$y[i+hr_after_sun] + rnorm(1,0,6.5) #change y UTM
  }
  raw.move.df <- rbind(raw.move.df, working.df) #add to original database
}

raw.move.df <- raw.move.df %>% arrange(ID, Timestamp)

# add a covariate for Julian Day
raw.move.df$YDay <- yday(raw.move.df$Timestamp)

# limits optimization from straying outside parameter bounds
ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)

# Create dataframe of nest locations to use with center
nestlocations_dd <- nestdetails.raw %>% 
  filter(Nest.Attempt == 1) %>%
  dplyr::select(BirdID = Alum.Band.ID, NestID, NestLat, NestLong) %>%
  filter(!is.na(NestLat))
nestllcoord <- SpatialPoints(nestlocations_dd[,4:3], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"))
nestutmcoord <- spTransform(nestllcoord,CRS("+proj=utm +zone=19 ellps=WGS84"))
# add UTM locations to data frame
nestlocations_dd$x <- attr(nestutmcoord,"coords")[,1]
nestlocations_dd$y <- attr(nestutmcoord,"coords")[,2]
nestlocations_dd <- nestlocations_dd %>%
  mutate(NestID = as.character(NestID)) %>%
  mutate(ID = paste(BirdID, substr(NestID, nchar(NestID)-5, nchar(NestID)-2), sep="_"))


### Need crawlWrap and prepData individual for each bird to have individual specific centers
crawllist <- unique(raw.move.df$ID)
remove(wintertonest1_prepped)
for(i in 1:length(crawllist)){
  crawlid <- crawllist[i]
  bird1 <- raw.move.df %>% filter(ID == crawlid)
  nestx <- nestlocations_dd[nestlocations_dd$ID == crawlid, 5]
  nesty <- nestlocations_dd[nestlocations_dd$ID == crawlid, 6]
  
  #Estimate mean center of winter points
  winterpoints <- bird1 %>% filter(YDay < 92 | YDay > 300)
  wintercenterx <- round(mean(winterpoints$x))
  wintercentery <- round(mean(winterpoints$y))
  
  #Create matrix with winter center and 1st nest location
  if(length(nestx)==0){ #If a bird did not nest
    nonnestpoints <- bird1 %>% filter(YDay > 121 & YDay < 214) #subset days from May1 through July 31
    nestx <- round(mean(nonnestpoints$x)) #average points to get centroid
    nesty <- round(mean(nonnestpoints$y))
  }
  centers <- matrix(c(wintercenterx, wintercentery, 
                     nestx, nesty), nrow = 2, byrow = T,
                   dimnames=list(c("wintercenter", "nest1"),c("x","y")))
  
  
  bird1crawl <- crawlWrap(bird1, Time.name = "Timestamp", timeStep = "hour",
                          attempts=30, fixPar = c(NA, NA), prior = ln.prior, retryFits=20)
  bird1prepped <- prepData(data=bird1crawl, center = centers, covNames = "YDay")
  
  if(exists("wintertonest1_prepped")){
    wintertonest1_prepped <- rbind(wintertonest1_prepped, bird1prepped)
  }else{
    wintertonest1_prepped <- bird1prepped
  }
}

#check output
summary(wintertonest1_prepped)
#If there was an error in crawlWrap that ended loop, can just set start in loop to 
#1+total observations in summary
#just make sure to set it back when you are done


# ###########################################################
# ####### Stationary, Searching, and Dispersing Models ######
# ###########################################################
# ### MODEL 1 - base model, no covariates
# ### momentuHMM Function Objects
# nSims <- 100 # number of imputatons
# retryFits <- 10 # number attempt to re-fit based on random perturbation
# nbStates <- 3 # Number of states
# stateNames <- c("stationary", "localized", "dispersal") # label states
# dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes
# 
# ## constrain step length parameters: 
# # Mean -> dispersal>searching>roosting
# # SD -> no relationship
# stepDM<-matrix(c(
#   1,0,0,0,0,0,0,0,0,
#   1,1,0,0,0,0,0,0,0,
#   1,1,1,0,0,0,0,0,0,
#   0,0,0,1,0,0,0,0,0,
#   0,0,0,0,1,0,0,0,0,
#   0,0,0,0,0,1,0,0,0,
#   0,0,0,0,0,0,1,0,0,
#   0,0,0,0,0,0,0,1,0,
#   0,0,0,0,0,0,0,0,1),
#   nrow = 3*nbStates,byrow=TRUE,
#   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
#                 c("mean_123:(Intercept)", "mean_2","mean_3",
#                   paste0("sd_",1:nbStates,":(Intercept)"),
#                   paste0("zero_",1:nbStates,":(Intercept)"))))
# 
# #define the directions of the differences
# stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
#                          dimnames=list(colnames(stepDM),c("lower","upper")))
# #Userbound constraint on step
# stepBounds <- matrix(c(0,Inf,
#                        0,Inf,
#                        250, Inf,
#                        0, Inf,
#                        0, Inf,
#                        0, Inf,
#                        .2,.4,
#                        .2,.4,
#                        .2,.4), nrow = 3*nbStates, byrow = T,
#                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# 
# ## constrain turning angle concentration parameters: 
# # Concentration -> searching < dispersal
# angleDM<-matrix(c(1,0,0,
#                   1,1,0,
#                   1,1,1),nrow = nbStates,byrow=TRUE,
#                 dimnames=list(paste0("concentration_",1:nbStates),
#                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# 
# #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# #Userbound contraint on angle
# angleBounds <- matrix(c(0,0.94, 
#                         0,0.94, 
#                         0,0.94),nrow = nbStates,
#                       byrow=TRUE, dimnames=list(rownames(angleDM),
#                                                 c("lower","upper"))) 
# #define direction of differences
# angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
#                           nrow = ncol(angleDM), 
#                           dimnames=list(colnames(angleDM),c("lower","upper")))
# 
# 
# #Bundle individual parameter DM and workbounds
# DM<-list(step=stepDM,angle=angleDM)
# workBounds<-list(step=stepworkBounds,
#                  angle=angleworkBounds)
# userBounds <- list(step = stepBounds, 
#                    angle = angleBounds)
# 
# #prevents working parameters from straying along boundary
# prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}
# 
# #Fixes "roost" for all, fixPar = fixPar in fitHMM
# #fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))
# 
# # initial parameters
# Par <- list(step=c(10, 150,1000, #mean
#                    6, 131,203, #sd
#                    .3, .3, .3), #zero mass
#             angle = c(.01, 0.2 ,0.8))
# 
# 
# Par0.wtn1 <- getParDM(data = wintertonest1_prepped,
#                        nbStates = nbStates, 
#                        dist = dist,
#                        Par = Par, 
#                        DM = DM, 
#                        workBounds = workBounds, 
#                        userBounds = userBounds,
#                        estAngleMean = list(angle = FALSE))
# 
# # fit model 1
# turkm1.wtn1 <- fitHMM(data = wintertonest1_prepped, 
#                      nSims = nSims,
#                      nbStates = nbStates, 
#                      dist = dist, 
#                      Par0 = Par0.wtn1, 
#                      DM = DM, 
#                      workBounds = workBounds, 
#                      userBounds = userBounds,
#                      estAngleMean = list(angle=FALSE), 
#                      prior = prior,
#                      stateNames = stateNames)
# 
# turkm1.wtn1


##############################################################################################
#################################
####### Seasonal Breakdown ######
#################################
### MODEL 1
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("winter", "dispersal", "prenesting", "nesting") # label states
state_abb <- c("W", "D", "P", "N") #abbreviations for column names
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

## Constrain transition probabilities between seasons
# when you don't allow winter to go directly to prenesting...(run this now and see what happens)

fixbeta <- matrix(c(NA, NA, -100, #Winter
                    -100, NA, -100, #Dispersal
                    -100, -100, NA, #Searching
                    -100, -100, -100)) #Nesting

## constrain step length parameters: 
# Mean -> W < D, N < P < D
# SD -> P < D
# Zero -> no relationship
stepDM<-matrix(c(
  1,1,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,0,1,1,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",state_abb),paste0("sd_",state_abb), paste0("zero_",state_abb)),
                c("mean_N:(Intercept)", "mean_W","mean_D", "mean_P",
                  "sd_W", "sd_P:(Intercept)","sd_D", "sd_N",
                  paste0("zero_",state_abb,":(Intercept)"))))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf,0,0,0,
                           -Inf, -Inf, 0, -Inf,
                            rep(-Inf,4),
                            rep(Inf,12)),nrow = ncol(stepDM),
                          dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(100,Inf, #100,250,175,5 has given best results
                        250,Inf,
                        175,Inf,
                        5, Inf,
                        0, Inf,
                        0, Inf,
                        0, Inf,
                        0, Inf,
                        .2,.4,
                        .2,.4,
                        .2,.4,
                        .2,.4), nrow = 3*nbStates, byrow = T,
                      dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> N < W < P < D
angleDM<-matrix(c(1,1,0,0,
                  1,1,1,1,
                  1,1,1,0,
                  1,0,0,0),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration:(Intercept)","concentration_W", "concentration_P", "concentration_D")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,0.5, 
                        0,0.5, 
                        0,0.5, 
                        0,0.5),nrow = nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper"))) 
#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,0,
                            rep(Inf,4)),nrow = ncol(angleDM), dimnames=list(colnames(angleDM),c("lower","upper")))


#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,angle=angleworkBounds)
userBounds <- list(step = stepBounds, angle = angleBounds)

#prevents working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#Fixes "roost" for all, fixPar = fixPar in fitHMM
#fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))

# initial distribution
delta0 <- matrix(c(.99997, .00001, .00001, .00001),1,
                 dimnames = list(NULL,c("winter", "dispersal","prenesting", "nesting")))

# initial parameters
Par <- list(step=c(150,700,300,20, 
                   100,300,100,5, 
                   .3, .3, .3, .3),
            angle = c(.2, 0.45 , 0.35, .01))

Par0_m1s.wtn <- getParDM(data = wintertonest1_prepped,
                        nbStates = nbStates, 
                        dist = dist,
                        Par = Par, 
                        DM = DM, 
                        workBounds = workBounds, 
                        userBounds = userBounds,
                        estAngleMean = list(angle = FALSE)) #Try changing this to TRUE

# fit model 1
turk_m1s.wtn <- fitHMM(data = wintertonest1_prepped, 
                      nSims = nSims,
                      nbStates = nbStates, 
                      dist = dist, 
                      Par0 = Par0_m1s.wtn, 
                      DM = DM, 
                      workBounds = workBounds, 
                      userBounds = userBounds,
                      estAngleMean = list(angle=FALSE), #Try changing this to TRUE
                      prior = prior,
                      stateNames = stateNames,
                      fixPar=list(beta=fixbeta),
                      delta0 = delta0)


turk_m1s.wtn
states1 <- viterbi(turk_m1s.wtn)
table(states1)/nrow(wintertonest1_prepped)
plot(turk_m1s.wtn)

############################################
### DO NOT CHANGE ABOVE FOR CHRIST'S SAKE ###
############################################





############################################################################################
############################################################################################
### MODEL 2 - Include centers of movement
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("winter", "dispersal", "prenesting", "nesting") # label states
state_abb <- c("W", "D", "P", "N") #abbreviations for column names
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

## Constrain transition probabilities between seasons
# when you don't allow winter to go directly to prenesting...(run this now and see what happens)

fixbeta <- matrix(rep(c(NA, NA, -100, #Winter
                    -100, NA, -100, #Dispersal
                    -100, -100, NA, #Searching
                    -100, -100, -100),2), 2, byrow = T) #Nesting

## constrain step length parameters: 
# Mean -> W < D, N < P < D
# SD -> P < D
# Zero -> no relationship
stepDM2 <- matrix(c(
  1,1,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,0,1,1,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",state_abb),paste0("sd_",state_abb), paste0("zero_",state_abb)),
                c("mean_N:(Intercept)", "mean_W","mean_D", "mean_P",
                  "sd_W", "sd_P:(Intercept)","sd_D", "sd_N",
                  paste0("zero_",state_abb,":(Intercept)"))))

#define the directions of the differences
stepworkBounds2 <- matrix(c(-Inf,0,0,0,
                           -Inf, -Inf, 0, -Inf,
                           rep(-Inf,4),
                           rep(Inf,12)),nrow = ncol(stepDM2),
                         dimnames=list(colnames(stepDM2),c("lower","upper")))
#Userbound constraint on step
stepBounds2 <- matrix(c(100,Inf, #100,250,175,5 has given best results
                       250,Inf,
                       175,Inf,
                       5, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       .2,.4,
                       .2,.4,
                       .2,.4,
                       .2,.4), nrow = 3*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM2), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> N < W < P < D
# angleDM2 <- matrix(c(1,1,0,0,
#                   1,1,1,1,
#                   1,1,1,0,
#                   1,0,0,0),nrow = nbStates,byrow=TRUE,
#                 dimnames=list(paste0("concentration_",1:nbStates),
#                               c("concentration:(Intercept)","concentration_W", "concentration_P", "concentration_D")))

angleDM2 <-   matrix(c("wintercenter.angle",0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,
                     0,0,"nest1.angle",0,0,0,0,0,
                     0,0,0,"nest1.angle",0,0,0,0,
                     0,0,0,0,1,1,0,0,
                     0,0,0,0,1,1,1,1,
                     0,0,0,0,1,1,1,0,
                     0,0,0,0,1,0,0,0),nrow = nbStates*2,byrow=TRUE,
                   dimnames=list(c(paste0("mean_",state_abb), paste0("concentration_",state_abb)),
                                 c("mean_W","mean_D", "mean_P", "mean_N",
                                   "concentration:(Intercept)","concentration_W", "concentration_P", "concentration_D")))


#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
# angleBounds2 <- matrix(c(0,0.5, 
#                         0,0.5, 
#                         0,0.5, 
#                         0,0.5),nrow = nbStates,
#                       byrow=TRUE, dimnames=list(rownames(angleDM2),
#                                                 c("lower","upper"))) 

angleBounds2 <- matrix(c(-pi, pi,
                         -pi, pi,
                         -pi, pi,
                         -pi, pi,
                         0,0.5, 
                         0,0.5, 
                         0,0.5, 
                         0,0.5),nrow = nbStates*2,
                       byrow=TRUE, dimnames=list(rownames(angleDM2),
                                                 c("lower","upper"))) 
#define direction of differences
# angleworkBounds2 <- matrix(c(-Inf,0,0,0,
#                             rep(Inf,4)),nrow = ncol(angleDM2), dimnames=list(colnames(angleDM2),c("lower","upper")))
angleworkBounds2 <- matrix(c(rep(-Inf,5),0,0,0,
                             rep(Inf,8)),nrow = ncol(angleDM2), dimnames=list(colnames(angleDM2),c("lower","upper")))


#Bundle individual parameter DM and workbounds
DM2 <- list(step=stepDM2,angle=angleDM2)
workBounds2 <- list(step=stepworkBounds2,angle=angleworkBounds2)
userBounds2 <- list(step = stepBounds2, angle = angleBounds2)

#prevents working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#Fixes "roost" for all, fixPar = fixPar in fitHMM
#fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))

# initial distribution
fixdelta <- exp(c(100,rep(0,3)))/sum(exp(c(100,rep(0,3))))

# Formula for transition probabilities
transFormula <- formula("~toState3(I(nest1.dist<1000))")

#Set initial state to W
aInd <- NULL
for(id in unique(wintertonest1_prepped$ID)) {
  idInd <- which(wintertonest1_prepped$ID==id)
  aInd <- c(aInd,idInd[1])
}
knownStates <- rep(NA,nrow(wintertonest1_prepped)) 
knownStates[aInd] <- 1


# initial parameters
Par2 <- list(step=c(150,700,300,20, 
                   100,300,100,5, 
                   .3, .3, .3, .3),
            angle = c(0,0,0,0,
                      .2, 0.45 , 0.35, .01))

Par0_m2s.wtn <- getParDM(data = wintertonest1_prepped,
                         nbStates = nbStates, 
                         dist = dist,
                         Par = Par2, 
                         DM = DM2, 
                         workBounds = workBounds2, 
                         userBounds = userBounds2,
                         estAngleMean = list(angle = TRUE),
                         circularAngleMean=list(angle=TRUE)) #Try changing this to TRUE

# fit model 1
turk_m2s.wtn <- fitHMM(data = wintertonest1_prepped, 
                       nSims = nSims,
                       nbStates = nbStates, 
                       dist = dist, 
                       Par0 = Par0_m2s.wtn, 
                       DM = DM2, 
                       workBounds = workBounds2, 
                       userBounds = userBounds2,
                       estAngleMean = list(angle=TRUE), #Try changing this to TRUE
                       circularAngleMean=list(angle=TRUE),
                       prior = prior,
                       stateNames = stateNames,
                       formula = transFormula,
                       fixPar=list(beta=fixbeta, delta = fixdelta),
                       knownStates = knownStates)


turk_m2s.wtn
states2 <- viterbi(turk_m2s.wtn)
table(states2)/nrow(wintertonest1_prepped)
plot(turk_m2s.wtn)

############################################
### DO NOT CHANGE ABOVE FOR CHRIST'S SAKE ###
############################################





############################################################################################
############################################################################################
### MODEL 3 - include Individual random effects
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("winter", "dispersal", "prenesting", "nesting") # label states
state_abb <- c("W", "D", "P", "N") #abbreviations for column names
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

## Constrain transition probabilities between seasons
# when you don't allow winter to go directly to prenesting...(run this now and see what happens)

fixbeta <- matrix(rep(c(NA, NA, -100, #Winter
                        -100, NA, -100, #Dispersal
                        -100, -100, NA, #Searching
                        -100, -100, -100),2), 2, byrow = T) #Nesting

## constrain step length parameters: 
# Mean -> W < D, N < P < D
# SD -> P < D
# Zero -> no relationship
stepDM3 <- matrix(c(
  1,1,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,0,1,1,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",state_abb),paste0("sd_",state_abb), paste0("zero_",state_abb)),
                c("mean_N:(Intercept)", "mean_W","mean_D", "mean_P",
                  "sd_W", "sd_P:(Intercept)","sd_D", "sd_N",
                  paste0("zero_",state_abb,":(Intercept)"))))

#define the directions of the differences
stepworkBounds3 <- matrix(c(-Inf,0,0,0,
                            -Inf, -Inf, 0, -Inf,
                            rep(-Inf,4),
                            rep(Inf,12)),nrow = ncol(stepDM3),
                          dimnames=list(colnames(stepDM3),c("lower","upper")))
#Userbound constraint on step
stepBounds3 <- matrix(c(100,Inf, #100,250,175,5 has given best results
                        250,Inf,
                        175,Inf,
                        5, Inf,
                        0, Inf,
                        0, Inf,
                        0, Inf,
                        0, Inf,
                        .2,.4,
                        .2,.4,
                        .2,.4,
                        .2,.4), nrow = 3*nbStates, byrow = T,
                      dimnames = list(rownames(stepDM3), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> N < W < P < D
# angleDM3 <- matrix(c(1,1,0,0,
#                   1,1,1,1,
#                   1,1,1,0,
#                   1,0,0,0),nrow = nbStates,byrow=TRUE,
#                 dimnames=list(paste0("concentration_",1:nbStates),
#                               c("concentration:(Intercept)","concentration_W", "concentration_P", "concentration_D")))

angleDM3 <-   matrix(c("wintercenter.angle",0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,
                       0,0,"nest1.angle",0,0,0,0,0,
                       0,0,0,"nest1.angle",0,0,0,0,
                       0,0,0,0,1,1,0,0,
                       0,0,0,0,1,1,1,1,
                       0,0,0,0,1,1,1,0,
                       0,0,0,0,1,0,0,0),nrow = nbStates*2,byrow=TRUE,
                     dimnames=list(c(paste0("mean_",state_abb), paste0("concentration_",state_abb)),
                                   c("mean_W","mean_D", "mean_P", "mean_N",
                                     "concentration:(Intercept)","concentration_W", "concentration_P", "concentration_D")))


#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
# angleBounds3 <- matrix(c(0,0.5, 
#                         0,0.5, 
#                         0,0.5, 
#                         0,0.5),nrow = nbStates,
#                       byrow=TRUE, dimnames=list(rownames(angleDM3),
#                                                 c("lower","upper"))) 

angleBounds3 <- matrix(c(-pi, pi,
                         -pi, pi,
                         -pi, pi,
                         -pi, pi,
                         0,0.5, 
                         0,0.5, 
                         0,0.5, 
                         0,0.5),nrow = nbStates*2,
                       byrow=TRUE, dimnames=list(rownames(angleDM3),
                                                 c("lower","upper"))) 
#define direction of differences
# angleworkBounds3 <- matrix(c(-Inf,0,0,0,
#                             rep(Inf,4)),nrow = ncol(angleDM3), dimnames=list(colnames(angleDM3),c("lower","upper")))
angleworkBounds3 <- matrix(c(rep(-Inf,5),0,0,0,
                             rep(Inf,8)),nrow = ncol(angleDM3), dimnames=list(colnames(angleDM3),c("lower","upper")))


#Bundle individual parameter DM and workbounds
DM3 <- list(step=stepDM3,angle=angleDM3)
workBounds3 <- list(step=stepworkBounds3,angle=angleworkBounds3)
userBounds3 <- list(step = stepBounds3, angle = angleBounds3)

#prevents working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#Fixes "roost" for all, fixPar = fixPar in fitHMM
#fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))

# initial distribution
fixdelta <- exp(c(100,rep(0,3)))/sum(exp(c(100,rep(0,3))))

# Formula for transition probabilities
transFormula <- formula("~toState3(I(nest1.dist<500))")

#Set initial state to W
aInd <- NULL
for(id in unique(wintertonest1_prepped$ID)) {
  idInd <- which(wintertonest1_prepped$ID==id)
  aInd <- c(aInd,idInd[1])
}
knownStates <- rep(NA,nrow(wintertonest1_prepped)) 
knownStates[aInd] <- 1


# initial parameters
Par3 <- list(step=c(150,700,300,20, 
                    100,300,100,5, 
                    .3, .3, .3, .3),
             angle = c(0,0,0,0,
                       .2, 0.45 , 0.35, .01))

Par0_m3s.wtn <- getParDM(data = wintertonest1_prepped,
                         nbStates = nbStates, 
                         dist = dist,
                         Par = Par3, 
                         DM = DM3, 
                         workBounds = workBounds3, 
                         userBounds = userBounds3,
                         estAngleMean = list(angle = TRUE),
                         circularAngleMean=list(angle=TRUE)) #Try changing this to TRUE

# fit model 1
turk_m3s.wtn <- fitHMM(data = wintertonest1_prepped, 
                       nSims = nSims,
                       nbStates = nbStates, 
                       dist = dist, 
                       Par0 = Par0_m3s.wtn, 
                       DM = DM3, 
                       workBounds = workBounds3, 
                       userBounds = userBounds3,
                       estAngleMean = list(angle=TRUE), #Try changing this to TRUE
                       circularAngleMean=list(angle=TRUE),
                       prior = prior,
                       stateNames = stateNames,
                       formula = transFormula,
                       fixPar=list(beta=fixbeta, delta = fixdelta),
                       knownStates = knownStates)


turk_m3s.wtn
states3 <- viterbi(turk_m3s.wtn)
table(states3)/nrow(wintertonest1_prepped)
plot(turk_m3s.wtn)


### Things left to add:
## Individual Heterogeneity?
## Covariates to explore
# Age
# Study Area