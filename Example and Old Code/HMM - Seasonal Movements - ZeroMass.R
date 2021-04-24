### momentuHMM vignette code link:
### https://github.com/bmcclintock/momentuHMM/tree/master/vignettes/examples

require(dplyr)
require(move)
require(momentuHMM)
require(lubridate)
library(rgdal)

#set working directory
setwd("E:/Maine Drive/Analysis/momentuHMM") #at home
setwd("C:/Users/Matt Gonnerman/Google Drive/Analysis/momentuHMM") #at work

listfemales <- read.csv("Trapping - Data.csv") %>%
  dplyr::select(ID = AlumBand, Sex) %>%
  filter(Sex == "F")

#################
### DATA PREP ###
#################

#Login to download data from Movebank directly
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

#### Load GPS data from turkeys
#2018
turkeygps2018 <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login,
                                 timestamp_start = "20180101000000000",timestamp_end = "20180801000000000")

#list of birds that made it through April 14 in 2018
April2018 <- c(280, 10376, 258, 260, 286, 361, 355, 354, 251)

raw.move.df2018 <-as.data.frame(coordinates(turkeygps2018)) %>%
  mutate(Timestamp = timestamps(turkeygps2018)) %>%
  mutate(Timestamp = round_date(Timestamp, unit = "hour")) %>%
  mutate(ID = paste(as.character(turkeygps2018@trackId), "2018", sep = "_")) %>%
  filter(ID %in% paste("X", April2018, "_2018", sep = "")) %>% #this limits the tracks to birds that survived through April 15
  filter(month(Timestamp) < 8 | month(Timestamp) > 10) #Limit to Nov-July, when taking hourly locations


#2019
turkeygps2019 <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login,
                                 timestamp_start = "20181101000000000", timestamp_end = "20190801000000000")

#list of birds that made it through April 14 in 2019
April2019 <- c(410, 409, 833, 426, 375, 393, 444, 447, 398, 436, 435, 420, 354)

raw.move.df2019 <-as.data.frame(coordinates(turkeygps2019)) %>%
  mutate(Timestamp = timestamps(turkeygps2019)) %>%
  mutate(Timestamp = round_date(Timestamp, unit = "hour")) %>%
  mutate(ID = paste(as.character(turkeygps2019@trackId), "2019", sep = "_")) %>%
  filter(ID %in% paste("X", April2019, "_2019", sep = "")) %>% #this limits the tracks to birds that survived through April 15
  filter(month(Timestamp) < 8 | month(Timestamp) > 10) #Limit to Nov-July, when taking hourly locations


#2020
turkeygps2020 <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login,
                                 timestamp_start = "20191101000000000", timestamp_end = "20200801000000000")

raw.move.df2020 <-as.data.frame(coordinates(turkeygps2020)) %>%
  mutate(Timestamp = timestamps(turkeygps2020)) %>%
  mutate(Timestamp = round_date(Timestamp, unit = "hour")) %>%
  mutate(ID = paste(as.character(turkeygps2020@trackId), "2020", sep = "_")) %>%
  filter(month(Timestamp) < 8 | month(Timestamp) > 10) #Limit to Nov-July, when taking hourly locations


raw.move.df <- raw.move.df2018
raw.move.df <- raw.move.df2019 # maybe for this only load birds that didn't die before nesting
raw.move.df <- raw.move.df2020
raw.move.df <- do.call("rbind", list(raw.move.df2018, raw.move.df2019, raw.move.df2020))

#Removes tracks that are very short (> 300 points pre roost additions)
checkcount <- raw.move.df %>%
  group_by(ID) %>%
  summarize(Total_Points = n())%>%
  filter(Total_Points > 300)
raw.move.df <- raw.move.df %>% filter(ID %in% unique(checkcount$ID))

#raw.move.df <- raw.move.df2020
# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(raw.move.df[,1:2], proj4string=CRS(projection(turkeygps2018)))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=19 ellps=WGS84"))

# add UTM locations to data frame
raw.move.df$x <- attr(utmcoord,"coords")[,1]
raw.move.df$y <- attr(utmcoord,"coords")[,2]

# Duplicate roost locations for all night time hours
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

#add a covariate for Julian Day
raw.move.df$YDay <- yday(raw.move.df$Timestamp)

#limits optimization from straying outside parameter bounds
ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)

#Use crawl function to fill in missing locations
Turkey.crawl.zm <- crawlWrap(raw.move.df, Time.name = "Timestamp", timeStep = "hour",
                          attempts=20, fixPar = c(NA, NA), prior = ln.prior, retryFits=10)
#plot(Turkey.crawl.zm)

# create momentuHMMData object from crwData object
turkeyData.zm <- prepData(data=Turkey.crawl.zm, covNames = "YDay")
#View(turkeyData.zm)

# add cosinor covariate based on hour of day
turkeyData.zm$hour <- as.numeric(strftime(turkeyData.zm$Timestamp, format = "%H", tz="GMT"))

#determine proportion of zero step lengths
whichzero <- which(turkeyData.zm$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(turkeyData.zm)

###########################################################
####### Stationary, Searching, and Dispersing Models ######
###########################################################
### MODEL 1
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 3 # Number of states
stateNames <- c("stationary", "localized", "dispersal") # label states
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

## constrain step length parameters: 
# Mean -> dispersal>searching>roosting
# SD -> no relationship
stepDM<-matrix(c(
                 1,0,0,0,0,0,0,0,0,
                 1,1,0,0,0,0,0,0,0,
                 1,1,1,0,0,0,0,0,0,
                 0,0,0,1,0,0,0,0,0,
                 0,0,0,0,1,0,0,0,0,
                 0,0,0,0,0,1,0,0,0,
                 0,0,0,0,0,0,1,0,0,
                 0,0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,0,0,1),
               nrow = 3*nbStates,byrow=TRUE,
               dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
                             c("mean_123:(Intercept)", "mean_2","mean_3",
                               paste0("sd_",1:nbStates,":(Intercept)"),
                               paste0("zero_",1:nbStates,":(Intercept)"))))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
                         dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(0,Inf,
                       0,Inf,
                       250, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       .2,.4,
                       .2,.4,
                       .2,.4), nrow = 3*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> searching < dispersal
angleDM<-matrix(c(1,0,0,
                  1,1,0,
                  1,1,1),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_1:(Intercept)","concentration_2","concentration_3")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,0.94, 
                        0,0.94, 
                        0,0.94),nrow = nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper"))) 
#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
                          nrow = ncol(angleDM), 
                          dimnames=list(colnames(angleDM),c("lower","upper")))


#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,
                 angle=angleworkBounds)
userBounds <- list(step = stepBounds, 
                   angle = angleBounds)

#prevents working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#Fixes "roost" for all, fixPar = fixPar in fitHMM
#fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))

# initial parameters
Par <- list(step=c(10, 150,1000, #mean
                   6, 131,203, #sd
                   .3, .3, .3), #zero mass
            angle = c(.01, 0.2 ,0.8))

Par0_m1.zm <- getParDM(data = turkeyData.zm,
                    nbStates = nbStates, 
                    dist = dist,
                    Par = Par, 
                    DM = DM, 
                    workBounds = workBounds, 
                    userBounds = userBounds,
                    estAngleMean = list(angle = FALSE))

# fit model 1
turk_m1.zm <- fitHMM(data = turkeyData.zm, 
                  nSims = nSims,
                  nbStates = nbStates, 
                  dist = dist, 
                  Par0 = Par0_m1.zm, 
                  DM = DM, 
                  workBounds = workBounds, 
                  userBounds = userBounds,
                  estAngleMean = list(angle=FALSE), 
                  prior = prior,
                  stateNames = stateNames)

turk_m1.zm
# plot(turk_m1.zm,legend.pos="right")
turk_m1_states <- viterbi(turk_m1.zm)


turkey_states <- turkeyData.zm
turkey_states$State <- turk_m1_states

turkey_states1 <- turkey_states %>%
  mutate(Date = as.Date(Timestamp)) %>%
  #filter(hour(Timestamp) %in% ) #Maybe consider filtering out nightime hours but will need to use suncalc
  group_by(ID, Date) %>%
  summarize(State_Day_Avg = mean(State))

require(ggplot2)

DailyStates <- ggplot(data=turkey_states1, aes(x = Date, y=State_Day_Avg)) +
  geom_point(size = 1.5) +
  xlab("Date") +
  ylab("Daily Average State") +
  theme_grey(base_size = 18) +
  geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1) + #best fit line
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
DailyStates

id = "X286_2018"
states_onebird <- turkey_states1 %>% filter( ID == id)
graph_onebird <- ggplot(data = states_onebird, aes(x = Date, y = State_Day_Avg)) +
  geom_point(size = 1.5) +
  labs(x = "Date", y = "State", title = id) +
  theme_grey(base_size = 18) +
  geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1)  #best fit line
graph_onebird  

#################################
####### Seasonal Breakdown ######
#################################
### MODEL 1
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("winter", "dispersal", "searching", "nesting") # label states
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

## constrain transition probabilities between seasons
fixbeta <- matrix(c(NA, NA, -100, #Winter
                    -100, NA, -100, #Dispersal
                    -100, -100, NA, #Searching
                    -100, -100, NA)) #Nesting

## constrain step length parameters: 
# Mean -> dispersal>searching>roosting
# SD -> no relationship
stepDM<-matrix(c(
  1,1,0,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
                c("mean_int:(Intercept)", "mean_WS","mean_D",
                  paste0("sd_",1:nbStates,":(Intercept)"),
                  paste0("zero_",1:nbStates,":(Intercept)"))))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf, 0,0, rep(-Inf,8),rep(Inf,11)),nrow = ncol(stepDM),
                         dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(50,Inf,
                       250,Inf,
                       50,Inf,
                       0, Inf,
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
# Concentration -> searching < dispersal
angleDM<-matrix(c(1,1,0,
                  1,1,1,
                  1,1,0,
                  1,0,0),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration:(Intercept)","concentration_WS", "concentration_D")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,0.94, 
                        0,0.94, 
                        0,0.94, 
                        0,0.94),nrow = nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper"))) 
#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0, rep(Inf,3)),nrow = ncol(angleDM), dimnames=list(colnames(angleDM),c("lower","upper")))


#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,angle=angleworkBounds)
userBounds <- list(step = stepBounds, angle = angleBounds)

#prevents working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#Fixes "roost" for all, fixPar = fixPar in fitHMM
#fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))

# initial parameters
Par <- list(step=c(100, 700, 100, 10, 25 ,150, 25, 5, .3, .3, .3, .3), angle = c(.3, 0.5 , 0.3, .01))

Par0_m1s.zm <- getParDM(data = turkeyData.zm,
                       nbStates = nbStates, 
                       dist = dist,
                       Par = Par, 
                       DM = DM, 
                       workBounds = workBounds, 
                       userBounds = userBounds,
                       estAngleMean = list(angle = FALSE))

# fit model 1
turk_m1s.zm <- fitHMM(data = turkeyData.zm, 
                     nSims = nSims,
                     nbStates = nbStates, 
                     dist = dist, 
                     Par0 = Par0_m1s.zm, 
                     DM = DM, 
                     workBounds = workBounds, 
                     userBounds = userBounds,
                     estAngleMean = list(angle=FALSE), 
                     prior = prior,
                     stateNames = stateNames,
                     fixPar=list(beta=fixbeta))

turk_m1s.zm

