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


# #2020
# turkeygps2020 <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login,
#                                  timestamp_start = "20191101000000000", timestamp_end = "20200801000000000")
# 
# raw.move.df2020 <-as.data.frame(coordinates(turkeygps2020)) %>%
#   mutate(Timestamp = timestamps(turkeygps2020)) %>%
#   mutate(Timestamp = round_date(Timestamp, unit = "hour")) %>%
#   mutate(ID = paste(as.character(turkeygps2020@trackId), "2020", sep = "_")) %>%
#   filter(month(Timestamp) < 8 | month(Timestamp) > 10) #Limit to Nov-July, when taking hourly locations


raw.move.df <- raw.move.df2018
raw.move.df <- raw.move.df2019
# raw.move.df <- raw.move.df2020
raw.move.df <- rbind(raw.move.df2018, raw.move.df2019)


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
    working.df$x[i] <- working.df$x[i] + rnorm(1,0,6.5) #change x UTM
    working.df$y[i] <- working.df$y[i] + rnorm(1,0,6.5) #change y UTM
  }
  for(i in 1:hr_before_sun){
    working.df$Timestamp[i+hr_after_sun] <- working.df$Timestamp[i+hr_after_sun] + (i*60*60) #change time
    working.df$x[i+hr_after_sun] <- working.df$x[i+hr_after_sun] + rnorm(1,0,6.5) #change x UTM
    working.df$y[i+hr_after_sun] <- working.df$y[i+hr_after_sun] + rnorm(1,0,6.5) #change y UTM
  }
  raw.move.df <- rbind(raw.move.df, working.df) #add to original database
}

raw.move.df <- raw.move.df %>% arrange(ID, Timestamp)

#add a covariate for Julian Day
raw.move.df$YDay <- yday(raw.move.df$Timestamp)

#remove the possibility of zero-mass step issues
for(i in 2:nrow(raw.move.df)){
  if(raw.move.df$x[i] == raw.move.df$x[i-1] & raw.move.df$y[i] == raw.move.df$y[i-1]){
    raw.move.df$x[i] <- raw.move.df$x[i] + 1
  }
}

#limits optimization from straying outside parameter bounds
ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)

#Use crawl function to fill in missing locations
Turkey.crawl <- crawlWrap(raw.move.df, Time.name = "Timestamp", timeStep = "hour",
                          attempts=20, fixPar = c(NA, NA), prior = ln.prior, retryFits=10)
#plot(Turkey.crawl)

# create momentuHMMData object from crwData object
turkeyData <- prepData(data=Turkey.crawl, covNames = "YDay")
#View(turkeyData)

# add cosinor covariate based on hour of day
turkeyData$hour <- as.numeric(strftime(turkeyData$Timestamp, format = "%H", tz="GMT"))


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
stepDM<-matrix(c(1,0,0,0,0,0,
                 1,1,0,0,0,0,
                 1,1,1,0,0,0,
                 0,0,0,1,0,0,
                 0,0,0,1,1,0,
                 0,0,0,1,1,1),2*nbStates,6,byrow=TRUE,
               dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates)),
                             c("mean_123:(Intercept)","mean_2","mean_3",
                               "sd_123:(Intercept)","sd_2","sd_3")))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf,0,0,-Inf,0,0,rep(Inf,6)),6,2,
                         dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(0,Inf,
                       0,Inf,
                       300, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf), nrow = 2*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> searching < dispersal
angleDM<-matrix(c(1,0,0,
                  1,1,0,
                  1,1,1),nbStates,3,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_123:(Intercept)","concentration_2","concentration_3")))
#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,Inf,Inf,Inf),3,2,dimnames=list(colnames(angleDM),c("lower","upper")))

#Userbound contraint on angle
angleBounds <- matrix(c(0,0.96, 
                        0,0.96, 
                        0,0.96),nrow=nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper"))) 

#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,angle=angleworkBounds)
userBounds <- list(step = stepBounds, angle = angleBounds)

#prevents working parameters from straing along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

# initial parameters
Par <- list(step=c(10, 150,700,6, 131,203), angle = c(.01, 0.4 ,0.7))
                    
Par0_m1 <- getParDM(nbStates = nbStates, dist = dist,
                 Par = Par, DM = DM, workBounds = workBounds, userBounds = userBounds,
                 estAngleMean = list(angle = FALSE))

# fit model 1
turk_m1 <- fitHMM(data = turkeyData, nSims = nSims, nbStates = nbStates, 
                  dist = dist, Par0 = Par0_m1, DM = DM, workBounds = workBounds, userBounds = userBounds,
                  estAngleMean = list(angle=FALSE), 
                  prior = prior,stateNames = stateNames)

turk_m1
# plot(turk_m1,legend.pos="right")
turk_m1_states <- viterbi(turk_m1)


turkey_states <- turkeyData
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

id = "X286_2018"
states_onebird <- turkey_states1 %>% filter( ID == id)
graph_onebird <- ggplot(data = states_onebird, aes(x = Date, y = State_Day_Avg)) +
  geom_point(size = 1.5) +
  labs(x = "Date", y = "State", title = id) +
  theme_grey(base_size = 18) +
  geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1)  #best fit line
graph_onebird  

###################################################################################
### MODEL 2
#include cosinor function to address daily patterns

formula <- ~ cosinor(hour, period = 24)

Par0_m2 <- getPar0(model = turk_m1, formula = formula, 
                   nbStates = nbStates, dist = dist,
                   Par = Par, DM = DM, workBounds = workBounds, userBounds = userBounds,
                   estAngleMean = list(angle = FALSE))

turk_m2 <- fitHMM(data = turkeyData, nSims = nSims, nbStates = nbStates, 
                  dist = dist, Par0 = Par0_m2$Par, DM = DM, workBounds = workBounds, userBounds = userBounds,
                  estAngleMean = list(angle=FALSE), 
                  prior = prior,stateNames = stateNames, formula = formula)


turk_m2
plot(turk_m2)



#############################################################################################
####### Seasonal Models - 4 States ##########################################################
#############################################################################################
#### Do not distinguish between summer and prenesting here.
#### Maybe in future see if you can use angle concentration D > P > W = S > N 
### MODEL 1
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("Winter", "Dispersal", "Prenesting/Summer", "Nesting") # label states
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

###NEED TO MERGE A JULIAN DAY COVARIATE ###

## constrain transition probabilities between seasons
fixbeta <- matrix(c(NA, NA, -100, #Winter
                    -100, NA, -100, #Dispersal
                    -100, -100, NA, #Prenesting
                    -100, NA, NA)) #Nesting

## constrain step length parameters: 
# Mean -> Nesting<Wintering<Prenesting=Summer<Dispersal
# SD -> no relationship
stepDM<-matrix(c(1,1,0,0,0,0,0,0,
                 1,1,1,1,0,0,0,0,
                 1,1,1,0,0,0,0,0,
                 1,0,0,0,0,0,0,0,
                 0,0,0,0,1,0,0,0,
                 0,0,0,0,0,1,0,0,
                 0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,0,1),2*nbStates,8,byrow=TRUE,
               dimnames=list(c(paste0("mean_",c("W", "D","PS","N")),
                               paste0("sd_",c("W", "D","PS","N"))),
                             c("mean_ALL:(Intercept)","mean_W","mean_PS", "mean_D", 
                               paste0("sd_",c("W", "D","PS","N"),":(Intercept)"))))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf,0,0,0,-Inf,-Inf,-Inf,-Inf,rep(Inf,8)),8,2,
                         dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(0,Inf,
                       600,Inf,
                       0,Inf,
                       0,Inf,
                       0,Inf,
                       0,Inf,
                       0, Inf,
                       0, Inf), nrow = 2*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> N<W=P=S<D
angleDM<-matrix(c(1,1,0,
                  1,1,1,
                  1,1,0,
                  1,0,0
                  ),nbStates,3,byrow=TRUE,
                dimnames=list(paste0("concentration_",c("W", "D","PS","N")),
                              c("concentration_N:(Intercept)","concentration_WPS","concentration_D")))

###MAYBE CAN DISTINGUISH SUMMER FROM PRENESTING USING CONCENTRATION

#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,Inf,Inf,Inf),3,2,dimnames=list(colnames(angleDM),c("lower","upper")))

#Userbound contraint on angle
angleBounds <- matrix(c(0,0.94, 
                        0,0.94, 
                        0,0.94, 
                        0,0.94),nrow=nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper"))) 

#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,
         angle=angleDM)
workBounds<-list(step=stepworkBounds,
                 angle=angleworkBounds)
userBounds <- list(step = stepBounds, 
                   angle = angleBounds)

#prevents working parameters from straing along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

# initial parameters
Par <- list(step=c(150, 1000,200,10,100,300,100,10), 
            angle = c(.2, 0.8 ,0.2, .01))

Par0_m1 <- getParDM(nbStates = nbStates,
                    dist = dist,
                    Par = Par, DM = DM, 
                    workBounds = workBounds, 
                    userBounds = userBounds,
                    estAngleMean = list(angle = FALSE))

# fit model 1
turk_seas_m1 <- fitHMM(data = turkeyData, 
                       nSims = nSims,
                       nbStates = nbStates, 
                       dist = dist, Par0 = Par0_m1, 
                       DM = DM, 
                       workBounds = workBounds, 
                       userBounds = userBounds,
                       estAngleMean = list(angle=FALSE),
                       fixPar = list(beta = fixbeta),
                       prior = prior,
                       stateNames = stateNames)


#plot(turk_seas_m1,legend.pos="right")










#####################################
### FIND BEST STARTING PARAMETERS ###
#####################################
#Recurring loop to find best starting values for initial parameters
niter <- 30
allm <- list()
allpar <- list()
for(i in 1:niter){
  stepMean0 <- runif(3, min = c(5, 30, 300), max = c(20, 200, 700))
  stepSD0 <- stepMean0/runif(3, min = .5, max = 3)
  
  angleCon0 <- runif(3, min = c(0.01, .1, .6), max = c(.1, .6, .90))
  
  
  stepPar0 <- c(stepMean0, stepSD0) 
  anglePar0 <- c(angleCon0) 
  
  Par <- list(step=stepPar0, angle = anglePar0)
  
  Par0_m1 <- getParDM(nbStates = nbStates, dist = dist,
                      Par = Par, DM = DM, workBounds = workBounds,
                      estAngleMean = list(angle = FALSE))

  # fit null model
  allm[[i]] <- fitHMM(data = turkeyData, nSims = nSims, ncores = ncores, nbStates = nbStates, 
                    dist = dist, Par0 = Par0_m1, DM = DM, workBounds = workBounds,
                    estAngleMean = list(angle=FALSE), 
                    prior = prior,stateNames = stateNames)
  allpar[[i]] <- Par
  print(i)

}
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
whichbest <- which.min(allnllk)
mbest <- allm[[whichbest]] 
mbest


#inster below in fitHMM to try and get a better fit through multiple tries.
#retryFits = retryFits, 


# I am pretty sure I can limit the seasonal movements purely by the transition probabilities.
# Only real issue would be prenesting->nesting and nesting->summer vs prenesting
