require(dplyr)
require(move)
require(momentuHMM)
require(lubridate)
library(rgdal)

#Login to download data from Movebank directly
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

### Download through 1st nest ###
#Only use birds that survived through April 14
April2018 <- c(286,355,354,251,361)
#Format as nest ids in database
nestids_2018 <- paste(April2018, 2018, 1, sep = "-") 

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
birdlist <- c(April2018)
nestlist <- c(nestids_2018)

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
llcoord <- SpatialPoints(raw.move.df[,1:2], proj4string=CRS(projection(indbird_mbdown)))
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

# add a covariate for Ordinal Day
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
