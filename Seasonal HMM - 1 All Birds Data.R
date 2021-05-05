### Load Necessary packages
lapply(c('dplyr', 'move', 'momentuHMM', 'lubridate', 'rgdal'), require, character.only = T)

### Login Credentials for MoveBanks
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

# NotIn operator for later
`%notin%` <- Negate(`%in%`)

### Bring in relevant databases
# Trapping
trap.raw <- read.csv("Trapping - Data.csv")

# Nest Monitoring
nestmonitor.raw <- read.csv("Nest Monitoring - Nest Status Checks.csv") 
nestdetails.raw <- read.csv("Nest Monitoring - Nest Info.csv")

#List of GPS Birds
gpshens <- trap.raw %>% 
  dplyr::select(BirdID = AlumBand, Trans.Type, Sex, Date) %>% 
  filter(Trans.Type == "GPS_Back", Sex == "F") %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
         Year = year(Date))
  

#List of 1st Nests
gpsnest1 <- nestdetails.raw %>% 
  dplyr::select(NestID, BirdID = Alum.Band.ID, Year, Nest.Attempt) %>% 
  filter(Nest.Attempt == 1) %>%
  filter(BirdID %in% gpshens$BirdID)


### Get Start and End Dates for movebank download based on nest database
#Start is Nov 1 for year prior to nest
#End is last day bird was on nest
lastonDates <- merge(nestmonitor.raw, gpsnest1, by = "NestID", all.y = T) %>%
  dplyr::select(-Comment) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(Nest.Attempt == 1) %>%
  filter(Status == 0) %>%
  group_by(NestID) %>%
  summarize(LastOnNest = max(Date)) %>%
  mutate(MovebankStartDate = paste(year(LastOnNest),"0101000000000",sep ="" )) %>%
  mutate(MovebankEndDate = paste(gsub("\\D","", LastOnNest), "000000000", sep = ""))
nest1Downloads <- as.data.frame(merge(gpsnest1, lastonDates, by = "NestID")) %>%
  dplyr::select(-NestID, -Nest.Attempt, -LastOnNest)

# ### I REMOVED THIS BIT SO THAT I COULD USE NEST CENTERS, IF YOU WANT TO INCLUED THEM LATER PROBABLY COULD USE LAST KNOWN LOCATION IN A TRACK
# #For birds without nests, we want a dataframe that matches lastonDates to rbind with to get a final download DF
# nonestDownloads <- data.frame(BirdID = gpshens$BirdID[which(gpshens$BirdID %notin% gpsnest1$BirdID)],
#                               Year = gpshens$Year[which(gpshens$BirdID %notin% gpsnest1$BirdID)],
#            MovebankStartDate = paste(gpshens$Year[which(gpshens$BirdID %notin% gpsnest1$BirdID)],"0101000000000",sep ="" ),
#            MovebankEndDate = paste(gpshens$Year[which(gpshens$BirdID %notin% gpsnest1$BirdID)],"0801000000000",sep ="" ))
# 
# downloads.df <- rbind(nest1Downloads, nonestDownloads)
downloads.df <- nest1Downloads

### Need a single dataframe for the crawlWrap function
# loop through the birds, download from movebank and combine into database for crawl
rm("raw.move.df")
rm("BMV.df")
for(i in 1:nrow(downloads.df)){
  animalname <- as.character(downloads.df$BirdID[i])
  timestart <- downloads.df$MovebankStartDate[i]
  timeend <- downloads.df$MovebankEndDate[i]
  year <- downloads.df$Year[i]
  
  turkeygps <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", 
                                    login = login,
                                    animal = animalname,
                                    timestamp_start = timestart,
                                    timestamp_end = timeend)
  t_turkeygps <- spTransform(turkeygps, center=T)
  turk_dBBMM <- brownian.bridge.dyn(t_turkeygps, raster = 30, location.error = 17, margin = 5, window.size = 15, ext = 1.3)
  ind.BMV.df <- data.frame(ID = paste(animalname, year, sep = "_"),
                       Timestamp = turkeygps@data$timestamp,
                       BMV = getMotionVariance(turk_dBBMM))
  turkeygps@data$ID <- paste(animalname, year, sep = "_")
  
  if(exists("raw.move.df")){
    raw.move.df <- rbind(raw.move.df, as.data.frame(turkeygps@data))
    BMV.df <- rbind(BMV.df, ind.BMV.df)
  }else{
    raw.move.df <- as.data.frame(turkeygps@data)
    BMV.df <- ind.BMV.df
  }
}

### Cleanup final download dataframe
move.df.final <- raw.move.df %>% 
  dplyr::select(lat = location_lat, lon = location_long, ID, Timestamp = timestamp)

### Need to have coordinates in UTM not LL
# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(move.df.final[,2:1], proj4string=CRS(projection(turkeygps)))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=19 ellps=WGS84"))

# add UTM locations to data frame
move.df.final$x <- attr(utmcoord,"coords")[,1]
move.df.final$y <- attr(utmcoord,"coords")[,2]


### Duplicate roost locations for all night time hours
# includes random error based on gps location uncertainty (~12.7m)
roosts <- move.df.final %>%
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
    dplyr::select(lon, lat, Timestamp, ID, x, y)
  
  for(i in 1:hr_after_sun){
    working.df$Timestamp[i] <- working.df$Timestamp[i] - (i*60*60) #change time
  }
  for(i in 1:hr_before_sun){
    working.df$Timestamp[i+hr_after_sun] <- working.df$Timestamp[i+hr_after_sun] + (i*60*60) #change time
  }
  move.df.final <- rbind(move.df.final, working.df) #add to original database
}

move.df.final <- move.df.final %>% arrange(ID, Timestamp) %>%
  mutate(YDay = yday(Timestamp))

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
crawllist <- unique(move.df.final$ID)
rm("wintertonest1_prepped")
for(i in 1:length(crawllist)){
  crawlid <- crawllist[i]
  bird1 <- move.df.final %>% filter(ID == crawlid)
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
  bird1prepped$YDay <- yday(bird1prepped$Timestamp)
  
  if(exists("wintertonest1_prepped")){
    wintertonest1_prepped <- rbind(wintertonest1_prepped, bird1prepped)
  }else{
    wintertonest1_prepped <- bird1prepped
  }
}

wintertonest1_cleaned <- wintertonest1_prepped %>%
  dplyr::select(ID, x, y, step, angle, Timestamp, wintercenter.dist, wintercenter.angle, nest1.dist, nest1.angle, YDay)

#check output
summary(wintertonest1_cleaned)

save(wintertonest1_cleaned, file = "hmminput.RData")
#################################################################################
### For Hierarchical Data, need to create levels within the data
## We want to average BMV by day

bmvmean.df <- BMV.df %>%
  mutate(Date = as.Date(Timestamp)) %>%
  group_by(ID, Date) %>%
  summarize(Timestamp2 = first(Timestamp),
            BMV_mean = mean(BMV))

wintertonestTS <- wintertonest1_cleaned %>%
  mutate(Date = as.Date(Timestamp)) %>%
  group_by(ID, Date) %>%
  summarize(Timestamp = first(Timestamp))

bmvmean.full.df <- as.data.frame(merge(wintertonestTS, bmvmean.df, by = c("ID", "Date"), all.x = T) %>% dplyr::select(-Date, -Timestamp2) %>%
  group_by(ID) %>%
  mutate(BMV_mean = ifelse(Timestamp == last(Timestamp) & is.na(BMV_mean), tail(BMV_mean, n=2L)[1], BMV_mean),
         BMV_mean = ifelse(is.na(BMV_mean), first(na.omit(BMV_mean)), BMV_mean))
  )

level1.df <- merge(bmvmean.full.df, wintertonest1_cleaned, by = c("ID", "Timestamp"), all.x = T) %>% mutate(level = "1")
level2i.df <- level1.df %>% mutate(BMV_mean = NA, level = "2i")
level2.df <- wintertonest1_cleaned %>% mutate(BMV_mean = NA, level = "2")

turk.hhmm.data <- rbind(level1.df, level2i.df, level2.df) %>%
  arrange(ID, Timestamp) %>%
  mutate(hour = hour(Timestamp))

summary(turk.hhmm.data)

unique(turk.hhmm.data$level)

write.csv(turk.hhmm.data, "HHMM Turkey Data.csv", row.names = F)


#################################################################################
### Temp Seasonal Dates for each bird based visually on BMV data and Plots
require(ggplot2)
for(i in 1:length(unique(BMV.df$ID))){
  indID <- unique(BMV.df$ID)[i]
  
  ind.BMV <- BMV.df %>% filter(ID == indID)
  
  BMV.quant <- quantile(ind.BMV$BMV, seq(0,1,.20),na.rm = T)
  
  ind.points <- raw.move.df %>% filter(ID == indID) %>% rename(Timestamp = timestamp)
  
  ind.points <-  merge(ind.points, ind.BMV, by = c("Timestamp", "ID"), all.x = T) %>%
    dplyr::select(ID, Timestamp, location_lat, location_long, BMV) %>%
    mutate(location_lat2 = lag(location_lat), location_long2 = lag(location_long))
  
  ind.BMV.plot <- ggplot(ind.BMV, aes(y = BMV, x = Timestamp)) +
    geom_point(aes(color = BMV))  +
    scale_color_gradientn(colours = rainbow(6), values = BMV.quant/max(BMV.quant))
  
  ind.points.plot <- ggplot(ind.points, aes(x = location_long, y = location_lat)) +
    geom_segment(aes(x = location_long, y = location_lat,
                     xend = location_long2, yend = location_lat2, color = BMV)) +
    geom_point(aes(color = BMV)) +
    scale_color_gradientn(colours = rainbow(6), values = BMV.quant/max(BMV.quant)) +
    theme_classic() +
    labs(title = paste(indID, " - Movement Track with BMV colors"))
  ggsave
    
}

















