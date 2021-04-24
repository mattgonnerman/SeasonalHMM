### Load Necessary packages
lapply(c('dplyr', 'move', 'momentuHMM', 'lubridate', 'rgdal'), require, character.only = T)

### Login Credentials for MoveBanks
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

### Download all Turkey Movement Data
raw.move.all <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login)

### Filter Data down to what will be used in analysis
# Group Individuals X Year
# 
move.all <- as.data.frame(coordinates(raw.move.all)) %>%
  mutate(Timestamp = timestamps(indbird_mbdown)) %>%
  mutate(Timestamp = round_date(Timestamp, unit = "hour")) %>%
  mutate(ID = paste(animalname, yearnest, sep = "_"))
  