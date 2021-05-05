### Load Necessary packages
lapply(c('data.tree', 'momentuHMM', 'parallel', 'dplyr'), require, character.only = T)

#Load data and prep it for momentuHMM package
raw.move.data <- read.csv("HHMM Turkey Data.csv", colClasses = c(level = "character")) %>%
  dplyr::select(-step, -angle)

raw.move.data.rough <- raw.move.data %>% filter(level == 1)

turkData <- prepData(raw.move.data.rough,
                     covNames=c("hour", 
                                "YDay", 
                                "wintercenter.dist", 
                                "wintercenter.angle", 
                                "nest1.dist", 
                                "nest1.angle")
                     )

### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("winter", "dispersal", "prenesting", "nesting") # label states
state_abb <- c("W", "D", "P", "N") #abbreviations for column names
dist = list(BMV_mean = "gamma") # distributions for observation processes

###################################
### MODEL 1 - Basic Model
## constrain step length parameters: 
# Mean -> W < D, N < P < D
# SD -> P < D
# Zero -> no relationship
BMVDM <- matrix(c(
  1,1,0,0,0,0,0,0,
  1,1,1,1,0,0,0,0,
  1,1,1,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,1),
  nrow = 2*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",state_abb),paste0("sd_",state_abb)),
                c("mean_N:(Intercept)", "mean_W","mean_P", "mean_D",
                  "sd_W", "sd_P","sd_D", "sd_N")
                )
  )

#define the directions of the differences
BMVworkBounds <- matrix(c(-Inf,0,0,0,
                           rep(-Inf,4),
                           rep(Inf,8)),nrow = ncol(BMVDM),
                         dimnames=list(colnames(BMVDM),c("lower","upper")))

#Bundle individual parameter DM and workbounds
DM<-list(BMV_mean=BMVDM)
workBounds<-list(BMV_mean=BMVworkBounds)


### TPM and Initial Distribution
fixbeta <- matrix(c(NA, NA, -100, #Winter
                    -100, NA, -100, #Dispersal
                    -100, -100, NA, #Searching
                    -100, -100, -100),#Nesting
                  4,12, byrow = T) 
# beta0 <- matrix(c(-.1,-.1,-100,
#                   -100,-.1,-100,
#                   -100,-100,-.1,
#                   -100,-100,-100),
#                 4,12, byrow = T)
beta0 <- matrix(c(-.1,-.1,-100,-100,-.1,-100,-100,-100,-.1,-100,-100,-100,
                  -.1,-.1,-100,-100,-.1,-100,-100,-100,-.1,-100,-100,-100,
                  -.1,-1,-100,-100,-.1,-100,-100,-100,-.1,-100,-100,-100,
                  -.1,-1,-100,-100,-.1,-100,-100,-100,-.1,-100,-100,-100),
                4,12, byrow = T)

fixdelta <- matrix(c(.99999997, .00000001, .00000001, .00000001))
delta0 <- matrix(c(.99997, .00001, .00001, .00001),1,
                 dimnames = list(NULL,c("winter", "dispersal","prenesting", "nesting")))

# Formula for transition probabilities
# transFormula <- formula("~betaCol1(I(YDay > 90)) + betaCol2(I(YDay > 90)) + betaCol1(I(wintercenter.dist>5000)) + betaCol2(I(nest1.dist<5000))")
transFormula <- formula("~YDay + wintercenter.dist + nest1.dist")
# transFormula <- formula("~YDay")
# initial parameters
Par <- list(BMV_mean=c(100,400,300,10, 
                   100,250,500,30)
            )
Par0_m1s.wtn <- getParDM(data = turkData,
                         nbStates = nbStates, 
                         dist = dist,
                         Par = Par, 
                         DM = DM,  
                         workBounds = workBounds,
                         delta0 = delta0,
                         beta0 = beta0,
                         estAngleMean = list(angle = FALSE),
                         fixPar=list(beta=fixbeta, delta=fixdelta),
                         formula = transFormula) #Try changing this to TRUE

checkPar0(turkData,
          nbStates = nbStates, 
          dist = dist, 
          Par0 = Par0_m1s.wtn,
          DM=DM, 
          delta0 = delta0,
          beta0 = beta0,
          workBounds = workBounds,
          fixPar=list(beta=fixbeta, delta=fixdelta),
          formula = transFormula
)



# fit model 1
turk_m1s.wtn <- fitHMM(data = turkData,
                       nbStates = nbStates, 
                       dist = dist, 
                       Par0 = Par0_m1s.wtn, 
                       DM = DM, 
                       workBounds = workBounds, 
                       estAngleMean = list(angle=FALSE), #Try changing this to TRUE
                       prior = prior,
                       stateNames = stateNames,
                       fixPar=list(beta=fixbeta, delta=fixdelta),
                       formula = transFormula,
                       retryFits = retryFits,
                       ncores = (detectCores()/2))


turk_m1s.wtn
states1 <- viterbi(turk_m1s.wtn)
table(states1)/nrow(turkData)
plot(turk_m1s.wtn)

turkData.states <- turkData
turkData.states$State <- states1
