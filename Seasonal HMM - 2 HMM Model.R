### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 0 # number attempt to re-fit based on random perturbation
nbStates <- 4 # Number of states
stateNames <- c("winter", "dispersal", "prenesting", "nesting") # label states
state_abb <- c("W", "D", "P", "N") #abbreviations for column names
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

###################################
### MODEL 1 - Basic Model
# Constrain transition probabilities between seasons

fixbeta <- matrix(c(NA, NA, -100, #Winter
                    -100, NA, -100, #Dispersal
                    -100, -100, NA, #Searching
                    -100, -100, -100)) #Nesting

## constrain step length parameters: 
# Mean -> W < D, N < P < D
# SD -> P < D
# Zero -> no relationship
stepDM <- matrix(c(
  1,1,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,1,1,1,0,0,0,0,0,
  0,0,0,0,1,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",state_abb),paste0("sd_",state_abb), paste0("zero_",state_abb)),
                c("mean_N:(Intercept)", "mean_W","mean_D", "mean_P",
                  "sd_W", "sd_P)","sd_D", "sd_N",
                  paste0("zero_",state_abb,":(Intercept)"))))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf,0,0,0,
                           rep(-Inf,8),
                           rep(Inf,12)),nrow = ncol(stepDM),
                         dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(0,Inf, #100,250,175,5 has given best results
                       0,Inf,
                       0,Inf,
                       0, Inf,
                       10, Inf,
                       10, Inf,
                       10, Inf,
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
                  1,0,1,0,
                  1,0,0,0),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration:(Intercept)","concentration_W", "concentration_P", "concentration_D")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,1, 
                        0,1, 
                        0,1, 
                        0,1),nrow = nbStates,
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
                       delta0 = delta0,
                       ncores = (detectCores()/2))


turk_m1s.wtn
states1 <- viterbi(turk_m1s.wtn)
table(states1)/nrow(wintertonest1_prepped)
plot(turk_m1s.wtn)
