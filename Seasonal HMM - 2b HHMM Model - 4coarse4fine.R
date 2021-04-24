### Load Necessary packages
lapply(c('data.tree', 'momentuHMM', 'parallel'), require, character.only = T)

#Load data and prep it for momentuHMM package
raw.move.data <- read.csv("HHMM Turkey Data.csv", colClasses = c(level = "character")) %>%
  dplyr::select(-step, -angle)

turkData <- prepData(raw.move.data,
                    covNames=c("hour", "YDay"),
                    hierLevels=c("1","2i","2"),
                    coordLevel = "2")

#############################################
### Define the hierarchical HMM structure ###
#############################################
### define hierarchical HMM: states 1-3 = coarse state 1 (resident/foraging); states 4-6 = coarse state 2 (mobile/foraging); states 7-9 = coarse state 3 (travelling/migrating)
hierStates <- data.tree::Node$new("turkey HHMM states")
hierStates$AddChild("Winter")   # resident/foraging
hierStates$Winter$AddChild("W1", state=1)
hierStates$Winter$AddChild("W2", state=2)
hierStates$Winter$AddChild("W3", state=3)
hierStates$Winter$AddChild("W4", state=4)
hierStates$AddChild("Dispersal")  # mobile/foraging
hierStates$Dispersal$AddChild("D1", state=5)
hierStates$Dispersal$AddChild("D2", state=6)
hierStates$Dispersal$AddChild("D3", state=7)
hierStates$Dispersal$AddChild("D4", state=8)
hierStates$AddChild("PreNesting")    # travelling/migrating
hierStates$PreNesting$AddChild("P1", state=9)
hierStates$PreNesting$AddChild("P2", state=10)
hierStates$PreNesting$AddChild("P3", state=11)
hierStates$PreNesting$AddChild("P4", state=12)
hierStates$AddChild("Nesting")    # travelling/migrating
hierStates$Nesting$AddChild("N1", state=13)
hierStates$Nesting$AddChild("N2", state=14)
hierStates$Nesting$AddChild("N3", state=15)
hierStates$Nesting$AddChild("N4", state=16)

print(hierStates,"state")

nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))

###################################################
### Parameter Distributions and Starting Values ###
###################################################
# data stream distributions: level 1 = coarse level (BMV_mean="gamma"); level 2 = fine level (step="gamma", angle = "wrpcauchy)
hierDist <- data.tree::Node$new("turkey HHMM dist")
hierDist$AddChild("level1")
hierDist$level1$AddChild("BMV_mean", dist="gamma")
hierDist$AddChild("level2")
hierDist$level2$AddChild("step", dist="gamma")
hierDist$level2$AddChild("angle", dist="wrpcauchy")

print(hierDist,"dist")

### defining start values based on those reported by Adam et al
BMV.mu0 <- c(75, 500, 300,1)
BMV.sigma0 <- c(50, 100, 100,50)
step.mu0 <- step.sigma0 <- step.zero <- angle.con <- list()
step.mu0[[1]] <- step.mu0[[2]] <- step.mu0[[3]] <- step.mu0[[4]] <-c(5, 20, 150, 300)
# step.mu0[[2]] <- c(5, 20, 150) #For if you want to specify different start values
# step.mu0[[3]] <- c(5, 20, 150)
# step.mu0[[4]] <- c(5, 20, 150)
step.sigma0[[1]] <- step.sigma0[[2]] <- step.sigma0[[3]] <- step.sigma0[[4]] <- c(20, 40, 100, 200)
# step.sigma0[[2]] <- c(0.043, 0.047, 0.342)
# step.sigma0[[3]] <- c(0.109, 0.462, 1.878)
# step.sigma0[[4]] <- c(0.109, 0.462, 1.878)
step.zero[[1]] <- step.zero[[2]] <- step.zero[[3]] <- step.zero[[4]] <- c(0.985, 0.006, 0.005, 0.004)
# step.zero[[2]] <- c(0.043, 0.047, 0.342)
# step.zero[[3]] <- c(0.109, 0.462, 1.878)
# step.zero[[4]] <- c(0.109, 0.462, 1.878)
angle.con[[1]] <- angle.con[[2]] <- angle.con[[3]] <- angle.con[[4]] <- c(0.1, 0.2, 0.3, 0.4)
# angle.con[[2]] <- c(1.791e-08, 0.035, 0.002)
# angle.con[[3]] <- c(0.012, 2.933e-04, 3.462e-09)
# angle.con[[4]] <- c(0.012, 2.933e-04, 3.462e-09)

Par0 <- list(BMV_mean=c(rep(BMV.mu0,each=4),rep(BMV.sigma0,each=4)),
             step=c(unlist(step.mu0),unlist(step.sigma0),unlist(step.zero)),
             angle=c(unlist(angle.con)))

####################################
### Design Matrix Specifications ###
####################################
BMV.DM <- matrix(c(rep(c(1,1,0,0,0,0,0,0),4),
                   rep(c(1,1,1,1,0,0,0,0),4),
                   rep(c(1,1,1,0,0,0,0,0),4),
                   rep(c(1,0,0,0,0,0,0,0),4),
                   rep(c(0,0,0,0,1,0,0,0),4),
                   rep(c(0,0,0,0,0,1,0,0),4),
                   rep(c(0,0,0,0,0,0,1,0),4),
                   rep(c(0,0,0,0,0,0,0,1),4)),
                 nrow = 2*nbStates, byrow = T,
                 ncol = 8,
                 dimnames = list(paste0(rep(c("mean_","sd_"),each=nbStates),1:nbStates),
                                 c(paste0(rep(c("mean_","sd_"),each=4),1:length(hierStates$children)))))

step.DM <- matrix(c(rep(c(1,0,0,0,0,0,0,0,0,0,0,0,
                          1,1,0,0,0,0,0,0,0,0,0,0,
                          1,1,1,0,0,0,0,0,0,0,0,0,
                          1,1,1,1,0,0,0,0,0,0,0,0),4),
                    rep(c(0,0,0,0,1,0,0,0,0,0,0,0,
                          0,0,0,0,0,1,0,0,0,0,0,0,
                          0,0,0,0,0,0,1,0,0,0,0,0,
                          0,0,0,0,0,0,0,1,0,0,0,0),4),
                    rep(c(0,0,0,0,0,0,0,0,1,1,1,1,
                          0,0,0,0,0,0,0,0,1,1,1,0,
                          0,0,0,0,0,0,0,0,1,1,0,0,
                          0,0,0,0,0,0,0,0,1,0,0,0),4)),
                 nrow = 3*nbStates, byrow = T,
                 ncol = 12,
                 dimnames = list(paste0(rep(c("mean_","sd_","zero_"),each=nbStates),1:nbStates),
                                 c(paste0(rep(c("mean_","sd_","zero_"),each=4),1:4))))

angle.DM <- matrix(c(rep(c(1,0,0,0,
                           1,1,0,0,
                           1,1,1,0,
                           1,1,1,1),4)),
                  nrow = 1*nbStates, byrow = T,
                  ncol = 4,
                  dimnames = list(paste0(rep(c("con_"),each=nbStates),1:nbStates),
                                  c(paste0(rep(c("con_"),each=4),1:4))))

  
DM <- list(BMV_mean = BMV.DM,
           step = step.DM,
           angle = angle.DM)


#############################################
### Specify Work Bounds for distributions ###
#############################################
#These correspond to the columns of the DM
#define the directions of the differences
BMVworkBounds <- matrix(c(-Inf,Inf, 
                           0,Inf,
                           0,Inf,
                           0, Inf,
                           -Inf, Inf,
                           -Inf, Inf,
                           -Inf, Inf,
                           -Inf, Inf), 
                         nrow = ncol(BMV.DM), byrow = T,
                         dimnames = list(colnames(BMV.DM), c("lower", "upper")))

stepworkBounds <- matrix(c(-Inf,Inf,
                           0,Inf,
                           0,Inf,
                           0, Inf,
                           -Inf, Inf,
                           -Inf, Inf,
                           -Inf, Inf,
                           -Inf, Inf,
                           -Inf, Inf,
                           0, Inf,
                           0, Inf,
                           0, Inf), 
                     nrow = ncol(step.DM), byrow = T,
                     dimnames = list(colnames(step.DM), c("lower", "upper")))

angleworkBounds <- matrix(c(-Inf,Inf,
                            -Inf,0,
                            -Inf,0,
                            -Inf,0),
                      nrow = ncol(angle.DM),
                      byrow=TRUE, dimnames=list(colnames(angle.DM),
                                                c("lower","upper"))) 

workBounds<-list(BMV_mean = BMVworkBounds,
                 step=stepworkBounds,
                 angle=angleworkBounds)

#####################################################
### Specify User Defined Bounds for distributions ###
#####################################################
#These correspond to the actual parameter values for each row
BMVBounds <- matrix(c(rep(c(0,Inf),4),
                      rep(c(0,Inf),4),
                      rep(c(0,Inf),4),
                      rep(c(0,Inf),4),
                      rep(c(0,Inf),4),
                      rep(c(0,Inf),4),
                      rep(c(0,Inf),4),
                      rep(c(0,Inf),4)), 
                     nrow = 4*ncol(BMV.DM), byrow = T,
                     dimnames = list(rep(colnames(BMV.DM),each=4), c("lower", "upper")))

stepBounds <- matrix(c(rep(c(0,Inf,
                             0,Inf,
                             0,Inf,
                             0,Inf),4),
                       rep(c(0,Inf,
                             0,Inf,
                             0,Inf,
                             0,Inf),4),
                       rep(c(0,.97,
                             0,.03,
                             0,.03,
                             0,.03),4)), 
                     nrow = 4*ncol(step.DM), byrow = T,
                     dimnames = list(rep(colnames(step.DM),each=4), c("lower", "upper")))

angleBounds <- matrix(c(rep(c(0,1,
                              0,1,
                              0,1,
                              0,1),4)),
                      nrow = 4*ncol(angle.DM), byrow=TRUE, 
                      dimnames=list(rep(colnames(angle.DM),each=4), c("lower","upper"))) 

userBounds <- list(BMV_mean = BMVBounds,
                   step = stepBounds,
                   angle = angleBounds)

################################
### Initial parameter Values ###
################################
# get initial parameter values for data stream probability distributions on the working scale
Par <- getParDM(turkData,
                hierStates=hierStates,
                hierDist=hierDist,
                Par=Par0,
                DM=DM)

##################################################
### Define Hierarchical Transition Probability ###
##################################################
hierFormula <- data.tree::Node$new("turkey HHMM formula")
hierFormula$AddChild("level1", formula=~YDay)
hierFormula$AddChild("level2", formula=~cosinor(hour, period=24))


#number of repeats of the betas will depend on the number of covariates in the above formulas, repeats will be level specific
hierBeta <- data.tree::Node$new("turkey beta")
hierBeta$AddChild("level1",beta=matrix(c(0,0,-100,
                                         -100,0,-100,
                                         -100,-100,0,
                                         -100,-100,-100),2,length(hierStates$children)*(length(hierStates$children)-1), byrow = T))
hierBeta$AddChild("level2")
hierBeta$level2$AddChild("Winter",beta=matrix(rep(c(0,0,0,
                                                    0,0,0,
                                                    0,0,0,
                                                    0,0,0),3),3,length(hierStates$Winter$children)*(length(hierStates$Winter$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("Dispersal",beta=matrix(rep(c(0,0,0,
                                                       0,0,0,
                                                       0,0,0,
                                                       0,0,0),3),3,length(hierStates$Dispersal$children)*(length(hierStates$Dispersal$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("PreNesting",  beta=matrix(rep(c(0,0,0,
                                                          0,0,0,
                                                          0,0,0,
                                                          0,0,0),3),3,length(hierStates$PreNesting$children)*(length(hierStates$PreNesting$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("Nesting",  beta=matrix(rep(c(0,0,0,
                                                       0,0,0,
                                                       0,0,0,
                                                       0,0,0),3),3,length(hierStates$Nesting$children)*(length(hierStates$Nesting$children)-1),byrow=TRUE))

#Fix the betas for the coarse scale to prevent switching to previous states
hierFixBeta <- data.tree::Node$new("turkey beta")
hierFixBeta$AddChild("level1",beta=matrix(c(NA,NA,-100,
                                         -100,NA,-100,
                                         -100,-100,NA,
                                         -100,-100,-100),2, length(hierStates$children)*(length(hierStates$children)-1), byrow = T))
hierFixBeta$AddChild("level2")
hierFixBeta$level2$AddChild("Winter",beta=matrix(rep(c(NA,NA,NA,
                                                    NA,NA,NA,
                                                    NA,NA,NA,
                                                    NA,NA,NA),3),3,length(hierStates$Winter$children)*(length(hierStates$Winter$children)-1),byrow=TRUE))
hierFixBeta$level2$AddChild("Dispersal",beta=matrix(rep(c(NA,NA,NA,
                                                       NA,NA,NA,
                                                       NA,NA,NA,
                                                       NA,NA,NA),3),3,length(hierStates$Dispersal$children)*(length(hierStates$Dispersal$children)-1),byrow=TRUE))
hierFixBeta$level2$AddChild("PreNesting",  beta=matrix(rep(c(NA,NA,NA,
                                                          NA,NA,NA,
                                                          NA,NA,NA,
                                                          NA,NA,NA),3),3,length(hierStates$PreNesting$children)*(length(hierStates$PreNesting$children)-1),byrow=TRUE))
hierFixBeta$level2$AddChild("Nesting",  beta=matrix(rep(c(NA,NA,NA,
                                                       NA,NA,NA,
                                                       NA,NA,NA,
                                                       NA,NA,NA),3),3,length(hierStates$Nesting$children)*(length(hierStates$Nesting$children)-1),byrow=TRUE))

###################################
### Define Initial Distribution ###
###################################
hierDelta <- data.tree::Node$new("turkey HHMM delta")
hierDelta$AddChild("level1",delta=matrix(c(-100, -100, -100),1))
hierDelta$AddChild("level2")
hierDelta$level2$AddChild("Winter",delta=matrix(c(0,0,0),1))
hierDelta$level2$AddChild("Dispersal",delta=matrix(c(0,0,0),1))
hierDelta$level2$AddChild("PreNesting",delta=matrix(c(0,0,0),1))
hierDelta$level2$AddChild("Nesting",delta=matrix(c(0,0,0),1))

#Fix the initial distribution to start in Winter
hierFixDelta <- data.tree::Node$new("turkey HHMM delta")
hierFixDelta$AddChild("level1",delta=matrix(c(-100, -100, -100),1))
hierFixDelta$AddChild("level2")
hierFixDelta$level2$AddChild("Winter",delta=matrix(c(NA,NA,NA),1))
hierFixDelta$level2$AddChild("Dispersal",delta=matrix(c(NA,NA,NA),1))
hierFixDelta$level2$AddChild("PreNesting",delta=matrix(c(NA,NA,NA),1))
hierFixDelta$level2$AddChild("Nesting",delta=matrix(c(NA,NA,NA),1))


##################################
### Check Model Specifications ###
##################################
# check hierarchical model specification and parameters
checkPar0(turkData,
          hierStates=hierStates,
          hierDist=hierDist,
          Par0=Par,
          hierFormula=hierFormula,
          DM=DM, 
          workBounds = workBounds, 
          userBounds = userBounds,
          hierBeta=hierBeta,
          hierDelta=hierDelta,
          fixPar=list(beta=hierFixBeta, delta=hierFixDelta))

#################################
### Run the Hierarchical HHMM ###
#################################
# compute hierarchical state transition probabilities based on initial values
iTrProbs <- getTrProbs(turkData,
                       hierStates=hierStates,
                       hierDist=hierDist,
                       Par0=Par,
                       hierFormula=hierFormula,
                       DM=DM, 
                       workBounds = workBounds, 
                       userBounds = userBounds,
                       hierBeta=hierBeta,
                       hierDelta=hierDelta,
                       fixPar=list(beta=hierFixBeta, delta=hierFixDelta)
                       )
iTrProbs$level1$gamma[,,1] # tpm at first time step for level1
lapply(iTrProbs$level2,function(x) x$gamma[,,1]) # tpm at first time step for level2

print(Sys.time())
turk.hhmm.4.4 <- fitHMM(codData,
                        hierStates=hierStates,
                        hierDist=hierDist,
                        hierFormula=hierFormula,
                        Par0=Par,
                        hierBeta=hierBeta,
                        hierDelta=hierDelta,
                        DM=DM, 
                        workBounds = workBounds, 
                        userBounds = userBounds,
                        fixPar=list(beta=hierFixBeta, delta=hierFixDelta),
                        ncores = detectCores()/2
                        )
turk.hhmm.4.4
print(Sys.time())













