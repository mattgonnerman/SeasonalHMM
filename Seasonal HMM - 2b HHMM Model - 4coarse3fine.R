### Load Necessary packages
lapply(c('data.tree', 'momentuHMM'), require, character.only = T)

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
hierStates$AddChild("Dispersal")  # mobile/foraging
hierStates$Dispersal$AddChild("D1", state=4)
hierStates$Dispersal$AddChild("D2", state=5)
hierStates$Dispersal$AddChild("D3", state=6)
hierStates$AddChild("PreNesting")    # travelling/migrating
hierStates$PreNesting$AddChild("P1", state=7)
hierStates$PreNesting$AddChild("P2", state=8)
hierStates$PreNesting$AddChild("P3", state=9)
hierStates$AddChild("Nesting")    # travelling/migrating
hierStates$Nesting$AddChild("N1", state=10)
hierStates$Nesting$AddChild("N2", state=11)
hierStates$Nesting$AddChild("N3", state=12)

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
step.mu0[[1]] <- step.mu0[[2]] <- step.mu0[[3]] <- step.mu0[[4]] <-c(5, 20, 150)
# step.mu0[[2]] <- c(5, 20, 150) #For if you want to specify different start values
# step.mu0[[3]] <- c(5, 20, 150)
# step.mu0[[4]] <- c(5, 20, 150)
step.sigma0[[1]] <- step.sigma0[[2]] <- step.sigma0[[3]] <- step.sigma0[[4]] <- c(20, 40, 200)
# step.sigma0[[2]] <- c(0.043, 0.047, 0.342)
# step.sigma0[[3]] <- c(0.109, 0.462, 1.878)
# step.sigma0[[4]] <- c(0.109, 0.462, 1.878)
step.zero[[1]] <- step.zero[[2]] <- step.zero[[3]] <- step.zero[[4]] <- c(0.99, 0.006, 0.004)
# step.zero[[2]] <- c(0.043, 0.047, 0.342)
# step.zero[[3]] <- c(0.109, 0.462, 1.878)
# step.zero[[4]] <- c(0.109, 0.462, 1.878)
angle.con[[1]] <- angle.con[[2]] <- angle.con[[3]] <- angle.con[[4]] <- c(0.1, 0.2, 0.3)
# angle.con[[2]] <- c(1.791e-08, 0.035, 0.002)
# angle.con[[3]] <- c(0.012, 2.933e-04, 3.462e-09)
# angle.con[[4]] <- c(0.012, 2.933e-04, 3.462e-09)

Par0 <- list(BMV_mean=c(rep(BMV.mu0,each=3),rep(BMV.sigma0,each=3)),
             step=c(unlist(step.mu0),unlist(step.sigma0),unlist(step.zero)),
             angle=c(unlist(angle.con)))

####################################
### Design Matrix Specifications ###
####################################
BMV.DM <- matrix(c(rep(c(1,1,0,0,0,0,0,0),3),
                   rep(c(1,1,1,1,0,0,0,0),3),
                   rep(c(1,1,1,0,0,0,0,0),3),
                   rep(c(1,0,0,0,0,0,0,0),3),
                   rep(c(0,0,0,0,1,0,0,0),3),
                   rep(c(0,0,0,0,0,1,0,0),3),
                   rep(c(0,0,0,0,0,0,1,0),3),
                   rep(c(0,0,0,0,0,0,0,1),3)),
                 nrow = 2*nbStates, byrow = T,
                 ncol = 8,
                 dimnames = list(paste0(rep(c("mean_","sd_"),each=nbStates),1:nbStates),
                                 c(paste0(rep(c("mean_","sd_"),each=4),1:length(hierStates$children)))))

step.DM <- matrix(c(rep(c(1,0,0,0,0,0,0,0,0,
                          1,1,0,0,0,0,0,0,0,
                          1,1,1,0,0,0,0,0,0),4),
                    rep(c(0,0,0,1,0,0,0,0,0,
                          0,0,0,0,1,0,0,0,0,
                          0,0,0,0,0,1,0,0,0),4),
                    rep(c(0,0,0,0,0,0,1,1,1,
                          0,0,0,0,0,0,1,1,0,
                          0,0,0,0,0,0,1,0,0),4)),
                 nrow = 3*nbStates, byrow = T,
                 ncol = 9,
                 dimnames = list(paste0(rep(c("mean_","sd_","zero_"),each=nbStates),1:nbStates),
                                 c(paste0(rep(c("mean_","sd_","zero_"),each=3),1:3))))

angle.DM <- matrix(c(rep(c(1,0,0,
                           1,1,0,
                           1,1,1),4)),
                  nrow = 1*nbStates, byrow = T,
                  ncol = 3,
                  dimnames = list(paste0(rep(c("con_"),each=nbStates),1:nbStates),
                                  c(paste0(rep(c("con_"),each=3),1:3))))

  
DM <- list(BMV_mean = BMV.DM,
           step = step.DM,
           angle = angle.DM)

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



hierBeta <- data.tree::Node$new("turkey beta")
hierBeta$AddChild("level1",beta=matrix(c(-18.585, -2.86, -2.551, -1.641, -2.169, -2.415),1,length(hierStates$children)*(length(hierStates$children)-1)))
hierBeta$AddChild("level2")
hierBeta$level2$AddChild("Winter",beta=matrix(c(rep(c(, 2,  2.765, -1.607,  2.273,  4.842),4),4,length(hierStates$resForage$children)*(length(hierStates$resForage$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("Dispersal",beta=matrix(c(-2.156, -3.662,   3.01,  0.597, -0.313,  2.897, 
                                                   0.067,  -1.22, -0.799, -0.797,   0.15,  0.379, 
                                                   -0.112, -0.195, -0.269, -0.215,  1.539,  0.728),4,length(hierStates$mobForage$children)*(length(hierStates$mobForage$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("PreNesting",  beta=matrix(c( -2.53, -4.279,  2.507, -0.228, 10.803, 12.873, 
                                                    -0.04,  1.221, -0.301,  0.284, -0.106, -0.077, 
                                                    0.629, -0.226, -0.253, -0.303,  0.011,  0.036),4,length(hierStates$transit$children)*(length(hierStates$transit$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("Nesting",  beta=matrix(c( -2.53, -4.279,  2.507, -0.228, 10.803, 12.873, 
                                                       -0.04,  1.221, -0.301,  0.284, -0.106, -0.077, 
                                                       0.629, -0.226, -0.253, -0.303,  0.011,  0.036),4,length(hierStates$transit$children)*(length(hierStates$transit$children)-1),byrow=TRUE))


























