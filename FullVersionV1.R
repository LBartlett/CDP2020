#### Simulations for anticipating when an observed cohort of developing honey bee brood will perform certain tasks in future
#### Intended to guide experiments and aid in analysis of experimental data where brood are subjected to acute treatment events
#### Associated Manuscript:
#### Written by:
#### Lewis J Bartlett lewis.bartlett@uga.edu

## Here, uses an example involving brood transplantation 


#######################
##### User Inputs #####
#######################


### Set directory where reference data, parameter sheet, assessments etc. are stored - probably the wd for this pulled repo

RefDir <- getwd()

## Set a storage directory for simulations outputs, better not to be on 'cloud storage' as can cause problems
OpDir <- 'D:/BroodMixAnalaysis/SimOutputs/'

## Set location from which to base development time
# Create own by measuring bee development times
# Or choose from those provided as part of this publication:
# 'Alabama', 'Georgia', ...

SITE <- 'Alabama'

## Choose task allocation reference
# Options: 'Delaplane' [recommended]; 'Huang', ...

TaskRef<- 'Delaplane'

## All other parameters should be specified in 'Initialisation.csv'


#######################
#######  Main  ########
#######################


# Set parameters to be explored / varied / set

# Read in initialisation file where parameters can be edited
Init <- read.csv(paste0(RefDir,'/Initialisation.csv'),header=T)

#Ages
#Select location to draw development times from
DevLoc <- SITE

#Pull development data for that location

DevData <- read.csv(file = paste0(RefDir,"/DevelopmentData/",DevLoc,".csv"),
                    header = T)


# Set sample density of development time distribution
# Set repeat number (the above can probably better increase repeats, to a degree)

ANum <- Init$Value[which(Init$Parameter == 'DTSamp')]
RepT <- Init$Value[which(Init$Parameter == 'Reps')]

# Set artificially extreme min/mix for development time distribution generation sequence
XDist <- seq(15,30, length.out = 80)

# Create a normal distribution based off location's development time data
DevDist <- dnorm(XDist, mean = mean(DevData$EDays), sd = sd(DevData$EDays))

ERange <- round(sample(x = XDist, size = ANum, replace = T, prob = DevDist))

# Death rates per day
# Vary based on literature if known toxins or pathogens in colony
# or if large number of inviable eggs (diploid males)
#Egg
EMort <- Init$Value[which(Init$Parameter == 'EggMort')]
#Larvae
LMort <- Init$Value[which(Init$Parameter == 'LarvMort')]
#Pupa
PMort <- Init$Value[which(Init$Parameter == 'PupMort')]
#Adult base rate plus additional per day of age
AMortBase <- Init$Value[which(Init$Parameter == 'BaseAdultMort')]
AMortScale <- Init$Value[which(Init$Parameter == 'DailyAddMort')]

#Maximum number of cells on a brood frame
# 5500 recommended based on a langstroth deep
MaxSize <- Init$Value[which(Init$Parameter == 'FrameSize')]

#How long to run each repeat for
RunT <- Init$Value[which(Init$Parameter == 'RunTime')]

# Get list of frame assessment sheets
# Separate assessment sheet for each assessor
# Matches/relies on using location [SITE] in the file name

FAsses <- list.files(path = paste0(RefDir,"/FrameAssessments/"), pattern = SITE,
                     recursive = F, full.names = T)

# Begin individual replicates

RepTrack <-0

for (EmA in ERange){
  for(Rep in 1:RepT){
    
    #Set pupation and hatch times based on emergence time
    PuA <- round(EmA*0.5)
    HaA <- round(EmA*0.15)
    
    # Select a random initialisation data sheet from those provided
    
    FrameD <- read.csv(file = FAsses[sample(x = 1:length(FAsses), size = 1)], 
                       header = T)
    
    
    # Spin up main data frame with variably assigned bee ages guided by brood data
    # Initialises model
    
    #Pull reference
    ExpRef <- FrameD
    ExpRef$CTot <- MaxSize*FrameD$Size
    
    # Matrix with row for each frame, columns for numbers of individuals at each age from laying
    AMat <- matrix(0, nrow = NROW(ExpRef), ncol = (EmA+RunT))
    
    # Assign initial ages according to classifications in read-in file. Uses a uniform random number and then scales
    # Along with some moderately efficient R matrix magic
    
    # Random numbers
    AMat[,1:EmA] <- runif(length(AMat[,1:EmA]), min = 0, max = 1)
    
    # Eggs
    AMat[,1:HaA] <- round( AMat[,1:HaA] * (matrix(rep(x = (FrameD$Eggs*ExpRef$CTot/rowSums(AMat[,1:HaA])), times = HaA), nrow = NROW(AMat), ncol = HaA)) )
    
    #Larvae
    AMat[,(1+HaA):PuA] <- round( AMat[,(1+HaA):PuA] * (matrix(rep(x = (FrameD$Larvae*ExpRef$CTot/rowSums(AMat[,(1+HaA):PuA])), times = PuA - HaA), nrow = NROW(AMat), ncol = PuA - HaA)) )
    
    #Pupae
    AMat[,(1+PuA):EmA] <- round( AMat[,(1+PuA):EmA] * (matrix(rep(x = (FrameD$Capped*ExpRef$CTot/rowSums(AMat[,(1+PuA):EmA])), times = EmA - PuA), nrow = NROW(AMat), ncol = EmA - PuA)) )
    
    # Total number of bees on each frame according to this simulation at time of assessment
    ExpRef$BTot <- rowSums(AMat)
    
    # Create death rate vector
    # Rows correspond to bee age, with appropriate death likelihood for that day
    AD <- (1:(NCOL(AMat)-EmA))
    MortV <- c(rep(EMort, times = HaA), 
               rep(LMort, times = PuA - HaA), 
               rep(PMort, times = EmA - PuA), 
               1/(1+ ( ((1/AMortBase)-1) * (exp(-AMortScale*AD^2)  ) ) )
    )
    
    # Convert to survival vector for speedier code
    SurV <- 1-MortV
    
    # Create tracking array (starts at 'day 0' i.e. initialisation, day of assessment)
    PolyTrack <- array(data = NA, dim = c(NROW(AMat), NCOL(AMat), RunT+1))
    
    # Store day '0' (1)
    
    PolyTrack[,,1] <- AMat
    
    # Run for number of days specified
    for(D in 1:RunT){
      
      # Calculate survivors
      SMat <- matrix(mapply(rbinom, n = c(1), size = as.vector(t(AMat)), prob = SurV), ncol = NCOL(AMat), byrow = T)
      
      # Replace AMat with survivors, but shift by 1 day
      
      # Trim off last column of survivor matrix
      SMat <- SMat[,-NCOL(SMat)]
      
      # Replace all but 'youngest' column of AMat
      AMat[,2:NCOL(AMat)] <- SMat
      
      # Replace youngest column of AMat with 0s
      AMat[,1] <- 0
      
      # Update Tracker
      PolyTrack[,,D+1] <- AMat
      
    }
    
    # Individual run finished, save.
    RepTrack <- RepTrack+1
    save(ExpRef, PolyTrack, EmA, PuA, HaA,
         file = paste0(OpDir,'Run',as.character(RepTrack),'.RData')
    )   
    
    # Monitor progress
    print(RepTrack)
    
  }
}

## Clean out all the simulation workspace. All data now 'hard saved' in the above directory
rm(list = 
     ls()[- match(c('OpDir', 'RefDir', 'Init','TaskRef'), ls())]
)

###############################
### Part Two: task allocation##
###############################

# Decide what resolutions plots are to be at
# options are 'Colony',''Experiment' ,or 'Both'
TRes <- 'Colony'

#### For approximating age groupings & survivorship of bee brood cohorts

# Set Directory to read in runs from
RunsDir <- OpDir

# Get number and filepaths of runs
TRuns <- list.files(path = RunsDir, pattern = '.RData', full.names = T)

NRuns <- NROW(TRuns)

# Pull runtime of sims
RT <- Init$Value[which(Init$Parameter == 'RunTime')]

# Create lists to hold reference data frames
# tracking arrays, and changing variables from each run

RefL <- vector('list', NRuns)
TrackL <- vector('list', NRuns)
EmL <- vector('list', NRuns)
PuL <- vector('list', NRuns)
HaL <- vector('list', NRuns)

# Read in all data

for(Run in 1:NRuns){
  
  load(TRuns[Run])
  RefL[[Run]] <- ExpRef
  TrackL[[Run]] <- PolyTrack
  
  EmL[[Run]] <- EmA
  PuL[[Run]] <- PuA
  HaL[[Run]] <- HaA
  
}

# Clear out the no longer required parts of the workspace
rm(list = setdiff(ls(), c('RefL', 'TrackL', 'EmL', 'PuL', 'HaL', 'NRuns', 'RT', 'RefDir','RunsDir', 'TaskRef', 'TRes')))

## Create a function to collapse totals of bees across frames into recipient colonies

CollapseF <- function(YardA, RT, EmA, YardRef){
  
  Hives <- unique(cbind(YardRef$Colony, YardRef$Treatment))
  HivesRef <- cbind(YardRef$Colony, YardRef$Treatment)
  
  # Create collapse array
  
  OutArr <- array(NA, dim = c(NROW(Hives),RT+EmA,RT+1))
  
  for(N in 1:NROW(Hives)){
    
    OutArr[N,,] <- apply(X = YardA[which(HivesRef[,1]==Hives[N,1] & HivesRef[,2]== Hives[N,2]),,], 
                         MARGIN = c(2,3), FUN = sum)
    
  }
  
  return(OutArr)
  
}


## Create list of collapsed data

CollapseL <- vector('list', length = NRuns)

CollapseL <- mapply(FUN = CollapseF, YardA = TrackL, EmA = EmL, YardRef = RefL, RT = RT, SIMPLIFY = F)


# Now we have what was by-frame data collapsed down to by-hive

# We can now do 'whole colony' task partitioning

# Define age distributions of different tasks where 1 = day of eclosure (emergence)
# Read in median ages for different tasks

TaskAges <- read.csv(paste0(RefDir,'/TaskAges/',TaskRef,'.csv'), header = T)
TaskAges$Task <- as.character(TaskAges$Task)
TaskAges$Category <- as.character(TaskAges$Category)

# Create generic function to return a poisson density distribution for the task

TaskDP <- function(Task, TaskList, MaxAge){
  
  x <- 1:(MaxAge)
  
  DP <- dpois(x, lambda = TaskList$MedA[which(TaskList$Task == Task)])
  
  
}

## Create a function to add up task contribution for the target cohort

TaskAll <- function(Task, TaskList, MaxAge, EmA, TrackA){
  
  DistV <- rep(0, times = EmA)
  
  DistV <- append(DistV, TaskDP(Task = Task, TaskList = TaskList, MaxAge = MaxAge))
  
  TaskEmph <- apply(X = sweep(TrackA, MARGIN = 2, DistV, FUN = '*'),
                    MARGIN = 3, FUN = sum)
  
  return(TaskEmph)
}


# Create function to give a day for a specific quantile value for a vector of these trackers

QuantDay <- function(PTV, Quant){
  
  QCUM <- cumsum(PTV)/sum(PTV)
  
  QD <- which(
    abs(QCUM - Quant)
    ==
      min(abs(QCUM - Quant))
  )
  
  
  return(QD)
  
}

# Create function for average emphasis for each day across all replicates
# creates the 'solid' plot lines

AEmph <- function(TaskEmphases){
  
  Days <- 1:NROW(TaskEmphases[[1]])
  
  ME <- vector(length = NROW(Days))
  
  for (D in Days){
    TETrack <- 0
    for(X in 1:length(TaskEmphases)){
      
      TETrack <- TETrack + TaskEmphases[[X]][D]
      
    }
    
    ME[D] <- TETrack/(length(TaskEmphases))
    
  }
  
  return(ME)
  
}
# Create a function modified from above where the 'task' is simply existing as an adult bee

AdultAll <- function(MaxAge, EmA, TrackA){
  
  DistV <- rep(0, times = EmA)
  
  DistV <- append(DistV, rep(1, times = MaxAge))
  
  TaskEmph <- apply(X = sweep(TrackA, MARGIN = 2, DistV, FUN = '*'),
                    MARGIN = 3, FUN = sum)
  
  return(TaskEmph)
}

# Create function to subset the full data set to just the desired colonies
# This can be easily cloned and modified to create bespoke plotting / subsets based on e.g. treatments, if treatments were per-frame and split across colonies

CFilt <- function(TrackA, TrackRef, ColFil){
  
  Filt <- vector(mode = 'list', length = length(TrackA))
  
  for(N in 1:(length(TrackA))){
    
    Filt[[N]] <- TrackA[[N]][which(TrackRef[[1]]$Colony == ColFil),,]
    
  }
  
  return(Filt)
  
}


# Here is where we resolve to colony and/or experiment if we want by colony or across whole experiment graphical plots.

# Set some universal plotting parameters
# Colour w/ opacity
CC <- rgb(0.6, 0.6, 0.6, 0.6)
# Suitable margins
par(mar=c(5,6,8,2))

if(TRes == 'Colony'  |  TRes == 'Both'){
  
  # Pull number of colonies
  Colonies <- unique(RefL[[1]]$Colony)
  NC <- NROW(Colonies)
  #NC <- NROW(CollapseL[[1]][,1,1])
  
  # Create storage vectors for each colony
  ColonyTB <- vector(mode = "list", length = NC)
  ColonyAB <- vector(mode = "list", length = NC)
  
  # Cycle through colonies
  
  for(TC in 1:NC){
    
    # Colony name
    TCName <- Colonies[TC]
    
    ## Create list to hold task emphases for cohort for each task
    TaskEs <- vector(mode = "list", length = NROW(TaskAges))
    
    ## Subset full tracking list to only data from this colony [track list colony subset]
    
    TLCS <- CFilt(TrackA = TrackL, TrackRef = RefL, ColFil = TCName)
    
    for(N in 1:NROW(TaskAges)){
      
      CurTask <- TaskAges$Task[N]
      
      # Create list of relative contribution of target cohort to this task across experiment through time
      
      TaskAL <- mapply(FUN = TaskAll, TrackA = TLCS, EmA = EmL, 
                       MoreArgs = list( Task = CurTask, TaskList = TaskAges, MaxAge = RT),
                       SIMPLIFY = F)
      
      PeakDays <- lapply(TaskAL, QuantDay, Quant = 0.5)
      
      PDM <- round(mean(unlist(PeakDays)), digits = 0)
      
      Q16Days <- lapply(TaskAL, QuantDay, Quant = 0.16)
      Q84Days <- lapply(TaskAL, QuantDay, Quant = 0.84)
      
      MQ16 <- round(mean(unlist(Q16Days)), digits = 0)
      MQ84 <- round(mean(unlist(Q84Days)), digits = 0)
      
      #Plot through the list
      #Make plot
      plot(x = 0:NROW(TaskAL[[1]]), 
           y = seq(from = 0, to = (max(unlist(TaskAL))*1.05), length.out =(NROW(TaskAL[[1]])+1) ),
           yaxt = "n", 
           type = "n", xlab = 'Days  Post-Assessment', 
           ylab = 'Cohort Contribution\n(Arbitrary Units)',
           main = paste0(CurTask,', Colony ',TCName,'\n Approx. Peak Day: ',PDM,'\n Approx. Q16 Day: ',MQ16, '\n Approx. Q84 Day: ',MQ84)
      )
      
      # Each replicate
      
      for(X in 1: length(TaskAL)){
        
        points(x = 1:NROW(TaskAL[[X]]), y = TaskAL[[X]], 
               type = 'l', lwd = 1.2, col = CC)
        
      }
      
      # Summative average
      
      points(x = 1:NROW(TaskAL[[X]]), y = AEmph(TaskAL), 
             type = 'l', lwd = 5, col = 'black')
      
      
      # Append emphasis list to storage list
      
      TaskEs[[N]] <- TaskAL
      
    }
    
    
    # For total number emerged
    
    AdultAL <- mapply(FUN = AdultAll, TrackA = TLCS, EmA = EmL, 
                      MoreArgs = list(MaxAge = RT),
                      SIMPLIFY = F)
    
    PeakAAL <- lapply(AdultAL, QuantDay, Quant = 0.5)
    
    PDMAAL <- round(mean(unlist(PeakAAL)), digits = 0)
    
    
    AALQ16Days <- lapply(AdultAL, QuantDay, Quant = 0.16)
    AALQ84Days <- lapply(AdultAL, QuantDay, Quant = 0.84)
    
    AALMQ16 <- round(mean(unlist(AALQ16Days)), digits = 0)
    AALMQ84 <- round(mean(unlist(AALQ84Days)), digits = 0)
    
    #Plot through the list
    
    #Make plot
    plot(x = 0:NROW(AdultAL[[1]]), 
         y = seq(from = 0, to = (max(unlist(AdultAL))*1.05), length.out =(NROW(AdultAL[[1]])+1) ),
         yaxt = "n", 
         type = "n", xlab = 'Days  Post-Assessment', 
         ylab = 'Number of Alive & Emerged Cohort\n(Whole Experiment)',
         main = paste0('Adults from Cohort, Colony ',TCName,'\n Approx. Peak Day: ',PDMAAL, '\n Approx. Q16 Day: ',AALMQ16, '\n Approx. Q84 Day: ',AALMQ84)
    )
    
    for(X in 1: length(AdultAL)){
      
      points(x = 1:NROW(AdultAL[[X]]), y = AdultAL[[X]], 
             type = 'l', lwd = 1.2, col = CC)
      
    }
    
    # Summative average
    points(x = 1:NROW(AdultAL[[X]]), y = AEmph(AdultAL), 
           type = 'l', lwd = 5, col = 'black')
    
    
    # Append to per-colony storage lists
    
    ColonyTB[[TC]] <- TaskEs
    ColonyAB[[TC]] <- AdultAL
    
  }
  
  
}

##### Next resolution #####

if (TRes == 'Experiment'  |  TRes == 'Both'){
  
  ## Create list to hold task emphases for cohort for each task
  
  TaskEs <- vector(mode = "list", length = NROW(TaskAges))
  
  for(N in 1:NROW(TaskAges)){
    
    CurTask <- TaskAges$Task[N]
    
    # Create list of relative contribution of target cohort to this task across experiment through time
    
    TaskAL <- mapply(FUN = TaskAll, TrackA = TrackL, EmA = EmL, 
                     MoreArgs = list( Task = CurTask, TaskList = TaskAges, MaxAge = RT),
                     SIMPLIFY = F)
    
    PeakDays <- lapply(TaskAL, QuantDay, Quant = 0.5)
    
    PDM <- round(mean(unlist(PeakDays)), digits = 0)
    
    Q16Days <- lapply(TaskAL, QuantDay, Quant = 0.16)
    Q84Days <- lapply(TaskAL, QuantDay, Quant = 0.84)
    
    MQ16 <- round(mean(unlist(Q16Days)), digits = 0)
    MQ84 <- round(mean(unlist(Q84Days)), digits = 0)
    
    #Plot through the list
    
    #Make plot
    plot(x = 0:NROW(TaskAL[[1]]), 
         y = seq(from = 0, to = (max(unlist(TaskAL))*1.05), length.out =(NROW(TaskAL[[1]])+1) ),
         yaxt = "n", 
         type = "n", xlab = 'Days  Post-Assessment', 
         ylab = 'Cohort Contribution\n(Arbitrary Units)',
         main = paste0(CurTask,', All Colonies\n','Approx. Peak Day: ',PDM,'\n Approx. Q16 Day: ',MQ16, '\n Approx. Q84 Day: ',MQ84)
    )
    
    # Each replicate
    
    for(X in 1: length(TaskAL)){
      
      points(x = 1:NROW(TaskAL[[X]]), y = TaskAL[[X]], 
             type = 'l', lwd = 1.2, col = CC)
      
    }
    
    # Summative average
    
    points(x = 1:NROW(TaskAL[[X]]), y = AEmph(TaskAL), 
           type = 'l', lwd = 5, col = 'black')
    
    
    # Append emphasis list to storage list
    
    TaskEs[[N]] <- TaskAL
    
  }
  
  
  # For total number emerged
  
  AdultAL <- mapply(FUN = AdultAll, TrackA = TrackL, EmA = EmL, 
                    MoreArgs = list(MaxAge = RT),
                    SIMPLIFY = F)
  
  PeakAAL <- lapply(AdultAL, QuantDay, Quant = 0.5)
  
  PDMAAL <- round(mean(unlist(PeakAAL)), digits = 0)
  
  AALQ16Days <- lapply(AdultAL, QuantDay, Quant = 0.16)
  AALQ84Days <- lapply(AdultAL, QuantDay, Quant = 0.84)
  
  AALMQ16 <- round(mean(unlist(AALQ16Days)), digits = 0)
  AALMQ84 <- round(mean(unlist(AALQ84Days)), digits = 0)
  
  #Plot through the list
  
  #Make plot
  plot(x = 0:NROW(AdultAL[[1]]), 
       y = seq(from = 0, to = (max(unlist(AdultAL))*1.05), length.out =(NROW(AdultAL[[1]])+1) ),
       yaxt = "n",
       type = "n", xlab = 'Days  Post-Assessment', 
       ylab = 'Number of Alive & Emerged Cohort\n(Whole Experiment)',
       main = paste0('Adults from Cohort, All Colonies \n', 'Approx. Peak Day: ',PDMAAL, '\n Approx. Q16 Day: ',AALMQ16, '\n Approx. Q84 Day: ',AALMQ84)
  )
  
  for(X in 1: length(AdultAL)){
    
    points(x = 1:NROW(AdultAL[[X]]), y = AdultAL[[X]], 
           type = 'l', lwd = 1.2, col = CC)
    
  }
  
  # Summative average
  points(x = 1:NROW(AdultAL[[X]]), y = AEmph(AdultAL), 
         type = 'l', lwd = 5, col = 'black')
  
  
}

#
### Storage lists can be used with easily modified code above to retrieve numbers embedded in
### graphics for further analysis / easier  exportation.

# Code by Dr. Lewis J. Bartlett. lewis.bartlett@uga.edu