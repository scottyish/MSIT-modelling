# Setup -------------------------------------------------------------------
rm(list=ls())

# Set-up directories
if(Sys.getenv('HOME')=='/home/scotti') {
  ## means it is on a server
  main_dir <- "/home/scotti/2022/behaveModelling/pmwg" # LISA
  cores <- 16
} else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
  main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/LISA/2022/behaveModelling/pmwg'
  cores <- 4
}

print(main_dir)
setwd(main_dir)

library(pmwg)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

source("utils/dataPrep.R")
source("utils/diagnostics.R")
source("utils/postPredict.R")
source("utils/run.R")

source("models/RW/RD.R")
source("models/RW/dists_ARD.R")
source("models/RW/dists.R")

source("models/SST/EXG.R")
source("models/SST/dists.R")

# source("pmwg/sampling_factors.R")
source("pmwg/sampling.R")
source("jointUtils/jointRun.R")


# library(dmcAdapt)


#Name: data name to load
#SubNCol: Where we'll look for the subject identifier
#accCol: Response correctness column in data, used to filter data
#accCriterion: Filter data based on accuracy per subject
#RL: Whether to employ RL learning rules in drift rate calculations (for now assumes an advantage framework)
#RTLimits: RTs above and below these will be filtered out
#Factors: Columns of the original data that should be included in the sampler, should also include the condition on which response correctness is determined if response-coded
#Match: a list with as name the parameter that should have matched entries in the sampler. If response coded first vector entry is the response, second vector entry is whether that response is correct.
#Match: if accuracy coded, only provide accuracy column. 
#Constants: Constants in the parameters that you want to add, e.g. to satisfy scaling constraints
#ParNames: The names of the parameters you want to include in the sampler. This has quite a specific format:
#Format explained: First parameter is the 'base' parameter, then if included in match _1 means incorrect response, _2 correct (DDM convention).
#continued, For all parameters add underscore for potential other factors and to which level of this factor it belongs to. Add additional factors with additional underscore. Separate factor and level with .
#If factor is the first entry of the match vector, it will assume that there should be separate values for the accumulators. e.g. b_response.1 means separate boundaries for the accumulators, this one for accumulator 1.
#llFunc, the likelihood function to use, for now only racing diffusion (Wald) included.
#priorMean: vector of lenght parnames with priors in order. Assumes normal priors and pars so transform in LL if necessary.
#startPoints: Sensible startpoints to start sampling. (I was lazy)

modFunc <- function(a, b){(1 + a)*b}
plusFunc <- function(a,b){a+b}
simplexFunc <- function(b){1 -b}

experiments <- list(
  "MSIT" = list(name = "MSIT",
                subNCol = "subject",
                accCol = "correct",
                rtCol = 'rt',
                accCriterion = 0.5,
                "RL" = F,
                RTLimits = c(0.15, 2),
                factors = c("pos", "resp", "condit", "uniq", "flank", "stimuli"),
                transFunc = list(func = driftsMSIT, earlyTransform = F),
                constants = list(s = 1, A = 0, vPos.3 = 0),
                match = list(vPos = c("resp", "uniq")),
                parNames = c("t0","vFlank", "vSimon", "vPos.1", "vPos.2", "vMatch", "v", "B"),
                llFunc = likelihood.RD,
                priorMean =  c(0.15, rep(0.1, 7)),
                startpoints = rep(2,8), # used to be NULL but wasnt running
                modelName = "MSIT_s_A_vPos3"),
  "SST" = list(name = "SST", 
               rtCol = 'RT', 
               subNCol = "s", 
               accCol = "correct", 
               accCriterion = 0.5,
               "RL" = F, 
               RTLimits = c(0.15, 2), 
               factors = c("SSD", "R", 'SS', "direction"),
               transFunc = list(func = parsSST, earlyTransform = F),
               constants = list(mu_1=log(1e6),sigma_1=log(0.001),tau_1=log(0.001)), #gf = qnorm(.000001)), 
               match = list(mu = c("R", "direction"), sigma = c('R', 'direction'), tau = c('R', 'direction')),
               parNames = c("mu_2", "sigma_2", "tau_2", "muS", "sigmaS", "tauS", "tf", "gf"),
               llFunc = likelihood.EXG, 
               priorMean = rep(0, 8), 
               startpoints = rep(2, 8),
               modelName = "SST_gftf_normstartpoints")
)
experiments <- prepDataJoint_SSTMSIT(experiments, splitSess = F, hardCutAcc = T, hardCutN = T, n_sessions = 2) #If splitsess, we can run a joint model comparing session one with session two for the current task.

# sharedPars <- list(t0 = list(priors = -1, startpoints = NULL))
set.seed(1234)
pmwgJointRun(experiments, pstar = .7, n_factors = 5, n_cores = cores) #, sharedPars = sharedPars)

