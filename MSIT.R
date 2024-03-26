# Setup -------------------------------------------------------------------
rm(list=ls())

# setwd('/home/scotti/Downloads/behaveModelling')
# res_dir <- '/home/scotti/Downloads/behaveModelling/SST2/rdatas'
# cores = 8

# Set-up directories
if(Sys.getenv('HOME')=='/home/scotti') {
  ## means it is on a server
  main_dir <- "/home/scotti/2022/behaveModelling/pmwg" # LISA

  cores <- 16
} else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
  main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/LISA/2022/behaveModelling/pmwg'

  cores <- 4
}

setwd(main_dir)

library(dplyr, warn.conflicts = FALSE)
library(bayesplot)
options(dplyr.summarise.inform = FALSE)
# setwd("~/RL-PMWG/")

source("utils/dataPrep.R")
source("utils/diagnostics.R")
source("utils/postPredict.R")
source("utils/run.R")
source("models/RW/RD.R")
source("models/RW/dists_ARD.R")
source("models/RW/dists.R")
source("pmwg/sampling.R")

### PRIORS FROM ONLINE EXP ####
# t0          vFlank     vSimon     vPos.1     vPos.2     vMatch          v          B 
# 0.21365281 1.11009427 1.11507915 0.23113828 0.08102416 2.46968668 0.65054216 1.72542060 

# nieks data
# exp1 <- list("MSIT" = list(name = "MSIT", rtCol = 'rt', subNCol = "subjectNumber", accCol = "correct", accCriterion = 0.5,
#                                    "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq", "flank", "nstims", "stimuli"),
#                                    transFunc = list(func = driftsMSIT, earlyTransform = F),
#                                    constants = list(s = 1, A = 0, vPos.3 = 0), match = list(v = c("resp", "uniq")),
#                                    parNames = c("t0", "vFlank", "vSimon", "vPos.1", "vPos.2", "vMatch", "v", "B"),
#                                    llFunc = likelihood.RD, priorMean = c(-0.5, rep(0.1, 7)), startpoints = NULL,
#                                    modelName = "MSIT_procs"))
# 
# exp1 <- prepDataJoint(exp1, splitSess = F, hardCutAcc = T, hardCutN = T, n_sessions = 1) #If splitsess, we can run a joint model comparing session one with session two for the current task.

# debug(runSampler)
# pmwgRun(exp1, n_cores = 4)


# # ntnu data, original model            
# exp <- list("MSIT" = list(name = "MSIT", rtCol = 'rt', subNCol = "subject", accCol = "correct", accCriterion = 0.5,
#                           "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq", "flank", "nstims", "stimuli"),
#                           transFunc = list(func = driftsMSIT, earlyTransform = F),
#                           constants = list(s = 1, A = 0, vPos.3 = 0), match = list(v = c("resp", "uniq")),
#                           parNames = c("t0", "vFlank", "vSimon", "vPos.1", "vPos.2", "vMatch", "v", "B"),
#                           llFunc = likelihood.RD, priorMean = c(-0.5, rep(0.1, 7)), startpoints = NULL,
#                           modelName = "MSIT_procs")
# )

# # ntnu data, fix v baseline to 1 and remove A=0 as and v.pos3 as constants
# exp <- list("MSIT" = list(name = "MSIT", rtCol = 'rt', subNCol = "subject", accCol = "correct", accCriterion = 0.5,
#                           "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq", "flank", "nstims", "stimuli"),
#                           transFunc = list(func = driftsMSIT, earlyTransform = F),
#                           constants = list(s = 1, v = 1), match = list(v = c("resp", "uniq")),
#                           parNames = c("t0", "vFlank", "vSimon", "vPos.1", "vPos.2", "vPos.3", "vMatch", "B","A"),
#                           llFunc = likelihood.RD, priorMean = c(-0.5, rep(0.1, 8)), startpoints = NULL,
#                           modelName = "MSIT_procs")
# )

# ntnu data, fix v baseline to 1 and remove v.pos3 as constant
exp <- list("MSIT" = list(name = "MSIT", rtCol = 'rt', subNCol = "subject", accCol = "correct", accCriterion = 0.5,
                          "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq", "flank", "nstims", "stimuli"),
                          transFunc = list(func = driftsMSIT, earlyTransform = F),
                          constants = list(s = 1, v = 0, vPos.3 = 0, A = 0), match = list(v = c("resp", "uniq")),
                          parNames = c("t0", "vFlank", "vSimon", "vPos.1", "vPos.2", "vMatch", "B"),
                          llFunc = likelihood.RD, priorMean = c(0.15, rep(0.1, 6)), startpoints = NULL,
                          modelName = "MSIT_procs")
)

exp <- prepDataJointNTNU(exp, splitSess = F, hardCutAcc = T, hardCutN = T, n_sessions = 1) #If splitsess, we can run a joint model comparing session one with session two for the current task.

# debug(runSampler)
pmwgRun(exp, n_cores = cores)


# to make plots after
# postCheck(sampled$experiments, PP = T, n = 50)
