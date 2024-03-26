rm(list=ls())


# remove 23, 25, 29

library(dplyr)
library(ggplot2)
library(stringr)
library(gridGraphics)
library(grid)
library(gridExtra)
library(plotrix)
library(IDPmisc)

options(max.print=1000000)
options(warn=-1)


# define function to calculate standard error
SE <- function(x) {
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

# round dataframe columns function
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

#### MSIT
# General info: 200 trials (22 null trials)
# 7s trials

main_dir <- getwd() # where is the data?
res_dir <- file.path(getwd(), 'data') # where do you want to save things?

# msit_paths <- Sys.glob(file.path(main_dir,  'zipdata', 'sub-*','ses-sstmsit','sub-*_ses-sstmsit_data_tmpunzip',
#                                  'sub-*_ses-sstmsit_data','RESOURCES','task_data','*_task-MSIT_*_block-*_events.tsv'))
msit_paths <- Sys.glob(file.path(main_dir, 'all_behavioural_data_NTNU', '*_task-MSIT_*_block-*_events.tsv'))

msit_dat <- c()
for (i in 1:length(msit_paths)){
  
  dat_tsv <- read.table(msit_paths[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
  if ('rt' %in% colnames(dat_tsv)){
    msit_dat <- rbind(msit_dat, dat_tsv)}
  
}

msit_dat <- msit_dat %>% dplyr::select(trial_nr, block_nr, subject, null_trial, event_type, stimuli, response,
                                condition, choice_key, rt, correct_response, block_nr)


msit_dat$subject <- as.factor(msit_dat$subject)                               # make each subjectNumber a factor
msit_dat$trial_nr <- msit_dat$trial_nr+1
msit_dat <- msit_dat[msit_dat$null_trial==0,]

# extract one line per trial (this is probably overcomplicated... but it works.. )
msit_dat <- msit_dat[msit_dat$event_type == 'stimulus' | (msit_dat$event_type=='response' & !is.na(msit_dat$rt)),] # take only stimlus and reponse even types
ind=grep("stimulus",msit_dat$event_type) # index rows that are event type stimulus
msit_dat <- msit_dat[sort(c(ind[c(0,abs(diff(ind))==1)==1]-1,grep("response",msit_dat$event_type))),] # only keep stimulus rows when no response is given
msit_dat$response <- as.character(msit_dat$response)
msit_dat$response[msit_dat$event_type=='stimulus'] <- 'None'

# change some other stuff
msit_dat$correct_response <- as.character(msit_dat$correct_response)
msit_dat$condition <- as.factor(msit_dat$condition)                           # make condit a factor 
msit_dat <- msit_dat[msit_dat$choice_key!='NULL',]           # remove trials where the subjectNumber did not respond

msit_dat$correct <- 0                                                         # create accuracy colmun
msit_dat$correct[(msit_dat$choice_key=='r' & msit_dat$correct_response=='1') | (msit_dat$choice_key=='g' & msit_dat$correct_response=='2') | 
                   (msit_dat$choice_key=='y' & msit_dat$correct_response=='3')] <- 1 # fill accuracy column

# 
msit_dat <- msit_dat %>% dplyr::select(trial_nr, subject, stimuli, condition, correct_response, rt, choice_key, correct, block_nr, response)
msit_dat$choiceTime <- choiceTime <- 0.6 # soft deadline

str(msit_dat)
# 'data.frame':	5008 obs. of  10 variables:
#   $ trial_nr     : num  1 2 3 4 5 6 7 8 9 10 ...
# $ subjectNumber: Factor w/ 21 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ stimuli      : int  30 223 212 30 322 223 211 322 211 223 ...
# $ condit       : Factor w/ 4 levels "1","2","3","4": 2 3 4 2 4 3 4 4 4 3 ...
# $ uniq         : chr  "3" "3" "1" "3" ...
# $ rt           : num  0.673 0.819 0.73 0.638 0.641 0.475 0.6 0.591 0.964 0.603 ...
# $ key_press    : chr  "2" "3" "1" "3" ...
# $ correct      : num  0 1 1 1 1 1 0 1 1 1 ...
# $ block        : int  1 1 1 1 1 1 1 1 1 1 ...
# $ session      : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...

msit_dat <- unique(msit_dat) # remove duplcaited rows caused by saving of behavioural file twice
msit_dat <- droplevels(msit_dat)

# subject R    RT pos resp condit uniq flank nstims stimuli
# 291       2 2 0.666   2    2      1    2     2      1     020
# 292       2 1 0.866   1    1      1    1     1      1     100
# 293       2 2 0.698   1    2      2    2     2      1     200
# 294       2 2 0.708   1    2      4    2     1      3     211
# 295       2 2 0.585   2    2      4    1     2      3     212
# 296       2 1 0.602   1    1      1    1     1      1     100

# sort data for nieks code

msit_dat$uniq <- msit_dat$correct_response
msit_dat$response[msit_dat$response=='y'] <- '3'
msit_dat$response[msit_dat$response=='r'] <- '1'
msit_dat$response[msit_dat$response=='g'] <- '2'
msit_dat <- msit_dat[msit_dat$response!="None",]
msit_dat$R <- as.factor(msit_dat$response)
msit_dat$resp <- as.numeric(msit_dat$response)
msit_dat$RT <- round(msit_dat$rt,3)
msit_dat$condit <- msit_dat$condition
msit_dat$stimuli <- as.character(formatC(msit_dat$stimuli, width=3,format='d',flag='0'))

msit_dat$nstims <- nchar(sub("0", "", msit_dat$stimuli))
msit_dat$nstims[msit_dat$nstims=='2'] <- '1'

msit_dat$pos <- 0
for (i in 1:length(msit_dat$uniq)){
  # for (i in 1:5){
    # print(msit_dat$uniq[i])
    # print(msit_dat$stimuli[i])
msit_dat$pos[i] <- gregexpr(pattern = msit_dat$uniq[i], msit_dat$stimuli[i])[[1]][1]
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


msit_dat$flank <- 0
for (i in 1:length(msit_dat$uniq)){
  if (msit_dat$nstims[i] == 1){
    msit_dat$flank[i] <- msit_dat$uniq[i]
  } 
  else {
    msit_dat$flank[i] <- Mode(str_split(msit_dat$stimuli[i], pattern="", n = Inf, simplify = FALSE)[[1]])
  }
}

msit_dat <- msit_dat %>% dplyr::select(subject, R, RT, pos, resp, condit, uniq, flank, nstims, stimuli)

msit_dat <- msit_dat[!msit_dat$subject %in% c(23,25,29),]

msit_dat$pos <- as.integer(msit_dat$pos)
msit_dat$condit <- as.integer(msit_dat$condit)
msit_dat$uniq <- as.integer(msit_dat$uniq)
msit_dat$flank <- as.numeric(msit_dat$flank)
msit_dat$nstims <- as.numeric(msit_dat$nstims)

msit_dat <- droplevels(msit_dat)

save(msit_dat, file='./data/data_MSIT_NTNU_prepped_final.RData')
