# check and sort/analyze SST data
# This script find outliers and checks the data
# the output from this script should be passed to ss_mod_setup.R

rm(list=ls())

library(gamlss.dist)

# Set-up directories
# if(Sys.getenv('HOME')=='/home/scotti') {
#   ## means it is on a server
#   main_dir <- "/home/scotti/2022/behaveModelling" # LISA
#   data_dir <- file.path(main_dir, 'SST/rdatas')
#   res_dir <- file.path(main_dir, 'SST/samples')
#   pdf_dir <- file.path(main_dir, 'SST/pdfs')
#   cores <- 16
# } else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
#   main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/LISA/2022/behaveModelling'
#   data_dir <- file.path(main_dir, 'SST/rdatas')
#   res_dir <- file.path(main_dir, 'SST/samples')
#   pdf_dir <- file.path(main_dir, 'SST/pdfs')  
#   cores <- 4
# }

# setwd(main_dir)
main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/tux18/Public/trondheim/scripts/stopSignal/jointModelling/JointModelPMWG-master'


# setwd('/Users/scotti/surfdrive/Projects/7T_SST_MSIT_NTNU')

source("dmc/dmc.R")
load_model ("EXG-SS","exgSS.R")

load('data/dat_for_test_data.RData')

########## DATA STRUCTURE #########

# > head(dat)
### NEW
# R     S s  SSD        RT SS
# 11 RIGHT right 2 0.20 0.6176246 SS
# 21 RIGHT right 2  Inf 0.5684693 GO
# 30  LEFT  left 2  Inf 0.3685658 GO
# 41 RIGHT right 2  Inf 0.4989910 GO
# 50    NR right 2 0.15        NA SS
# 60  LEFT  left 2 0.20 0.5483389 SS

####### TEST ALL SUBJECTS FOR EXCLUSIONS AND OUTLIERS ######
dat$R <- factor(dat$R, levels = c("NR","LEFT","RIGHT")) # put factor levels in correct order
#sort_labs <- paste(sort(as.integer(levels(dat$s))))
#dat$s <- factor(dat$s, levels=sort_labs) # fix lexicographic sorting of numbers

remove_subs <- c() # create empty vector of subs to remove

str(dat)
## NEW
# 'data.frame':	6400 obs. of  6 variables:
#   $ R  : Factor w/ 3 levels "NR","LEFT","RIGHT": 3 3 2 3 1 2 2 1 1 3 ...
# $ S  : Factor w/ 2 levels "left","right": 2 2 1 2 2 1 1 2 1 2 ...
# $ s  : Factor w/ 33 levels "2","3","4","5",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ SSD: num  0.2 Inf Inf Inf 0.15 ...
# $ RT : num  0.618 0.568 0.369 0.499 NA ...
# $ SS : Factor w/ 2 levels "GO","SS": 2 1 1 1 2 2 1 2 2 1 ...

lapply(dat,levels)
# $R
# [1] "NR"    "LEFT"  "RIGHT"
# 
# $S
# [1] "left"  "right"
# 
# $s
# [1] "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18"
# [18] "19" "20" "21" "23" "24" "25" "26" "27" "29" "31" "32" "33" "34" "37" "38" "39"
# 
# $SSD
# NULL
# 
# $RT
# NULL
# 
# $SS
# [1] "GO" "SS"

# starting with .. subjects

dat <- dat[dat$RT >= 0.15 | is.na(dat$RT)==1,]   # Remove RTs < 200ms as these could not have be conscious decisions
dat <- dat[dat$RT <= 1.6 | is.na(dat$RT)==1,]
is.stop <- is.finite(dat$SSD) # find stop signal trials
gdat <- dat[!is.stop,]        # find go trial data
ssdat <- dat[is.stop,]        # find stop trial data

# Go omissions
is.nr <- gdat$R=="NR"         # which trials were go omissions on go trials 
# get rid of people with > 10 go omissions
sort(round(tapply(is.nr,gdat$s,mean)*100,2))
# 2     6    20    23    31    32    34     3     5    25    33    39    16    18    21 
# 0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.67  0.67  0.67  0.67  0.67  1.33  1.33  1.33 
# 24    37     7    12     9    11    15    17    26    29    14    38     4    10     8 
# 1.33  1.33  2.00  2.00  2.67  2.67  3.33  4.67  4.67  5.33  6.67  6.67  8.00  8.67 10.00 
# 19    13    27 
# 15.33 22.67 28.67 

remove_subs <- append(remove_subs, c(27, 19, 13)) # remove subs with > 8 go omissions

# Go accuracy
crct <- gdat$R == toupper(gdat$S)
round(sort(tapply(crct,gdat[,"s"],mean)),3)
# 27    13    19     8    10     4    14    38    29    17    26    11    15    18     9 
# 0.713 0.773 0.847 0.893 0.913 0.920 0.933 0.933 0.940 0.953 0.953 0.967 0.967 0.967 0.973 
# 7    12    16    23     5    21    24    33    37    39     3    25     2     6    20 
# 0.980 0.980 0.980 0.980 0.987 0.987 0.987 0.987 0.987 0.987 0.993 0.993 1.000 1.000 1.000 
# 31    32    34 
# 1.000 1.000 1.000 

remove_subs <- append(remove_subs, c(27,13,19)) # remove subs with < 85% accuracy

# check race model 
# remove subs if failed stops took longer than go trials
sort(round(tapply(gdat$RT,gdat$s,mean,na.rm=TRUE)-tapply(ssdat$RT,ssdat$s,mean,na.rm=TRUE),3))
# 13     27     23     34     16     10     19     32      8     39      6     20     33 
# -0.073 -0.024  0.011  0.019  0.028  0.029  0.032  0.047  0.052  0.054  0.059  0.071  0.075 
# 7     11      2     29     12     21      5     24     17     31      9     18      3 
# 0.076  0.079  0.082  0.083  0.084  0.086  0.092  0.095  0.097  0.103  0.104  0.104  0.119 
# 14     38     37     25     26     15      4 
# 0.120  0.121  0.125  0.137  0.139  0.140  0.176 

remove_subs <- append(remove_subs, c(13,27)) # remove subs with results < 0 (means that ss RTs were longer than go RTs)

# check stopping accuracy
resp <- ssdat$R =='NR'
round(sort(tapply(resp, ssdat[,"s"], mean)),2)
# 23   27   13   32   39   16   33   34    6    9   18    2    3   11   29   31    5   20 
# 0.42 0.44 0.48 0.48 0.48 0.50 0.50 0.50 0.52 0.52 0.52 0.54 0.54 0.54 0.54 0.54 0.56 0.56 
# 21   25   37   38   15    4   10   12   17   24    7    8   14   26   19 
# 0.56 0.56 0.56 0.56 0.58 0.60 0.60 0.60 0.60 0.60 0.62 0.62 0.62 0.62 0.66 

remove_subs <- append(remove_subs, c(19)) # remove subs either < 35% or > 65% accurate

###### EXTRA CHECKS #######

# plot all ihibition functions per participant
par(mfrow=c(2,3))
for (i in levels(dat$s)){
  tapply(!is.na(dat$RT[dat$s==i]),dat[dat$s==i,c("SS","SSD")],mean)
  plot_SS_if.dmc(dat[dat$s==i,], main = paste0("Inhibition function, sub - ", i))
}

# plot all staircases per participant
for (i in levels(dat$s)){
  plot(1:length(dat$SSD[dat$SS=="SS"&dat$s==i]), dat$SSD[dat$SS=="SS"&dat$s==i], type="l", main=paste0("sub - ", i, ", staircase over exp"))
}

# RTs over experiment per participant
for (i in levels(dat$s)){
  plot(1:length(dat$RT[dat$SS=="GO"&dat$s==i]), dat$RT[dat$SS=="GO"&dat$s==i], type="l", main=paste0("sub - ", i, ", RT over exp"))
}

# remove_subs <- append(remove_subs, c(40,43,52,56, 60, 13,22,38))

#####################################################################
######### REMOVE EXCLUSIONS FROM THE REST OF THE ANALYSIS ###########
dat <- dat[!dat$s %in% remove_subs,]
dat <- droplevels(dat)
#####################################################################
#####################################################################

par(mfrow=c(1,1))

# Go RT
head(sort(gdat$RT))
tail(sort(gdat$RT))
hist(gdat$RT,breaks="fd")

# Stop RT
head(sort(ssdat$RT))
tail(sort(ssdat$RT))
hist(ssdat$RT,breaks="fd")

# Show the different SSDs:
sort(tapply(as.character(dat$SSD),dat[,c("SS")],unique)$SS)
# [1] "0.05" "0.1"  "0.15" "0.2"  "0.25" "0.3"  "0.35" "0.4"  "0.45" "0.5"  "0.55" "0.6"  "0.65" "0.7"  "0.75" "0.8" 

# Show the number of trials for each SSD:
Ns = tapply(dat$RT,dat$SSD,length)
Ns
# 0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8  Inf 
# 19   36   76  110   93   85   88  100  106   95   82   63   51   44   24    3 3224 

# Show response rate:
tapply(!is.na(dat$RT),dat[,c("SS")],mean)
# GO        SS 
# 0.9776675 0.4455814 

# Response rate broken down by SSD & corresponding inhibition function:
tapply(!is.na(dat$RT),dat[,c("SS","SSD")],mean)
plot_SS_if.dmc(dat)  #P(Respond) increases as a function of SSD, as it should

# Plot median signal-respond RT per SSD:
tapply(dat$RT,dat$SSD,median,na.rm=TRUE) 
plot_SS_srrt.dmc(dat) # Median SRRT increases as a function of SSD, as it should

# Show number of signal-respond RTs per SSD:
Nr = tapply(!is.na(dat$RT),dat[,c("SS","SSD")],sum)[2,]
Nr
# 0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8  Inf 
# 6   13   24   53   38   37   33   41   49   47   40   34   21   22   18    3   NA 

# Signal-respond RTs (stop dist, red line) should be faster than go RTs:
hist(dat$RT[dat$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8)) # Go dist
lines(density(dat$RT[dat$SS=="SS"],na.rm=T),col="red",lwd=2) # Stop dist

# test the racce model on one subject
test_race <- function(s,ssdat,gdat) {
  grt <- gdat[gdat$s==s,"RT"]
  srt <- ssdat[ssdat$s==s,"RT"]
  print(t.test(grt[!is.na(grt)],srt[!is.na(srt)]))
}

test_race("2",ssdat,gdat) # test subject 2

# set it up for pmwg
dat$S <- as.character(dat$S)
dat$S[dat$S=='right'] <- 2
dat$S[dat$S=='left'] <- 1
dat$R <- as.character(dat$R)
dat$R[dat$R=='RIGHT'] <- 2
dat$R[dat$R=='LEFT'] <- 1
dat$R[dat$R=='NR'] <- 3

dat <- dat %>% rename(direction = S, SS = SS, subject = s, RT = RT, R = R, SSD = SSD)

save(dat, file='data/prepped_data_sst.RData')

# subjects 27 19 13 removed
#   dat$correct <- as.numeric(dat$S)==as.numeric(dat$R)-1 # make correct oclumn for ease

R direction subject  SSD        RT SS
11 2         2       2 0.20 0.6176246 SS
21 2         2       2  Inf 0.5684693 GO
30 1         1       2  Inf 0.3685658 GO
41 2         2       2  Inf 0.4989910 GO
50 3         2       2 0.15        NA SS
60 1         1       2 0.20 0.5483389 SS

'data.frame':	4299 obs. of  6 variables:
  $ R        : chr  "2" "2" "1" "2" ...
$ direction: chr  "2" "2" "1" "2" ...
$ subject  : Factor w/ 22 levels "2","3","4","5",..: 1 1 1 1 1 1 1 1 1 1 ...
$ SSD      : num  0.2 Inf Inf Inf 0.15 ...
$ RT       : num  0.618 0.568 0.369 0.499 NA ...
$ SS       : Factor w/ 2 levels "GO","SS": 2 1 1 1 2 2 1 2 2 1 ...




