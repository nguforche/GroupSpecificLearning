#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
print("run code as:  /infodev1/infotools/R/3.4.3/bin/Rscript DataSetup.R 3")
stop("Please supply input number of months interval")
}

### ./configure --prefix=$che_dir/R/R-3.5.1
##  make && make install 
## install.packages(c("Matrix", "arulesSequences", "plyr", "lubridate", "data.table"))

#install.packages(c("arules","arulesSequences","plyr", "ggplot2", "caret", "gbm", "randomForest", "e1071", "devtools", "glmet", "missForest", "PresenceAbsence"))
## start R  /infodev1/infotools/R/bin/R
## \\hsrntfs\projects\hcpr\activity\s205373.group.specific.learning\che_cluster_analysis\excel 

## mnode of a vector 
getMode = function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
library(plyr)
library(lubridate)
suppressMessages(require(data.table))

suppressMessages(require(parallel))
suppressMessages(require(doParallel))
cores <-floor(0.80*detectCores())
registerDoParallel(cores=cores)


##load(file = "/infodev1/non-phi-projects/che/Research/GroupSpecificLearning/data/AnalyticData.RData")
### three months data k = 3 

k = as.numeric(args[1])

dir <- "/infodev1/non-phi-data/che/GroupSpecificLearning/data/newdata/"
base <- paste0(k, "month_20180226.csv")

dia.file <- "diabetes_"
htn.file <- "hypertension_"
hcl.file <- "hypercholesterolemia_"
ob.file <- "obesity_"
vt.file <- "vitals_"
lab.file <- "labs_"



orig <- ymd("2004-01-01")
mid <- ymd("2010-12-31")
mid2 <- ymd("2011-01-01")
end <- ymd("2015-12-31")

### diabetes 
db <- fread(file = paste0(c(dir, dia.file, base), collapse = ""), na.strings = "")
dia.vars <-  c("DM1", "DM2_1",   "DM2_2",     "DM3_1",       "DM3_2",       "DM3_3",  "DM3_4", 
                "DM4_1",    "DM4_2",              "DM4_3",          "DM4_4")
                
db$date <- ymd(paste0(paste0(paste0(db$startyear, "-"), db$startmon), "-1"))
db <- db[(date >= orig) & (date <= mid)]
db$time <-  as.numeric(difftime(db$date, orig, units = "days"))               

### hypertension 
htn <- fread(file = paste0(c(dir, htn.file, base), collapse = ""), na.strings = "")
htn.vars <- c("HTN1", "HTN2_1",   "HTN2_2",     "HTN3_1",       "HTN3_2",       "HTN3_3",  "HTN3_4", 
               "HTN4_1",    "HTN4_2",    "HTN4_3",      "HTN4_4")


htn$date <- ymd(paste0(paste0(paste0(htn$startyear, "-"), htn$startmon), "-1"))
htn <- htn[(date >= orig) & (date <= mid)]
htn$time <-  as.numeric(difftime(htn$date, orig, units = "days"))               

### Hypercholesterolemia
hcl <- fread(file = paste0(c(dir, hcl.file, base), collapse = ""), na.strings = "")
hcl.vars <- c("L1", "L2_1",   "L2_2",     "L3_1",       "L3_2",       "L3_3",  "L3_4",  "L4_1",   
             "L4_2",              "L4_3",          "L4_4")
hcl$date <- ymd(paste0(paste0(paste0(hcl$startyear, "-"), hcl$startmon), "-1"))
hcl <- hcl[(date >= orig) & (date <= mid)]  
hcl$time <-  as.numeric(difftime(hcl$date, orig, units = "days"))               

#### Obesity 
ob <- fread(file = paste0(c(dir, ob.file, base), collapse = ""), na.strings = "")
ob.vars <- c("OB1", "OB2_1",   "OB2_2",     "OB3_1",       "OB3_2",       "OB3_3",  "OB3_4",  "OB3_5")
ob$date <- ymd(paste0(paste0(paste0(ob$startyear, "-"), ob$startmon), "-1"))
ob <- ob[(date >= orig) & (date <= mid)]  
ob$time <-  as.numeric(difftime(ob$date, orig, units = "days"))               

### labs and vitals 
#### read in vitals and continue 
uits <- paste0(k, " months")
vt <- fread(file = paste0(c(dir, vt.file, base), collapse = ""), na.strings = "")
lb <- fread(file = paste0(c(dir, lab.file, base), collapse = ""), na.strings = "")
vt$date <- ymd(paste0(paste0(paste0(vt$startyear, "-"), vt$startmon), "-1"))
lb$date <- ymd(paste0(paste0(paste0(lb$startyear, "-"), lb$startmon), "-1"))
### labs and vitals in 2004-2010 
vt <- vt[(date >= orig) & (date <= mid)] 
lb <- lb[(date >= orig) & (date <= mid)]
vt$time <-  as.numeric(difftime(vt$date, orig, units = "days"))               
lb$time <-  as.numeric(difftime(lb$date, orig, units = "days")) 

### join data by pid and date 

setkey(db, pid, time)
setkey(htn, pid, time)
setkey(hcl, pid, time)
setkey(ob, pid, time)
setkey(vt, pid, time)
setkey(lb, pid, time)


dd <- db[htn, nomatch=0]
dd <- dd[hcl, nomatch=0]
dd <- dd[ob, nomatch=0] 
dd <- dd[vt, nomatch=0] 
dd <- dd[lb, nomatch=0] 

#######################################################################################################
#######################################################################################################
################## Project Goal: Compare predictive abaility of Sequential Patterns and Baseline  Cormobidities
#######################################################################################################
#######################################################################################################
#######################################################################################################
######## 
#### build sequential patterns in 2004/01/01 - 2010/12/31 to predict survival in 2011/01/01 - 2015/12/31
#### patient is alive in 2004/01/01 - 2010/12/31  but dead or alive in 2011/01/01 - 2015/12/31   
#### 
#### Build baseline comorbidities in 2004/01/01 - 2010/12/31
#### Three types: (1) any status (2) last status, and (3) most frequent status in 2004-2010 

my.max <- function(x) ifelse( length(x) > 0, max(x, na.rm = TRUE), 0)  
risk <- c(1, 1, 2, 1, 2, 2, 4, 2, 4, 4, 8)
ob.risk <- c(1, 1, 2, 1, 2, 2, 3, 4) 

Assign.risk  <- function(xx){
n <- NROW(xx)
nd <- do.call(rbind, lapply(seq_len(n), function(ii) {
DM = xx[ii, "DM1"]; HTN = xx[ii, "HTN1"]; HCL = xx[ii, "L1"]; OB = xx[ii, "OB1"]
tm = xx[ii, "time"]
dm <- my.max(risk[as.logical(xx[ii, dia.vars])]) 
ht <- my.max(risk[as.logical(xx[ii, htn.vars])])
hc <- my.max(risk[as.logical(xx[ii, hcl.vars])])
oc <- my.max(ob.risk[as.logical(xx[ii, ob.vars])]) 
return(cbind(time = tm , DM.risk = dm, HTN.risk= ht, HCL.risk = hc, 
         OB.risk = oc, DM = DM , HTN = HTN, HCL = HCL, OB = OB))
}) )
tb <- xx[order(xx$time,  decreasing = TRUE), ] 
## DM 
dm <- tb$DM1; DM.last = rep(dm[1],n); DM.freq = rep(getMode(dm),n); DM.any <- rep(max(dm),n)
## HTN 
ht <- tb$HTN1; HTN.last = rep(ht[1],n); HTN.freq = rep(getMode(ht),n); HTN.any <- rep(max(ht),n)
## HCL
hc <- tb$L1; HCL.last = rep(hc[1],n); HCL.freq <- rep(getMode(hc),n); HCL.any <- rep(max(hc),n)
## Obesity 
oc <- tb$OB1; Obese.last = rep(oc[1],n); Obese.freq = rep(getMode(oc),n); Obese.any <- rep(max(oc),n)
dn <- cbind(DM.last, DM.freq, DM.any, HTN.last, HTN.freq, HTN.any, HCL.last, HCL.freq, 
           HCL.any, Obese.last, Obese.freq, Obese.any)
cn <- xx[ cont.vars] 
cbind(nd, dn, cn)  
}

### Nodes where thre risk score can be either 1 or 0 
## the risk score for these nodes is equivalent to thre risk 
## score of the root node 
dt <- as.data.frame(dd)
dt[, dia.vars]["DM2_1"] <- dt[, dia.vars]["DM1"]
dt[, dia.vars]["DM3_1"] <- dt[, dia.vars]["DM1"]
dt[, htn.vars]["HTN2_1"] <- dt[, htn.vars]["HTN1"]
dt[, htn.vars]["HTN3_1"] <- dt[, htn.vars]["HTN1"]
dt[, hcl.vars]["L2_1"] <- dt[, hcl.vars]["L1"]
dt[, hcl.vars]["L3_1"] <- dt[, hcl.vars]["L1"]
dt[, ob.vars]["OB2_1"] <- dt[, ob.vars]["OB1"]
dt[, ob.vars]["OB3_1"] <- dt[, ob.vars]["OB1"]

cont.vars <- c(setdiff(names(lb), c("pid","startyear","endyear","startmon", "endmon", "date", "time")), 
              setdiff(names(vt), c("pid","startyear","endyear","startmon", "endmon", "date", "time")))
vars <- c("pid", "time", dia.vars, htn.vars, hcl.vars, ob.vars, cont.vars)

dd1 <- ddply(dt[, vars], .variables = "pid", .fun = Assign.risk, .parallel = TRUE)


### merge with demographics 
dat.demo <- read.csv(file = "/infodev1/non-phi-data/che/GroupSpecificLearning/data/demographics_20180604.csv", na.strings = "")
######## Data structure: sequences in 2004-2010 outcomes in 2011-2015
### take all patients born before 2004 and alive after 2004 
#### patient is alive in 2004/01/01 - 2010/12/31  but dead or alive in 2011/01/01 - 2015/12/31

ix.alive <- is.na(dat.demo$death_dt) 
dead <- dat.demo[!ix.alive, ]
alive <- dat.demo[ix.alive, ]
demo <- rbind(subset(dead, death_dt >= mid2), alive)

demo$birth_dt <- as.Date(strptime(demo$birth_dt,format="%Y-%m-%d"))
demo <- subset(demo, birth_dt < orig)

demo$dead <- 1*(!is.na(demo$dead))
demo$mortality <- 1*(!is.na(demo$death_dt)) 
demo$death_dt[demo$mortality == 0] <- "2015-12-31"  ## censor alive patients 

demo$death_dt <- as.Date(strptime(demo$death_dt,format="%Y-%m-%d"))


## compute survival follow up time 
## date last seen is date of death or end point: this can be improved using date of last vital measurement 
demo$death.time <- as.numeric(difftime(demo$death_dt, orig, units = "days")) 
## age 
demo$age <-  as.numeric(0.00273973*as.numeric(difftime(orig, demo$birth_dt,units = "days")))               
demo$age.Jan2011 <-  as.numeric(0.00273973*as.numeric(difftime(mid2, demo$birth_dt,units = "days")))               

### major cardiovascular events: MCV  
### MCVD models: for patients with any cdv event, date last seen = date of cvd. 
### for patients with no cvd event and still alive the date last seen = death of cvd event end point. 
### for patients with no cvd event and not alive, the date last seen is the date of death  

dx.cvd <- names(demo)[grep("DX2_all_", names(demo))]
dx.cvd.dt <- dx.cvd[grep("*_dt", dx.cvd)]
dx.cvd.code <- dx.cvd[grep("*_code", dx.cvd)]
cvd.vars <- setdiff(dx.cvd, c(dx.cvd.dt, dx.cvd.code))
demo$CVD <- apply(demo[, cvd.vars], 1, max, na.rm = TRUE)

## if a dx.cvd.dt is missing, patient does not have condition and replace with death_dt 
demo[, c(dx.cvd.dt)] <- lapply(demo[, c(dx.cvd.dt)], function(xx){ 
xx <- as.Date(strptime(xx,format="%m/%d/%Y"))
ix <- is.na(xx)
xx[ix] <- demo[ix, "death_dt"] 
xx
})

#cvd.time <- paste0(cvd.vars, ".time")
#demo[, cvd.time] <- lapply(demo[, dx.cvd.dt[-1]], function(xx){ 
#as.numeric(difftime(xx, orig, units = "days"))
#})

demo$cvd.all.dt <- apply(demo[, dx.cvd.dt[-1]], 1, min, na.rm = TRUE)
demo$cvd.all.time <- as.numeric(difftime(demo$cvd.all.dt, orig, units = "days"))

### take patients who are in demo and also in dat 
dat <- join(demo, dd1, by = "pid", type = "inner")


# ensure that once someone dies or have a cvd event, he or she stays that way 
## variables to work with: mortality, CVD, death.time, and cvd.all.time 

Mortality <- function(xx){
## mortality 
n <- NROW(xx)
#ix <- which(xx$death.time[1] >= xx$time) 
xx$death.time1 = xx$time
xx$death.time2 <- c(xx$time[-1], xx$death.time[1])
if(unique(xx$mortality) == 1){
xx$mortality[1:(n-1)] = 0 
xx$mortality[n] <- 1
yy <- xx[n, ]
yy$time = xx$death.time[1]
yy$death.time1 <- yy$time
yy$death.time2 <- yy$time
xx <- rbind(xx, yy)
}
xx
}

### first ensure all cvd.all.time > 0 before running this code 
CVD <- function(xx){
## mortality 
n <- NROW(xx)
#ix <- which(xx$death.time[1] >= xx$time) 
xx$cvd.all.time1 = xx$time
xx$cvd.all.time2 <- c(xx$time[-1], xx$cvd.all.time[1])
if(unique(xx$CVD) == 1){
xx$CVD[1:(n-1)] = 0 
xx$CVD[n] <- 1
yy <- xx[n, ]
yy$time = xx$cvd.all.time[1]
yy$cvd.all.time1 <- yy$time
yy$cvd.all.time2 <- yy$time
xx <- rbind(xx, yy)
}
xx
}

### mortality data 
dat.mort <- ddply(dat, .variables = "pid", .fun = Mortality, .parallel = TRUE)

### cvd data 
dat.cvd <- subset(dat, cvd.all.time >=0)
dat.cvd <- ddply(dat.cvd, .variables = "pid", .fun = CVD, .parallel = TRUE)


######### Aggregated data set
########## Agregate labs to combine with demographics 
Aggregated.Data  <- function(xx){
n <- NROW(xx)
com <- cbind(DM.last=xx$DM.last[1], DM.freq=xx$DM.freq[1], DM.any=xx$DM.any[1],  
         HTN.last=xx$HTN.last[1], HTN.freq = xx$HTN.freq[1], HTN.any = xx$HTN.any[1], HCL.last=xx$HCL.last[1], 
         HCL.freq = xx$HCL.freq[1], HCL.any = xx$HCL.any[1], Obese.last = xx$Obese.last[1], 
         Obese.freq = xx$Obese.freq[1], Obese.any = xx$Obese.any[1])

tb <- xx[order(xx$time,  decreasing = TRUE), ] 
dm <- tb$DM.risk; DM.risk.last = dm[1]; DM.risk.sum = sum(dm); DM.risk.ave = mean(dm);
ht <- tb$HTN.risk; HTN.risk.last = ht[1]; HTN.risk.sum = sum(ht); HTN.risk.ave = mean(ht)
ht <- tb$HCL.risk; HCL.risk.last = ht[1]; HCL.risk.sum = sum(ht); HCL.risk.ave = mean(ht)
ht <- tb$OB.risk; OB.risk.last = ht[1]; OB.risk.sum = sum(ht); OB.risk.ave = mean(ht)
risk <- cbind(DM.risk.last, DM.risk.sum, DM.risk.ave, HTN.risk.last, HTN.risk.sum, HTN.risk.ave, 
             HCL.risk.last, HCL.risk.sum, HCL.risk.ave, OB.risk.last, OB.risk.sum, OB.risk.ave)
 
cnt <- apply(xx[, cnt.nme], 2, sum, na.rm = TRUE)
lst <- apply(xx[, last.nme], 2, mean, na.rm = TRUE)
mn <- apply(xx[, median.nme], 2, mean, na.rm = TRUE)
cn <- cbind.data.frame(No.obs = n, t(cnt), t(mn), t(lst))
cbind(risk, com, cn)  
}


cnt.nme <- cont.vars[grep("_n", cont.vars)]
median.nme <- cont.vars[grep("_median", cont.vars)]
last.nme <- cont.vars[grep("_last", cont.vars)]
 
############################################################################
demo.vars <- c("age", "age.Jan2011", "sex","main_race","Education")
seq.vars <- c("DM.risk", "HTN.risk", "HCL.risk", "OB.risk")
risk.vars <-c("DM.risk.last", "DM.risk.sum", "DM.risk.ave", "HTN.risk.last", "HTN.risk.sum", "HTN.risk.ave", 
             "HCL.risk.last", "HCL.risk.sum", "HCL.risk.ave", "OB.risk.last", "OB.risk.sum", "OB.risk.ave")
com.vars <- c("DM.last", "DM.freq", "DM.any", "HTN.last", "HTN.any", "HCL.last", "HCL.freq", 
               "HCL.any", "Obese.last", "Obese.freq", "Obese.any")
  
resp.vars = c("dead", "mortality", "CVD") 
cont.vars <- c(cnt.nme, median.nme, last.nme)

vars <- c("pid", "time", resp.vars, demo.vars, seq.vars, com.vars, resp.vars, cont.vars) 
          
dd.mort <- ddply(dat.mort, .variables = "pid", .fun = Aggregated.Data, .parallel = TRUE)
dat.mort2 <- join(demo, dd.mort, by = "pid", type = "inner")


dd.cvd <- ddply(dat.cvd, .variables = "pid", .fun = Aggregated.Data, .parallel = TRUE)
dat.cvd2 <- join(demo, dd.cvd, by = "pid", type = "inner")

#filename = paste0(c(dir, "SequentialData_", paste0(c(k, "_months_interval_", as.character(today())), collapse = ""), ".csv"), collapse = "")
#write.csv(dat, file = filename, row.names = FALSE)

filename = paste0(c(dir, "SequentialData_", paste0(c(k, "_months_interval_", as.character(today())), collapse = ""), ".RData"), collapse = "")
save(dat.mort, dat.mort2, dat.cvd, dat.cvd2, demo,  demo.vars,seq.vars, com.vars, cont.vars, risk.vars,resp.vars, file = filename)


































