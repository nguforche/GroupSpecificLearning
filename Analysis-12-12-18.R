suppressMessages(require(plyr))
suppressMessages(require(caret))
suppressMessages(require(PresenceAbsence))
suppressMessages(require(missForest))
 
suppressMessages(library(parallel))
suppressMessages(library(doParallel))



library("devtools")
library(roxygen2)

Rcode <- as.package("/data5/bsi/clip/s201827.clip_hpc/che/MLpipeline/MLpipeline")
load_all(Rcode)
document(Rcode) 

dir = "/data5/bsi/clip/s201827.clip_hpc/che/GroupSpecificLearning/Data/"

k = 12 
filename = paste0(c(dir, "SequentialData_", paste0(c(k, "_months_interval_2018-12-12"), collapse = ""), ".RData"), collapse = "")

load(file = filename) 

models <- c("glm", "gbm", "rf")[1]

vars.last <- names(dat.mort2)[grep(".risk.last", names(dat.mort2))]
vars.sum <- names(dat.mort2)[grep(".risk.sum", names(dat.mort2))]
vars.ave <- names(dat.mort2)[grep(".risk.ave", names(dat.mort2))]
seqvars = c(vars.last,vars.sum,vars.ave)

com.last <- c("DM.last", "HTN.last", "HCL.last", "Obese.last")
com.freq <- names(dat.mort2)[grep(".freq", names(dat.mort2))]
com.any <- names(dat.mort2)[grep(".any", names(dat.mort2))]
com.all <- c(com.last,com.freq, com.any)

dd <- dat.mort2[sample(nrow(dat.mort2), nrow(dat.mort2)), ]
dd <- subset(dd, age >= 18) 

pid <- sample(nrow(dd), 10000)
dd <- dd[pid, ]

vars <- c(resp.vars, seqvars, com.all) 
ix.na <- sapply(dd[, vars], function(xx) sum(is.na(xx)) > 0)

cc <- complete.cases(dd[, c(vars)])
dd <- dd[cc, ] 


dd$mortality <- factor(ifelse(dd$mortality == 1, "Yes", "No"))


form <- as.formula(paste0("mortality ~ ", paste0(seqvars, collapse = "+"))) 
mod <- trainMLpipeline(form=form, data=dd, method = models,  importance=FALSE)
tab <- do.call(rbind, lapply(models, function(xx) {
perf <- Performance(object=mod, method= xx, positive.class = "Yes")
perf$Resample <- NULL 
cbind(model = xx, getSummary(perf,  groups = NULL,  alpha = 0.05))
}))
tab 


form1 <- as.formula(paste0("mortality ~ ", paste0(com.all, collapse = "+"))) 
mod1 <- trainMLpipeline(form=form1, data=dd, method = models, sampling = "down", IMR = 10,importance=FALSE)
tab1 <- do.call(rbind, lapply(models, function(xx) {
perf1 <- Performance(object=mod1, method= xx, positive.class = "Yes")
perf1$Resample <- NULL 
cbind(model = xx, getSummary(perf1,  groups = NULL,  alpha = 0.05))
}))
tab1 


#### cvd events 


d1 <- dat.cvd2[sample(nrow(dat.cvd2), nrow(dat.cvd2)), ]
d1 <- subset(d1, age >= 18) 

pid <- sample(nrow(d1), 10000)
d1 <- d1[pid, ]

cc <- complete.cases(d1[, c(vars)])
dd1 <- d1[cc, ] 

dd1$CVD <- factor(ifelse(dd1$CVD == 1, "Yes", "No"))

form <- as.formula(paste0("CVD ~ ", paste0(seqvars, collapse = "+"))) 
mod2 <- trainMLpipeline(form=form, data=dd1, method = models, importance=FALSE)
tab2 <- do.call(rbind, lapply(models, function(xx) {
perf <- Performance(object=mod2, method= xx, positive.class = "Yes")
perf$Resample <- NULL 
cbind(model = xx, getSummary(perf,  groups = NULL,  alpha = 0.05))
}))
tab2


form <- as.formula(paste0("CVD ~ ", paste0(com.all, collapse = "+"))) 
mod3 <- trainMLpipeline(form=form, data=dd1, method = models, importance=FALSE)
tab3 <- do.call(rbind, lapply(models, function(xx) {
perf <- Performance(object=mod3, method= xx, positive.class = "Yes")
perf$Resample <- NULL 
cbind(model = xx, getSummary(perf,  groups = NULL,  alpha = 0.05))
}))
tab3






form1 <- as.formula(paste0("mortality ~ ", paste0(com.any, collapse = "+"))) 
mod1 <- trainMLpipeline(form=form1, data=dd, method = models, sampling = "down", IMR = 10,importance=FALSE)

tab1 <- do.call(rbind, lapply(models, function(xx) {
perf1 <- Performance(object=mod1, method= xx, positive.class = "Yes")
perf1$Resample <- NULL 
cbind(model = xx, getSummary(perf1,  groups = NULL,  alpha = 0.05))
}))
tab1 

form2 <- as.formula(paste0("mortality ~ ", paste0(com.last, collapse = "+"))) 
mod2 <- trainMLpipeline(form=form2, data=dd, method = models, sampling = "down", IMR = 10,importance=FALSE)
tab2 <- do.call(rbind, lapply(models, function(xx) {
perf2 <- Performance(object=mod2, method= xx, positive.class = "Yes")
perf2$Resample <- NULL 
cbind(model = xx, getSummary(perf2,  groups = NULL,  alpha = 0.05))
}))
tab2


tab1 <- do.call(rbind, lapply(models, function(xx) {
perf1 <- Performance(object=mod1, method= xx, positive.class = "Yes")
perf1$Resample <- NULL 
cbind(model = xx, getSummary(perf1,  groups = NULL,  alpha = 0.05))
}))
tab1 


tab2 <- do.call(rbind, lapply(models, function(xx) {
perf2 <- Performance(object=mod2, method= xx, positive.class = "Yes")
perf2$Resample <- NULL 
cbind(model = xx, getSummary(perf2,  groups = NULL,  alpha = 0.05))
}))
tab2


write.table(tab, file = "seq.csv", sep = ",", row.names = FALSE) 
write.table(tab1, file = "com_last.csv", sep = ",", row.names = FALSE) 
write.table(tab2, file = "com_any.csv", sep = ",", row.names = FALSE) 


require(xtable)
require(Hmisc)
require(arsenal)

nme <- c("model", "PCC", "AUC", "sensitivity", "specificity", "G.mean", "BER", "F1.score", "Pos Pred Value")
print(xtable(tab[, nme]), include.rownames = FALSE)




dd <- dat.mort2[sample(nrow(dat.mort2), nrow(dat.mort2)), ]
dd <- subset(dd, age >= 18) 
cc <- complete.cases(dd[, c("mortality", seqvars,com.any,com.freq,com.last)])
dd <- dd[cc, ] 

form <- as.formula(paste0("mortality ~ ", paste0(c(demo.vars[c(1,3)], seqvars, com.any,com.freq, com.last), collapse = "+"))) 
dd[, c("sex", "DM.risk.last", "HTN.risk.last", "HCL.risk.last", "OB.risk.last")] <- lapply(dd[, c("sex", "DM.risk.last", "HTN.risk.last", "HCL.risk.last", "OB.risk.last")], as.factor)
dd[, c(com.any, com.last, com.freq)] <- lapply(dd[, c(com.any, com.last, com.freq)], as.factor)

val <-  tableby(formula = form, data = dd, cat.simplify = TRUE, numeric.stats = "meansd", numeric.simplify = TRUE, 
              cat.stats = "countpct", stats.labels = list(meansd = "Mean (SD)", countpct = "Count (Pct)"), digits = 2, digits.p= 3, )


as.data.frame(val)

print(summary(val), format = "latex")




res  <- summary(form1,  method="reverse", data = dd, quant=0.5, test=TRUE, overall=FALSE, continuous = 10)
  
val <- HmscData(x = res, digits=3)


, prn = any(n != N), pctdig = 0, what = c("%", 
    "proportion"), npct = c("numerator", "both", "denominator", 
    "none"), exclude1 = TRUE, vnames = c("labels", "names"), 
    prUnits = TRUE, sep = "/", abbreviate.dimnames = FALSE, prefix.width = max(nchar(lab)), 
    min.colwidth, formatArgs = NULL, round = NULL, prtest = c("P", 
        "stat", "df", "name"), prmsd = FALSE, long = FALSE, pdig = 3, 
    eps = 0.001, ...)  
  
  
  ltx <- latex(res, file = "periHmiscSummary.tex", digits = 2, npct = "none",vnames="name",  outer.size = "tiny", 
               prUnits = FALSE,middle.bold = FALSE, landscape=FALSE, insert.bottom = TRUE)

















