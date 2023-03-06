library(dplyr)
library(ggplot2)
library(cfdecomp)
library(foreach)
library("doParallel")
library(doRNG)
library(nnet)
library("foreign")
library("lattice")
library(Hmisc)
library(MASS)
library(splines)
library(quantreg)
library(mice) # multiple imputation with chained equations
library(nnet)

load("df_inflation.adj.RData")
df_gform <- df_inflation.adj

#load functions
source("functions.R")

##load the dataset in which we imputed missing data for t1 (ages 50-51)----
load("df_gform.t1_imp.RData")

# PREPARE INPUT FOR LONGITUDINAL G-FORMULA --------------------------------

# formula for regression models
# we use output as Y here, we will replace it with the time-varying covariate of interest
#I EXCLUDED ETHN!
formula.outcome <- c("outcome ~ 
                     
                        education + RAFEDUC_cat + RAMEDUC_cat +
                       
                       psych_probl_E + l.marital_stat_cat2 + l.H_members + l.H_child_cat +
                       
                         l.depress + relevel(l.occupation_class_4cat, ref = '4.doesnt work') + 
                         
                          relevel(l.LBRF, ref = '7.Not in LbrF') * (ns(AGEY_E,3) + l.ind_earn_adj) + 
                         
                         l.N_cond_cat ")
                         

# get column index for non lagged column names
names.from <- c("depress", "H_members", "H_child_cat", "N_cond_cat","marital_stat_cat2", "LBRF", "ind_earn_adj","occupation_class_4cat")

#get the lagged variables
names.to <- paste0("l.", names.from)

# assign as time varying covariate names
time.var_names <- names.from

#get their col.index
col.index.from <- unlist(lapply(names.from, col.index, df_gform), use.names = FALSE)
col.index.to <- unlist(lapply(names.to, col.index, df_gform), use.names = FALSE)

#before we start the analysis we need to exclude some observations due to scarcity issues at ages 50-51
# exclude
# . non-Hisp Black 4.Partly retired 1.White collar/desk job N=7
# 3.Hispanic 4.Partly retired 1.White collar/desk job N=4
# find ids that have this combination at t=1 and exclude them from dataframe=df_gform; dataframe.t1=df_gform.t1_imp
 table <- df_gform.t1_imp %>% group_by(ethn,LBRF,occupation_class_4cat) %>%
   summarise(ind_earn_adj = mean(ind_earn_adj[RAGENDER=="male"]),
             n=n())

rows.to.exclude <- which(df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$ethn =="2. non-Hisp Black"
                         & df_gform.t1_imp$LBRF == "4.Partly retired" & df_gform.t1_imp$occupation_class_4cat =="2.Pink collar/Service-related job"|
                           df_gform.t1_imp$ethn =="3.Hispanic"
                         & df_gform.t1_imp$LBRF == "4.Partly retired" & df_gform.t1_imp$occupation_class_4cat =="1.White collar/desk job"|
                           df_gform.t1_imp$ethn =="4.Other"
                         & df_gform.t1_imp$LBRF == "2.Works PT" & df_gform.t1_imp$occupation_class_4cat =="3.Blue-collar/Manual job"|
                           df_gform.t1_imp$ethn =="4.Other"
                         & df_gform.t1_imp$LBRF == "4.Partly retired" & df_gform.t1_imp$occupation_class_4cat =="2.Pink collar/Service-related job")
#observations excluded: length(rows.to.exclude)
ids.to.exclude <- df_gform.t1_imp$HHIDPN[rows.to.exclude]

df_gform.t1_imp <- subset(df_gform.t1_imp, df_gform.t1_imp$HHIDPN %in% ids.to.exclude==FALSE)

#same for the other sample
df_gform <- subset(df_gform, df_gform$HHIDPN %in% ids.to.exclude==FALSE)


# RUN G-FORMULA IN PARALLEL -----------------------------------------------

## LBRF INTERVENTION (Intervention A) ####
cl <- makeCluster(15) #use 15 cores
registerDoParallel(cl) 

set.seed(1234)

start.time <- Sys.time()

## according to stability test, 60 iterations for both MC and BS are sufficient for preliminary results
## the specification below are to reproduce the results of the paper
bs_it = 499 #bootstrap iterations
mc_it = 100 #monte carlo iterations

results.LBRF.ethn <- foreach(bs = 1:bs_it, .inorder=FALSE, .errorhandling="pass",.packages=c("cfdecomp","nnet", "Hmisc", "splines","quantreg","dplyr")) %dorng% {
  long.gformula_ethn(dataframe=df_gform, id="HHIDPN",formula.outcome=formula.outcome, time.var_names=time.var_names,
                mcsize=mc_it, tunits=length(unique(df_gform$t)), intervention="LBRF",
                quant=c(0.25,0.5,0.75,0.90,0.95))#default mcsize is 50
}

stopCluster(cl)

end.time <- Sys.time()
end.time - start.time # because of all the model fitting one bs iteration is expected to take around 45min

## save so we have the MC iterations as well
#save(results.LBRF.ethn, file=paste0("output_ethn.subgroup_MC_LBRF.int_BS.",bs_it,"_MC.",mc_it,"_",Sys.Date(),".Rdata"))

# AGGREGATE BOOTSTRAP RESULTS ---------------------------------------------
tunits=length(unique(df_gform$t))

#make list to fill
bs.results.F.w <- bs.results.F.b <- bs.results.F.h <- list()
bs.results.F.2.w <- bs.results.F.2.b <- bs.results.F.2.h <- list()

bs.results.M.w <- bs.results.M.b <- bs.results.M.h <- list()
bs.results.M.2.w <- bs.results.M.2.b <- bs.results.M.2.h <- list()
  
  
for (i in 1:length(time.var_names)){
  if(length(levels(df_gform[,time.var_names[i]])) > 2){
    matrix  <- matrix(NA,(tunits*length(levels(df_gform[,time.var_names[i]]))),4+bs_it)   
    rownames(matrix) <- rep(levels(df_gform[,time.var_names[i]]),tunits)
    matrix[,1] <- rep(seq(1,tunits,1),each=length(levels(df_gform[,time.var_names[i]])))
  } else {
    matrix <- matrix(NA,tunits,4+bs_it)
    matrix[,1] <- seq(1,tunits,1)
  }
  
  ## females
  matrix[,2] <- 0 #female
  ## nc
  matrix[,3] <- 0 #nc
  matrix[,4] <- 1 #white
  bs.results.F.w[[i]]  <- matrix
  
  matrix[,4] <- 2 #hisp
  bs.results.F.b[[i]]  <- matrix
  
  matrix[,4] <- 3 #black
  bs.results.F.h[[i]]  <- matrix

  
  ##cf
  matrix[,3] <- 1 #cf
  matrix[,4] <- 1 #white
  bs.results.F.2.w[[i]]  <- matrix
  
  matrix[,4] <- 2 #hisp
  bs.results.F.2.b[[i]]  <- matrix
  
  matrix[,4] <- 3 #black
  bs.results.F.2.h[[i]]  <- matrix
  
  
  ## males
  matrix[,2] <- 1 #male
  ##nc
  matrix[,3] <- 0 #nc
  matrix[,4] <- 1 #white
  bs.results.M.w[[i]]  <- matrix
  
  matrix[,4] <- 2 #hisp
  bs.results.M.b[[i]]  <- matrix
  
  matrix[,4] <- 3 #black
  bs.results.M.h[[i]]  <- matrix
  

  ##cf
  matrix[,3] <-1 #cf
  matrix[,4] <- 1 #white
  bs.results.M.2.w[[i]]  <- matrix
  
  matrix[,4] <- 2 #hisp
  bs.results.M.2.b[[i]]  <- matrix
  
  matrix[,4] <- 3 #black
  bs.results.M.2.h[[i]]  <- matrix
  
  
}

names(bs.results.F.w) <- names(bs.results.F.b) <- names(bs.results.F.h) <- time.var_names
names(bs.results.F.2.w) <- names(bs.results.F.2.b) <- names(bs.results.F.2.h) <- time.var_names

names(bs.results.M.w) <- names(bs.results.M.b) <- names(bs.results.M.h) <-time.var_names
names(bs.results.M.2.w) <- names(bs.results.M.2.b) <- names(bs.results.M.2.h) <-  time.var_names

for (i in 1:length(time.var_names)){
  for (b in 1:bs_it){
    bs.results.F.w[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[1]][[i]][,2:mc_it])
    bs.results.F.b[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[2]][[i]][,2:mc_it])
    bs.results.F.h[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[3]][[i]][,2:mc_it])

    bs.results.M.w[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[4]][[i]][,2:mc_it])
    bs.results.M.b[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[5]][[i]][,2:mc_it])
    bs.results.M.h[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[6]][[i]][,2:mc_it])

    bs.results.F.2.w[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[7]][[i]][,2:mc_it])
    bs.results.F.2.b[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[8]][[i]][,2:mc_it])
    bs.results.F.2.h[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[9]][[i]][,2:mc_it])

    bs.results.M.2.w[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[10]][[i]][,2:mc_it])
    bs.results.M.2.b[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[11]][[i]][,2:mc_it])
    bs.results.M.2.h[[i]][,4+b] <- rowMeans(results.LBRF.ethn[[b]][[12]][[i]][,2:mc_it])

  }}

results_output <- list(bs.results.F.w, bs.results.F.b, bs.results.F.h, 
                       bs.results.F.2.w, bs.results.F.2.b, bs.results.F.2.h, 
                       bs.results.M.w, bs.results.M.b, bs.results.M.h, 
                       bs.results.M.2.w, bs.results.M.2.b, bs.results.M.2.h)
names(results_output) <- c("natural.course_females_white","natural.course_females_black","natural.course_females_hisp",
                           "TE_int_females_white","TE_int_females_black","TE_int_females_hisp",
                           "natural.course_males_white","natural.course_males_black","natural.course_males_hisp",
                           "TE_int_males_white","TE_int_males_black","TE_int_males_hisp")

save(results_output, file=paste0("output_ethn.subgroup_LBRF.int_BS",bs_it,"_MC",mc_it,"_",Sys.Date(),".Rdata"))
rm(results_output);rm(results.LBRF.ethn)
gc()


#### REPEAT G-FORMULA RUN FOR THE OTHER THREE INTERVENTIONS ##############################
#### THE FUNCTION INPUT IS THE SAME ONLY THE INTERVENTION= INPUT NEEDS TO BE CHANGED #####
#### FOR LBRF and OCCUPATION: intervention = "LBRF+OCC" ################################## 
#### FOR LBRF, OCCUPATION and income: intervention = "LBRF+OCC+INCOME" ###################
#### FOR LBRF, OCCUPATION, income and previous income: intervention = "ALL + PAST INCOME"#
