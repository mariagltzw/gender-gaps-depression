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

formula.outcome <- c("outcome ~ 
                     
                        ethn + education + RAFEDUC_cat + RAMEDUC_cat +
                       
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



# RUN G-FORMULA IN PARALLEL -----------------------------------------------

## LBRF INTERVENTION (Intervention A) ####
cl <- makeCluster(15) #use 15 cores
registerDoParallel(cl) 

set.seed(1234)

start.time <- Sys.time()

## according to stability test, 60 iterations for both MC and BS are sufficient for preliminary results
bs_it = 499 #bootstrap iterations
mc_it = 100 #monte carlo iterations

results.LBRF.int <- foreach(bs = 1:bs_it, .inorder=FALSE, .errorhandling="pass", .packages=c("cfdecomp","nnet", "Hmisc", "splines","quantreg","dplyr")) %dorng% {
  long.gformula_int(dataframe=df_gform, id="HHIDPN",formula.outcome=formula.outcome, time.var_names=time.var_names,
                    mcsize=mc_it, tunits=length(unique(df_gform$t)), intervention="LBRF", 
                    quant=c(0.25,0.5,0.75,0.90,0.95)) # default mcsize is 50
}

stopCluster(cl)

end.time <- Sys.time()
end.time - start.time # because of all the model fitting one bs iteration is expected to take around 15min


## save so we have the MC iterations as well
#save(results.LBRF.int, file=paste0("output_MC_LBRF.int_BS.",bs_it,"_MC.",mc_it,"_",Sys.Date(),".Rdata"))



# AGGREGATE BOOTSTRAP RESULTS ---------------------------------------------
tunits=length(unique(df_gform$t))

# make empty list to fill
bs.results.F <- bs.results.F.2 <- bs.results.F.3 <- bs.results.M <- bs.results.M.2 <- bs.results.M.3 <- list()

for (i in 1:length(time.var_names)){
  if(length(levels(df_gform[,time.var_names[i]])) > 2){
    matrix  <- matrix(NA,(tunits*length(levels(df_gform[,time.var_names[i]]))),3+bs_it)   
    rownames(matrix) <- rep(levels(df_gform[,time.var_names[i]]),tunits)
    matrix[,1] <- rep(seq(1,tunits,1),each=length(levels(df_gform[,time.var_names[i]])))
  } else {
    matrix <- matrix(NA,tunits,3+bs_it)
    matrix[,1] <- seq(1,tunits,1)
  }
  
  matrix[,2] <- 0 #female
  matrix[,3] <- 0 #nc
  bs.results.F[[i]]  <- matrix
  
  matrix[,3] <- 1 #cf - TE
  bs.results.F.2[[i]]  <- matrix

  matrix[,2] <- 1 #male
  matrix[,3] <- 0 #nc
  bs.results.M[[i]]  <- matrix
  
  matrix[,3] <- 1 #cf - TE
  bs.results.M.2[[i]]  <- matrix

}

names(bs.results.F) <- names(bs.results.M) <- names(bs.results.F.2) <- names(bs.results.M.2) <- time.var_names

#aggregate
for (i in 1:length(time.var_names)){
  for (b in 1:bs_it){
    bs.results.F[[i]][,3+b] <- rowMeans(results.LBRF.int[[b]][[1]][[i]][,2:mc_it])
    bs.results.M[[i]][,3+b] <- rowMeans(results.LBRF.int[[b]][[2]][[i]][,2:mc_it])
    
    bs.results.F.2[[i]][,3+b] <- rowMeans(results.LBRF.int[[b]][[3]][[i]][,2:mc_it])
    bs.results.M.2[[i]][,3+b] <- rowMeans(results.LBRF.int[[b]][[4]][[i]][,2:mc_it])
  }}

results_output <- list(bs.results.F, bs.results.M,
                       bs.results.F.2, bs.results.M.2)

names(results_output) <- c("natural.course_females","natural.course_males",
                           "TE_int_females","TE_int_males")

# save bootstrap aggregated results
save(results_output, file=paste0("output_LBRF.int_BS",bs_it,"_MC",mc_it,"_",Sys.Date(),".Rdata"))

rm(results_output);rm(results.LBRF.int)
gc()



#### REPEAT G-FORMULA RUN FOR THE OTHER THREE INTERVENTIONS ##############################
#### THE FUNCTION INPUT IS THE SAME ONLY THE INTERVENTION= INPUT NEEDS TO BE CHANGED #####
#### FOR LBRF and OCCUPATION: intervention = "LBRF+OCC" ################################## 
#### FOR LBRF, OCCUPATION and income: intervention = "LBRF+OCC+INCOME" ###################
#### FOR LBRF, OCCUPATION, income and previous income: intervention = "ALL + PAST INCOME"#
