#########################################
########### FUNCTIONS FOR     ########### 
########### GENDER GAPS PAPER ########### 
########### GUELTZOW ET AL    ########### 
#########################################

# FUNCTIONS NEEDED TO PREPARE FOR G-FORMULA -------------------------------

## gets column index ##
col.index <- function(names,df){
  return(grep(paste0('^',names,'$'), colnames(df)))
}

## lag and lead function ##
## lags variable by one wave ## 
lag1.func <- function(variabletolag, idvariable, timevariable) {
  
  lag <- c(NA,ifelse(idvariable[2:length(idvariable)] == idvariable[1:(length(idvariable)-1)] &
                       timevariable[2:length(timevariable)] - timevariable[1:(length(timevariable)-1)] ==1,
                     variabletolag[1:(length(idvariable)-1)],NA))
  if (class(variabletolag)=="factor"){
    lag <- as.factor(lag)
    levels(lag) <- levels(variabletolag)
  }
  
  lag
}

## leads variable by one wave ## 
lead1.func <- function(variabletolead, idvariable, timevariable) {
  
  lead <- c(ifelse(idvariable[1:length(idvariable)-1] == idvariable[2:(length(idvariable))] &
                     timevariable[1:length(timevariable)-1] - timevariable[2:(length(timevariable))] == -1,
                   variabletolead[2:(length(idvariable))],NA),NA)
  if (class(variabletolead)=="factor"){
    lead <- as.factor(lead)
    levels(lead) <- levels(variabletolead)
  }
  
  lead
}



# FUNCTIONS USED WITHIN THE LONGITUDINAL G-FORMULA ------------------------------------

## predicts and then draws random values from model predicted prob. ##
## can accomodate glm, rqs and multinom ##
predict_function <- function(model1, dataset, quantiles=NA,residuals=NA) {
    if(class(model1)[1] =="glm"){
    x1 <- predict(model1,dataset,type='response')
    x1c <- as.factor(rbinom(n=length(x1),1,prob=x1))
    return(x1c)
    
  } else if (class(model1)[1] =="rqs"){ #quantile regression
    
    boundary <- bound.func(quantiles)
    
    if (length(quantiles) == 6){
      x1c  <- quant.assign.6(dataset, boundary, model1)
    } else if (length(quantiles) == 5) {
      x1c  <- quant.assign.5(dataset, boundary, model1)
    } else if (length(quantiles) == 4) {
      x1c  <- quant.assign.4(dataset, boundary, model1)
    } else if (length(quantiles) == 3) {
      x1c  <- quant.assign.3(dataset, boundary, model1)
    }

    return(x1c)
    
  } else if (class(model1)[1] =="multinom"){
    x1  <- predict(model1, type="probs",newdata = dataset)
    x1c <- as.factor(rMultinom(m = 1, probs = x1))
    return(x1c)
  } else {
    return("Model is not glm or multinomial!")
  }
}

## an additional function called by the predict_function ##
## used for predicting multiple quantiles from quantile regression ## 
## determine boundary values ##
bound.func <- function(x) {
  boundary <- NULL
  for(i in 2:length(x)) {
    boundary[i-1] <- mean(c(x[i-1],x[i]))
  }
  return(boundary)
}

## predict based on 6, 5, 4, or 3 quantiles ##
quant.assign.6 <- function(data,boundary,
                           fit_allQ) {
  
  pred_all <- predict(fit_allQ,type='response',data)
  
  pred1 <- pred_all[,1]
  pred2 <- pred_all[,2]
  pred3 <- pred_all[,3]
  pred4 <- pred_all[,4]
  pred5 <- pred_all[,5]
  pred6 <- pred_all[,6]
  
  n <- length(pred1)
  r <- runif(n)
  
  output <- ifelse(r <= boundary[1],pred1,
                   ifelse(r <= boundary[2],pred2,
                          ifelse(r <= boundary[3],pred3,
                                 ifelse(r <= boundary[4],pred4,
                                        ifelse(r <= boundary[5],pred5,
                                        pred6)))))
  return(output)
}
quant.assign.5 <- function(data,boundary,
                              fit_allQ) {
  
  pred_all <- predict(fit_allQ,type='response',data)
  
  pred1 <- pred_all[,1]
  pred2 <- pred_all[,2]
  pred3 <- pred_all[,3]
  pred4 <- pred_all[,4]
  pred5 <- pred_all[,5]
  
  n <- length(pred1)
  r <- runif(n)
  
  output <- ifelse(r <= boundary[1],pred1,
                   ifelse(r <= boundary[2],pred2,
                          ifelse(r <= boundary[3],pred3,
                                 ifelse(r <= boundary[4],pred4,
                                        pred5))))
  return(output)
}
quant.assign.4 <- function(data,boundary,
                           fit_allQ) {
  
  pred_all <- predict(fit_allQ,type='response',data)
  
  pred1 <- pred_all[,1]
  pred2 <- pred_all[,2]
  pred3 <- pred_all[,3]
  pred4 <- pred_all[,4]

  n <- length(pred1)
  r <- runif(n)
  
  output <- ifelse(r <= boundary[1],pred1,
                   ifelse(r <= boundary[2],pred2,
                          ifelse(r <= boundary[3],pred3,pred4)))
  return(output)
}
quant.assign.3 <- function(data,boundary,
                           fit_allQ) {
  
  pred_all <- predict(fit_allQ,type='response',data)
  
  pred1 <- pred_all[,1]
  pred2 <- pred_all[,2]
  pred3 <- pred_all[,3]
  
  n <- length(pred1)
  r <- runif(n)
  
  output <- ifelse(r <= boundary[1],pred1,
                   ifelse(r <= boundary[2],pred2, pred3))
  return(output)
}



# LONGITUDINAL G FORMULA FUNCTION FOR WHOLE SAMPLE ------------------------

## runs longitudinal gformula ##
## quantile regression for income with default quantiles specified in function ##
## needs input like this (for example):
## dataframe=df;id="HHIDPN"; mcsize=10; tunits=16 ##

## formula.outcome is the model specification as as character string with "outcome" as y: "outcome ~ x1+x2+..." ##
## it will be turned into a formula later

## time.var_names is a character string of the variables that will be modeled as time varying ##

## !!!! the way it is coded right now, income needs to be the last variable name in time.var_names!!!!
## current function can deal with categorical, binary and zero-inflated income variable ##
## !!!! needs to be adapted if other cont variables should be included!!!!!
long.gformula_int <- function(dataframe,dataframe.t1=df_gform.t1_imp,id,formula.outcome, time.var_names, mcsize=50,
                              tunits, quant=c(0.25,0.5,0.75,0.90,0.95),intervention) { 
  
  ## sample
  # sample individuals (with replacement) from ex.dat
  re.size <- length(unique(dataframe[,id]))
  ex.sample.F <- cluster.resample(dataframe[dataframe$RAGENDER=="female",], cluster.name = id, 
                                  size = re.size)
  ex.sample.M <- cluster.resample(dataframe[dataframe$RAGENDER=="male",], cluster.name = id, 
                                  size = re.size) 
  
  #### IMPUTATION CODE - GLM 
  ## IMPUTE PARENTS EDUCATION AND OCCUPATION ##
  # RAFEDUC
  fit_RAFED_F <- multinom(RAFEDUC_cat ~ ethn + education + RAMEDUC_cat + psych_probl_E +
                            marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat + occupation_class_3cat +
                            LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                            l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                            l.occupation_class_3cat + 
                            l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F)
  
  fit_RAFED_M    <- update(fit_RAFED_F,data=ex.sample.M)
  
  # RAMEDUC
  fit_RAMED_F <- multinom(RAMEDUC_cat ~ ethn + education + RAFEDUC_cat  + psych_probl_E +
                            marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat + occupation_class_3cat +
                            LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                            l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                            l.occupation_class_3cat + 
                            l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F)
  
  fit_RAMED_M    <- update(fit_RAMED_F,data=ex.sample.M)
  
  # occupation
  fit_occ_F <- multinom(occupation_class_3cat ~ ethn + education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E +
                          marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat +
                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                          l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                          l.occupation_class_3cat + 
                          l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F)
  
  fit_occ_M    <- update(fit_occ_F,data=ex.sample.M)
  
  ## predict only for after t1 ##
  ## for t1 we don't need to predict anymore, we can just use the already imputed dataframe.t1
  
  #RAFEDUC
  ex.sample.F$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F$RAFEDUC_cat) & ex.sample.F$t==1,
                                              as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female"]),
                                              ifelse(is.na(ex.sample.F$RAFEDUC_cat) & ex.sample.F$t!=1, as.character(predict(fit_RAFED_F,type = "class")),
                                                     as.character(ex.sample.F$RAFEDUC_cat))))
  ex.sample.M$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M$RAFEDUC_cat) & ex.sample.M$t==1,
                                              as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male"]),
                                              ifelse(is.na(ex.sample.M$RAFEDUC_cat) & ex.sample.M$t!=1, as.character(predict(fit_RAFED_M,type = "class")),
                                                     as.character(ex.sample.M$RAFEDUC_cat))))
  #RAMEDUC
  ex.sample.F$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F$RAMEDUC_cat) & ex.sample.F$t==1,
                                              as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female"]),
                                              ifelse(is.na(ex.sample.F$RAMEDUC_cat) & ex.sample.F$t!=1, as.character(predict(fit_RAMED_F,type = "class")),
                                                     as.character(ex.sample.F$RAMEDUC_cat))))
  ex.sample.M$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M$RAMEDUC_cat) & ex.sample.M$t==1,
                                              as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male"]),
                                              ifelse(is.na(ex.sample.M$RAMEDUC_cat) & ex.sample.M$t!=1, as.character(predict(fit_RAMED_M,type = "class")),
                                                     as.character(ex.sample.M$RAMEDUC_cat))))
  #occupation class 
  ex.sample.F$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F$occupation_class_4cat) & ex.sample.F$t==1,
                                                        as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female"]),
                                                        ifelse(is.na(ex.sample.F$occupation_class_4cat) & ex.sample.F$t!=1, as.character(predict(fit_occ_F,type = "class")),
                                                               as.character(ex.sample.F$occupation_class_4cat))))
  ex.sample.M$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M$occupation_class_4cat) & ex.sample.M$t==1,
                                                        as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male"]),
                                                        ifelse(is.na(ex.sample.M$occupation_class_4cat) & ex.sample.M$t!=1, as.character(predict(fit_occ_M,type = "class")),
                                                               as.character(ex.sample.M$occupation_class_4cat))))
  
  ## remake 3 class occ and l.occupation ##
  ex.sample.F$occupation_class_3cat <- as.factor(ifelse(ex.sample.F$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F$occupation_class_4cat)))
  ex.sample.M$occupation_class_3cat <- as.factor(ifelse(ex.sample.M$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M$occupation_class_4cat)))
  
  ex.sample.F$l.occupation_class_4cat <- lag1.func(ex.sample.F$occupation_class_4cat, ex.sample.F$HHIDPN, ex.sample.F$wave)
  ex.sample.M$l.occupation_class_4cat <- lag1.func(ex.sample.M$occupation_class_4cat, ex.sample.M$HHIDPN, ex.sample.M$wave)
  
  #### IMPUTATION CODE - GLM
  
  ex.sample.F.3 <- ex.sample.F.2 <- ex.sample.F
  ex.sample.M.3 <- ex.sample.M.2 <- ex.sample.M
  
  ## (re)fit models to ex.sample
  models_females <- vector(mode = "list", length = length(time.var_names))
  models_males   <- vector(mode = "list", length = length(time.var_names))
  
  for (i in 1:length(time.var_names)){
    
    formula <- as.formula(sub('outcome',paste0(time.var_names[i]),formula.outcome))
    
    if (class(dataframe[,time.var_names[i]]) == "numeric") {
      
      #multiple quantiles
      
      ## we will use quantile regression here because the distribution of income looks pretty bad (zero inflated and outliers)
      ## For larger problems it is advantageous to use the Frisch-Newton interior point method "fn"
      
      models_females[[i]] <- rq(formula, data=ex.sample.F, tau=quant, method="fn") #replace quant with 0.5 for median regression only
      models_males[[i]]   <- rq(formula, data=ex.sample.M, tau=quant, method="fn") #replace quant with 0.5 for median regression only
      
    } else if (length(unique(dataframe[,time.var_names[i]])) == 2){
      ## for binary ##
      models_females[[i]] <- glm(formula, family=binomial, data=ex.sample.F)
      models_males[[i]]   <- glm(formula, family=binomial, data=ex.sample.M)
      
    } else if (class(dataframe[,time.var_names[i]]) == "factor" & length(unique(dataframe[,time.var_names[i]])) > 2){         
      ## for categorical ##
      models_females[[i]] <- multinom(formula, data=ex.sample.F)
      models_males[[i]]   <- multinom(formula, data=ex.sample.M)
      
    } else if (time.var_names == "occupation_class"){  
      
      ##occupation class is a spatial case because the "doesnt work" category is deterministically based on LBRF
      models_females[[i]] <- multinom(occupation_class_3cat ~ ethn + education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.F)
      models_males[[i]] <- multinom(occupation_class_3cat ~ ethn + education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                      l.H_child_cat + l.depress +
                                      LBRF * (ns(AGEY_E,3) + ind_earn_adj)  + l.N_cond_cat,
                                    data=ex.sample.M)
      
    } else {print("Time-varying variable has no unique levels/values!")} 
    ## currently no implementation for continuous variable other than quantile regression ##
  }
  
  ## save models ##
  names(models_females) <- lapply(time.var_names, FUN = paste0,".model.F")
  names(models_males)   <- lapply(time.var_names, FUN = paste0,".model.M")
  
  ## create list made of empty matrices to save results ##
  mc.results.arr.F  <- list()
  
  for (i in 1:length(time.var_names)){
    if (length(levels(dataframe[,time.var_names[i]])) > 2){
      matrix  <- matrix(NA,(tunits*length(levels(dataframe[,time.var_names[i]]))),1+mcsize)   
      matrix[,1] <- rep(seq(1,tunits,1),each=length(levels(dataframe[,time.var_names[i]])))
      rownames(matrix) <- rep(levels(dataframe[,time.var_names[i]]),tunits)
      mc.results.arr.F[[i]] <- matrix
    } else {
      matrix <- matrix(NA,tunits,1+mcsize)
      matrix[,1] <- seq(1,tunits,1)
      mc.results.arr.F[[i]]  <- matrix
    }
  }
  names(mc.results.arr.F) <- time.var_names
  
  mc.results.arr.F.2 <- mc.results.arr.F.3 <- mc.results.arr.F
  mc.results.arr.M.2 <- mc.results.arr.M.3 <- mc.results.arr.M <- mc.results.arr.F
  
  ## start Monte Carlo Error Reduction Loop ## 
  for(m in 1:mcsize) {


    
    ## take individuals at time 1 and discard the other observations ##
    
    ######## start simulation with imputation at age 50 
    ex.sample.F <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female",]
    ex.sample.M <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male",]
    
    ## give idnr
    ex.sample.F$idnr <- 1:length(ex.sample.F[,id])
    ex.sample.M$idnr <- 1:length(ex.sample.M[,id])
    
    ## make copies for intervention and mediation simulations
    ex.sample.F.3 <- ex.sample.F.2 <- ex.sample.F
    ex.sample.M.3 <- ex.sample.M.2 <- ex.sample.M
    
    if (intervention=="ALL + PAST INCOME"){
      ## preparation for l.income intervention (intervention D)
      ## take income levels at t1 (age 50, 51) in males 
      ## and calculate mean levels based on lbrf status and occupational group
      ## (if we use more grouping variables we need to consider using a regression model
      ## because a prop table will have scarcity issues)
      ## and then give females new income levels at t1 based on that prop.table
      
      df_calc_mean.inc <- rbind(ex.sample.F.2, ex.sample.M.2)
      df_calc_mean.inc <- df_calc_mean.inc %>% group_by(LBRF, occupation_class_4cat) %>% mutate(ind_earn_adj = mean(ind_earn_adj[RAGENDER=="male"])) %>% 
        ungroup()
      
      ex.sample.F.2 <-  df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female",]   
      ex.sample.F.2 <- as.data.frame(ex.sample.F.2)
    }
    
    

# NATURAL COURSE APPROXIMATION LOOP ---------------------------------------
    ## start a loop that moves through the follow-up time units
    ## this part of the g-formula tries to reproduce the empirical data
    ## and is known as the 'natural course'
    
    ## for each time unit ##
    for(t in 2:tunits) {
      
      ## female
      ## make a copy of the previous row and then update year and age
      ## but don't make a copy if someone was censored in the previous year
      ex.sample.temp.F <- ex.sample.F[ex.sample.F$t==(t-1),]
      ex.sample.temp.F$t <- ex.sample.temp.F$t+1
      ex.sample.temp.F$AGEY_E <- ex.sample.temp.F$AGEY_E+2
      
      ## lag values
      ex.sample.temp.F[,col.index.to] <- ex.sample.temp.F[,col.index.from]
      
      ## males
      ex.sample.temp.M <- ex.sample.M[ex.sample.M$t==(t-1),]
      ex.sample.temp.M$t <- ex.sample.temp.M$t+1
      ex.sample.temp.M$AGEY_E <- ex.sample.temp.M$AGEY_E+2
      
      ## lag values
      ex.sample.temp.M[,col.index.to] <- ex.sample.temp.M[,col.index.from]
      
      ## predict time varying variables ##
      ## occupation is predicted slightly differently ##
      ## so we can't include it in this coming for loop ##
      l <- length(time.var_names)-2
      
      for (i in 1:l){
        ex.sample.temp.F[,time.var_names[i]] <- predict_function(model1=models_females[[i]], dataset=ex.sample.temp.F)
        ex.sample.temp.M[,time.var_names[i]] <- predict_function(model1=models_males[[i]], dataset=ex.sample.temp.M)
      }
      
      ex.sample.temp.F$ind_earn_adj <- predict_function(model1=models_females$ind_earn_adj, dataset=ex.sample.temp.F,
                                                        quantiles=quant)
      ex.sample.temp.M$ind_earn_adj <- predict_function(model1=models_males$ind_earn_adj, dataset=ex.sample.temp.M,
                                                        quantiles=quant)
      
      
      ## predict occupation ##
      ex.sample.temp.F$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F$LBRF == "3.Unemployed"|
                                                                   ex.sample.temp.F$LBRF == "5.Retired"|
                                                                   ex.sample.temp.F$LBRF == "6.Disabled"|
                                                                   ex.sample.temp.F$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                 as.character(predict_function(model1=models_females$occupation_class_4cat.model.F, dataset=ex.sample.temp.F))))
      ex.sample.temp.M$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M$LBRF == "3.Unemployed"|
                                                                   ex.sample.temp.M$LBRF == "5.Retired"|
                                                                   ex.sample.temp.M$LBRF == "6.Disabled"|
                                                                   ex.sample.temp.M$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                 as.character(predict_function(model1=models_males$occupation_class_4cat.model.M, dataset=ex.sample.temp.M))))
      
      ## join this newly created data with the old data ##
      ex.sample.F <- rbind(ex.sample.F,ex.sample.temp.F)
      rm(ex.sample.temp.F)
      
      ## order by ID variable - not necessary, but make sure it doesnt mess anything up!!! ##
      ex.sample.F <- ex.sample.F[order(ex.sample.F$idnr,ex.sample.F$t),]
      
      ## males
      ex.sample.M <- rbind(ex.sample.M,ex.sample.temp.M)
      rm(ex.sample.temp.M)
      
      ## order by ID variable - not necessary, but make sure it doesnt mess anything up!!! ##
      ex.sample.M <- ex.sample.M[order(ex.sample.M$idnr,ex.sample.M$t),]
    }
    

# INTERVENTION LOOP -------------------------------------------------------
    ## now let's create a counterfactual dataset for the TOTAL EFFECT ##
    ## we will equalise LBRF distribution among males and females
    ## assuming that females have the same LBRF distribution as males
    
    ## for each time unit ##
    for(t in 2:tunits) {
      
      ## make a copy of the previous row and then update year and age
      ## female
      ex.sample.temp.F.2 <- ex.sample.F.2[ex.sample.F.2$t==(t-1),]
      ex.sample.temp.F.2$t <- ex.sample.temp.F.2$t+1
      ex.sample.temp.F.2$AGEY_E <- ex.sample.temp.F.2$AGEY_E+2
      
      ## lag values
      ex.sample.temp.F.2[,col.index.to] <- ex.sample.temp.F.2[,col.index.from]
      
      ## males
      ex.sample.temp.M.2 <- ex.sample.M.2[ex.sample.M.2$t==(t-1),]
      ex.sample.temp.M.2$t <- ex.sample.temp.M.2$t+1
      ex.sample.temp.M.2$AGEY_E <- ex.sample.temp.M.2$AGEY_E+2
      
      ## lag values
      ex.sample.temp.M.2[,col.index.to] <- ex.sample.temp.M.2[,col.index.from]
      
      ## predict time varying variables ##
      ## occupation is predicted slightly differently ##
      ## so we can't include it in this coming for loop ##
      l <- length(time.var_names)-2
      for (i in 1:l){
        ex.sample.temp.F.2[,time.var_names[i]] <- predict_function(model1=models_females[[i]], dataset=ex.sample.temp.F.2)
        ex.sample.temp.M.2[,time.var_names[i]] <- predict_function(model1=models_males[[i]], dataset=ex.sample.temp.M.2)
      }
      
      if (intervention=="LBRF"){ ## equalize employment opportunities

        ex.sample.temp.F.2$ind_earn_adj <- predict_function(model1=models_females$ind_earn_adj, dataset=ex.sample.temp.F.2,
                                                            quantiles=quant)
        ex.sample.temp.M.2$ind_earn_adj <- predict_function(model1=models_males$ind_earn_adj, dataset=ex.sample.temp.M.2,
                                                            quantiles=quant)
        
      ## intervene here so that l.LBRF is changed ##
      ## take the observed distribution of employment in males and randomly assign it to females ##
        ex.sample.temp.F.2$LBRF <- predict_function(models_males$LBRF.model.M, ex.sample.temp.F.2)
        
      ## predict occupation ##
        ex.sample.temp.F.2$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2$LBRF == "3.Unemployed"|
                                                                       ex.sample.temp.F.2$LBRF == "5.Retired"|
                                                                       ex.sample.temp.F.2$LBRF == "6.Disabled"|
                                                                       ex.sample.temp.F.2$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                     as.character(predict_function(model1=models_females$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2))))
        ex.sample.temp.M.2$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2$LBRF == "3.Unemployed"|
                                                                       ex.sample.temp.M.2$LBRF == "5.Retired"|
                                                                       ex.sample.temp.M.2$LBRF == "6.Disabled"|
                                                                       ex.sample.temp.M.2$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                     as.character(predict_function(model1=models_males$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2))))
        
      } else if (intervention=="LBRF+OCC"){ ## equalize employment and occupation opportunities
        
        ex.sample.temp.F.2$ind_earn_adj <- predict_function(model1=models_females$ind_earn_adj, dataset=ex.sample.temp.F.2,
                                                            quantiles=quant)
        ex.sample.temp.M.2$ind_earn_adj <- predict_function(model1=models_males$ind_earn_adj, dataset=ex.sample.temp.M.2,
                                                            quantiles=quant)
        
        ## intervene here so that l.LBRF is changed ##
        ## take the observed distribution of employment in males and randomly assign it to females ##
        ex.sample.temp.F.2$LBRF <- predict_function(models_males$LBRF.model.M, ex.sample.temp.F.2)
        
        ## predict occupation ##
        ## intervene on occupation ##
        ex.sample.temp.F.2$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2$LBRF == "3.Unemployed"|
                                                                       ex.sample.temp.F.2$LBRF == "5.Retired"|
                                                                       ex.sample.temp.F.2$LBRF == "6.Disabled"|
                                                                       ex.sample.temp.F.2$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                     as.character(predict_function(model1=models_males$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2))))
        ex.sample.temp.M.2$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2$LBRF == "3.Unemployed"|
                                                                       ex.sample.temp.M.2$LBRF == "5.Retired"|
                                                                       ex.sample.temp.M.2$LBRF == "6.Disabled"|
                                                                       ex.sample.temp.M.2$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                     as.character(predict_function(model1=models_males$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2))))
      } else if (intervention=="LBRF+OCC+INCOME" | intervention=="ALL + PAST INCOME"){
        ## equalize employment, occupation and income opportunities. can also intervene on past income ##
        
        ##intervene on income ##
        ex.sample.temp.F.2$ind_earn_adj <- predict_function(model1=models_males$ind_earn_adj, dataset=ex.sample.temp.F.2,
                                                            quantiles=quant)
        ex.sample.temp.M.2$ind_earn_adj <- predict_function(model1=models_males$ind_earn_adj, dataset=ex.sample.temp.M.2,
                                                            quantiles=quant)
        
        ## intervene here so that l.LBRF is changed ##
        ## take the observed distribution of employment in males and randomly assign it to females ##
        ex.sample.temp.F.2$LBRF <- predict_function(models_males$LBRF.model.M, ex.sample.temp.F.2)
        
        ## predict occupation ##
        ## intervene occupation ##
        ex.sample.temp.F.2$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2$LBRF == "3.Unemployed"|
                                                                       ex.sample.temp.F.2$LBRF == "5.Retired"|
                                                                       ex.sample.temp.F.2$LBRF == "6.Disabled"|
                                                                       ex.sample.temp.F.2$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                     as.character(predict_function(model1=models_males$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2))))
        ex.sample.temp.M.2$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2$LBRF == "3.Unemployed"|
                                                                       ex.sample.temp.M.2$LBRF == "5.Retired"|
                                                                       ex.sample.temp.M.2$LBRF == "6.Disabled"|
                                                                       ex.sample.temp.M.2$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                     as.character(predict_function(model1=models_males$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2))))
      } 
      
      else {return("Intervention is not correctly specified!")}
      
      ## join this newly created data with the old data #
      ex.sample.F.2 <- rbind(ex.sample.F.2,ex.sample.temp.F.2)
      rm(ex.sample.temp.F.2)
      
      ## order by ID variable #
      ex.sample.F.2 <- ex.sample.F.2[order(ex.sample.F.2$idnr,ex.sample.F.2$t),]
      
      ## males ##
      ex.sample.M.2 <- rbind(ex.sample.M.2,ex.sample.temp.M.2)
      rm(ex.sample.temp.M.2)
      
      ## order by ID variable ##
      ex.sample.M.2 <- ex.sample.M.2[order(ex.sample.M.2$idnr,ex.sample.M.2$t),]
    }
    
    ## save MC results ##
    ## before we save our results, we will subset != retired ##
    ## we do this because we are not interested in the fact of retirement on this but the effect of the other lbrf variables! ##
    
    ex.sample.F.noRET <- ex.sample.F[ex.sample.F$LBRF!="5.Retired" & ex.sample.F$LBRF!="6.Disabled",]
    ex.sample.M.noRET <- ex.sample.M[ex.sample.M$LBRF!="5.Retired" & ex.sample.M$LBRF!="6.Disabled",]
    
    ex.sample.F.2.noRET <- ex.sample.F.2[ex.sample.F.2$LBRF!="5.Retired" & ex.sample.F.2$LBRF!="6.Disabled",]
    ex.sample.M.2.noRET <- ex.sample.M.2[ex.sample.M.2$LBRF!="5.Retired" & ex.sample.M.2$LBRF!="6.Disabled",]
    
    ## aggregate by t ##
    for (i in 1:length(time.var_names)){
      if(time.var_names[i] == "ind_earn_adj"){
        tab.F <- ex.sample.F.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F[[i]][,m+1] <- tab.F$mean
        
        tab.M <- ex.sample.M.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M[[i]][,m+1] <- tab.M$mean
        
        tab.F2 <- ex.sample.F.2.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.2[[i]][,m+1] <- tab.F2$mean
        
        tab.M2 <- ex.sample.M.2.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.2[[i]][,m+1] <- tab.M2$mean
        
      } else if(length(levels(dataframe[,time.var_names[i]]))==2){
        mc.results.arr.F[[i]][,m+1] <- prop.table(table(ex.sample.F.noRET$t, ex.sample.F.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M[[i]][,m+1] <- prop.table(table(ex.sample.M.noRET$t, ex.sample.M.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.F.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.noRET$t, ex.sample.F.2.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.noRET$t, ex.sample.M.2.noRET[,time.var_names[i]]),margin=1)[,2]

      } else {
        mc.results.arr.F[[i]][,m+1] <- t(prop.table(table(ex.sample.F.noRET$t, ex.sample.F.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M[[i]][,m+1] <- t(prop.table(table(ex.sample.M.noRET$t, ex.sample.M.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.F.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.noRET$t, ex.sample.F.2.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.noRET$t, ex.sample.M.2.noRET[,time.var_names[i]]),margin=1))

      }
    }
  }
  
  ## save results in a list ## 
  results <- list(mc.results.arr.F,mc.results.arr.M,
                  mc.results.arr.F.2,mc.results.arr.M.2)
  names(results) <- c("mc.results.arr.F","mc.results.arr.M", ## these are the natural course approximations
                      "mc.results.arr.F.2","mc.results.arr.M.2") ## this is the counterfactual/intervention scenario
  
  return(results)
}



# LONGITUDINAL G FORMULA FUNCTION FOR SUBGROUP ANALYSIS BY RACE/ETHNICITY ------------------------

## for more detailed comments within the function, see previous section!

long.gformula_ethn <- function(dataframe,dataframe.t1=df_gform.t1_imp,id,formula.outcome, time.var_names, mcsize=50,
                               tunits, quant=c(0.25,0.5,0.75,0.90,0.95),intervention) { 
  
  
  ## sample
  # sample individuals (with replacement) from ex.dat
  re.size <- length(unique(dataframe[,id]))
  ex.sample.F.w <- cluster.resample(dataframe[dataframe$RAGENDER=="female" & dataframe$ethn == "1. non-Hisp White",], cluster.name = id, 
                                    size = re.size)
  ex.sample.F.b <- cluster.resample(dataframe[dataframe$RAGENDER=="female"& dataframe$ethn == "2. non-Hisp Black",], cluster.name = id, 
                                    size = re.size)
  ex.sample.F.h <- cluster.resample(dataframe[dataframe$RAGENDER=="female"& dataframe$ethn == "3.Hispanic",], cluster.name = id, 
                                    size = re.size)
  
  ex.sample.M.w <- cluster.resample(dataframe[dataframe$RAGENDER=="male" & dataframe$ethn == "1. non-Hisp White",], cluster.name = id, 
                                    size = re.size)
  ex.sample.M.b <- cluster.resample(dataframe[dataframe$RAGENDER=="male"& dataframe$ethn == "2. non-Hisp Black",], cluster.name = id, 
                                    size = re.size)
  ex.sample.M.h <- cluster.resample(dataframe[dataframe$RAGENDER=="male"& dataframe$ethn == "3.Hispanic",], cluster.name = id, 
                                    size = re.size)
  
  ## IMPUTATION CODE - GLM
  # RAFEDUC
  fit_RAFED_F.w <- multinom(RAFEDUC_cat ~ education + RAMEDUC_cat + psych_probl_E +
                              marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat + occupation_class_3cat +
                              LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                              l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                              l.occupation_class_3cat + 
                              l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F.w)
  
  fit_RAFED_F.b    <- update(fit_RAFED_F.w,data=ex.sample.F.b)
  fit_RAFED_F.h    <- update(fit_RAFED_F.w,data=ex.sample.F.h)
  
  fit_RAFED_M.w    <- update(fit_RAFED_F.w,data=ex.sample.M.w)
  fit_RAFED_M.b    <- update(fit_RAFED_F.w,data=ex.sample.M.b)
  fit_RAFED_M.h    <- update(fit_RAFED_F.w,data=ex.sample.M.h)
  
  # RAMEDUC
  fit_RAMED_F.w <- multinom(RAMEDUC_cat ~ education + RAFEDUC_cat  + psych_probl_E +
                              marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat + occupation_class_3cat +
                              LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                              l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                              l.occupation_class_3cat + 
                              l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F.w)
  
  fit_RAMED_F.b    <- update(fit_RAMED_F.w,data=ex.sample.F.b)
  fit_RAMED_F.h    <- update(fit_RAMED_F.w,data=ex.sample.F.h)
  
  fit_RAMED_M.w    <- update(fit_RAMED_F.w,data=ex.sample.M.w)
  fit_RAMED_M.b    <- update(fit_RAMED_F.w,data=ex.sample.M.b)
  fit_RAMED_M.h    <- update(fit_RAMED_F.w,data=ex.sample.M.h)
  
  #occ
  fit_occ_F.w <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E +
                            marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat +
                            LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                            l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                            l.occupation_class_3cat + 
                            l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F.w)
  fit_occ_F.b    <- update(fit_occ_F.w,data=ex.sample.F.b)
  fit_occ_F.h    <- update(fit_occ_F.w,data=ex.sample.F.h)
  
  fit_occ_M.w    <- update(fit_occ_F.w,data=ex.sample.M.w)
  fit_occ_M.b    <- update(fit_occ_F.w,data=ex.sample.M.b)
  fit_occ_M.h    <- update(fit_occ_F.w,data=ex.sample.M.h)
  
  #only for t1 and more
  #for t1 we dont need to predict anymore, we can just use the already imputed dataframe.t1
  #RAFED
  ex.sample.F.w$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.w$RAFEDUC_cat) & ex.sample.F.w$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "1. non-Hisp White"]),
                                                ifelse(is.na(ex.sample.F.w$RAFEDUC_cat) & ex.sample.F.w$t!=1, as.character(predict(fit_RAFED_F.w,type = "class")),
                                                       as.character(ex.sample.F.w$RAFEDUC_cat))))
  
  ex.sample.F.b$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.b$RAFEDUC_cat) & ex.sample.F.b$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "2. non-Hisp Black"]),
                                                ifelse(is.na(ex.sample.F.b$RAFEDUC_cat) & ex.sample.F.b$t!=1, as.character(predict(fit_RAFED_F.b,type = "class")),
                                                       as.character(ex.sample.F.b$RAFEDUC_cat))))
  
  ex.sample.F.h$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.h$RAFEDUC_cat) & ex.sample.F.h$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "3.Hispanic"]),
                                                ifelse(is.na(ex.sample.F.h$RAFEDUC_cat) & ex.sample.F.h$t!=1, as.character(predict(fit_RAFED_F.h,type = "class")),
                                                       as.character(ex.sample.F.h$RAFEDUC_cat))))
  
  
  ex.sample.M.w$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.w$RAFEDUC_cat) & ex.sample.M.w$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "1. non-Hisp White"]),
                                                ifelse(is.na(ex.sample.M.w$RAFEDUC_cat) & ex.sample.M.w$t!=1, as.character(predict(fit_RAFED_M.w,type = "class")),
                                                       as.character(ex.sample.M.w$RAFEDUC_cat))))
  
  ex.sample.M.b$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.b$RAFEDUC_cat) & ex.sample.M.b$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "2. non-Hisp Black"]),
                                                ifelse(is.na(ex.sample.M.b$RAFEDUC_cat) & ex.sample.M.b$t!=1, as.character(predict(fit_RAFED_M.b,type = "class")),
                                                       as.character(ex.sample.M.b$RAFEDUC_cat))))
  
  ex.sample.M.h$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.h$RAFEDUC_cat) & ex.sample.M.h$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "3.Hispanic"]),
                                                ifelse(is.na(ex.sample.M.h$RAFEDUC_cat) & ex.sample.M.h$t!=1, as.character(predict(fit_RAFED_M.h,type = "class")),
                                                       as.character(ex.sample.M.h$RAFEDUC_cat))))
  
  #RAMED
  ex.sample.F.w$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.w$RAMEDUC_cat) & ex.sample.F.w$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "1. non-Hisp White"]),
                                                ifelse(is.na(ex.sample.F.w$RAMEDUC_cat) & ex.sample.F.w$t!=1, as.character(predict(fit_RAMED_F.w,type = "class")),
                                                       as.character(ex.sample.F.w$RAMEDUC_cat))))
  
  ex.sample.F.b$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.b$RAMEDUC_cat) & ex.sample.F.b$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "2. non-Hisp Black"]),
                                                ifelse(is.na(ex.sample.F.b$RAMEDUC_cat) & ex.sample.F.b$t!=1, as.character(predict(fit_RAMED_F.b,type = "class")),
                                                       as.character(ex.sample.F.b$RAMEDUC_cat))))
  
  ex.sample.F.h$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.h$RAMEDUC_cat) & ex.sample.F.h$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "3.Hispanic"]),
                                                ifelse(is.na(ex.sample.F.h$RAMEDUC_cat) & ex.sample.F.h$t!=1, as.character(predict(fit_RAMED_F.h,type = "class")),
                                                       as.character(ex.sample.F.h$RAMEDUC_cat))))
  
  
  ex.sample.M.w$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.w$RAMEDUC_cat) & ex.sample.M.w$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "1. non-Hisp White"]),
                                                ifelse(is.na(ex.sample.M.w$RAMEDUC_cat) & ex.sample.M.w$t!=1, as.character(predict(fit_RAMED_M.w,type = "class")),
                                                       as.character(ex.sample.M.w$RAMEDUC_cat))))
  
  ex.sample.M.b$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.b$RAMEDUC_cat) & ex.sample.M.b$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "2. non-Hisp Black"]),
                                                ifelse(is.na(ex.sample.M.b$RAMEDUC_cat) & ex.sample.M.b$t!=1, as.character(predict(fit_RAMED_M.b,type = "class")),
                                                       as.character(ex.sample.M.b$RAMEDUC_cat))))
  
  ex.sample.M.h$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.h$RAMEDUC_cat) & ex.sample.M.h$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "3.Hispanic"]),
                                                ifelse(is.na(ex.sample.M.h$RAMEDUC_cat) & ex.sample.M.h$t!=1, as.character(predict(fit_RAMED_M.h,type = "class")),
                                                       as.character(ex.sample.M.h$RAMEDUC_cat))))
  
  #occ
  ex.sample.F.w$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F.w$occupation_class_4cat) & ex.sample.F.w$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "1. non-Hisp White"]),
                                                          ifelse(is.na(ex.sample.F.w$occupation_class_4cat) & ex.sample.F.w$t!=1, as.character(predict(fit_occ_F.w,type = "class")),
                                                                 as.character(ex.sample.F.w$occupation_class_4cat))))
  
  ex.sample.F.b$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F.b$occupation_class_4cat) & ex.sample.F.b$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "2. non-Hisp Black"]),
                                                          ifelse(is.na(ex.sample.F.b$occupation_class_4cat) & ex.sample.F.b$t!=1, as.character(predict(fit_occ_F.b,type = "class")),
                                                                 as.character(ex.sample.F.b$occupation_class_4cat))))
  
  ex.sample.F.h$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F.h$occupation_class_4cat) & ex.sample.F.h$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$ethn == "3.Hispanic"]),
                                                          ifelse(is.na(ex.sample.F.h$occupation_class_4cat) & ex.sample.F.h$t!=1, as.character(predict(fit_occ_F.h,type = "class")),
                                                                 as.character(ex.sample.F.h$occupation_class_4cat))))
  
  
  
  ex.sample.M.w$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M.w$occupation_class_4cat) & ex.sample.M.w$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "1. non-Hisp White"]),
                                                          ifelse(is.na(ex.sample.M.w$occupation_class_4cat) & ex.sample.M.w$t!=1, as.character(predict(fit_occ_M.w,type = "class")),
                                                                 as.character(ex.sample.M.w$occupation_class_4cat))))
  
  ex.sample.M.b$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M.b$occupation_class_4cat) & ex.sample.M.b$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "2. non-Hisp Black"]),
                                                          ifelse(is.na(ex.sample.M.b$occupation_class_4cat) & ex.sample.M.b$t!=1, as.character(predict(fit_occ_M.b,type = "class")),
                                                                 as.character(ex.sample.M.b$occupation_class_4cat))))
  
  ex.sample.M.h$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M.h$occupation_class_4cat) & ex.sample.M.h$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$ethn == "3.Hispanic"]),
                                                          ifelse(is.na(ex.sample.M.h$occupation_class_4cat) & ex.sample.M.h$t!=1, as.character(predict(fit_occ_M.h,type = "class")),
                                                                 as.character(ex.sample.M.h$occupation_class_4cat))))
  
  
  
  # remake 3 class occ and l.occupation
  ex.sample.F.w$occupation_class_3cat <- as.factor(ifelse(ex.sample.F.w$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F.w$occupation_class_4cat)))
  ex.sample.F.b$occupation_class_3cat <- as.factor(ifelse(ex.sample.F.b$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F.b$occupation_class_4cat)))
  ex.sample.F.h$occupation_class_3cat <- as.factor(ifelse(ex.sample.F.h$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F.h$occupation_class_4cat)))
  
  ex.sample.M.w$occupation_class_3cat <- as.factor(ifelse(ex.sample.M.w$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M.w$occupation_class_4cat)))
  ex.sample.M.b$occupation_class_3cat <- as.factor(ifelse(ex.sample.M.b$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M.b$occupation_class_4cat)))
  ex.sample.M.h$occupation_class_3cat <- as.factor(ifelse(ex.sample.M.h$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M.h$occupation_class_4cat)))
  
  ex.sample.F.w$l.occupation_class_4cat <- lag1.func(ex.sample.F.w$occupation_class_4cat, ex.sample.F.w$HHIDPN, ex.sample.F.w$wave)
  ex.sample.F.b$l.occupation_class_4cat <- lag1.func(ex.sample.F.b$occupation_class_4cat, ex.sample.F.b$HHIDPN, ex.sample.F.b$wave)
  ex.sample.F.h$l.occupation_class_4cat <- lag1.func(ex.sample.F.h$occupation_class_4cat, ex.sample.F.h$HHIDPN, ex.sample.F.h$wave)
  
  ex.sample.M.w$l.occupation_class_4cat <- lag1.func(ex.sample.M.w$occupation_class_4cat, ex.sample.M.w$HHIDPN, ex.sample.M.w$wave)
  ex.sample.M.b$l.occupation_class_4cat <- lag1.func(ex.sample.M.b$occupation_class_4cat, ex.sample.M.b$HHIDPN, ex.sample.M.b$wave)
  ex.sample.M.h$l.occupation_class_4cat <- lag1.func(ex.sample.M.h$occupation_class_4cat, ex.sample.M.h$HHIDPN, ex.sample.M.h$wave)
  
  ## IMPUTATION CODE - GLM
  
  # (re)fit models to ex.sample
  models_females.w <- models_females.b <- models_females.h <- vector(mode = "list", length = length(time.var_names))
  models_males.w <- models_males.b <- models_males.h <- vector(mode = "list", length = length(time.var_names))
  
  for (i in 1:length(time.var_names)){
    
    formula <- as.formula(sub('outcome',paste0(time.var_names[i]),formula.outcome))
    
    if (class(dataframe[,time.var_names[i]]) == "numeric") {
      
      #multiple quantiles
      
      # we will use quantile regression here because the distribution
      # of income looks pretty bad (zero inflated and outliers)
      # For larger problems it is advantageous to use the Frisch-Newton interior point method "fn"
      models_females.w[[i]] <- rq(formula, data=ex.sample.F.w, tau=quant, method="fn") #replace quant with 0.5 for median regression only
      models_females.b[[i]] <- rq(formula, data=ex.sample.F.b, tau=quant, method="fn")
      models_females.h[[i]] <- rq(formula, data=ex.sample.F.h, tau=quant, method="fn")
      
      models_males.w[[i]] <- rq(formula, data=ex.sample.M.w, tau=quant, method="fn") #replace quant with 0.5 for median regression only
      models_males.b[[i]] <- rq(formula, data=ex.sample.M.b, tau=quant, method="fn")
      models_males.h[[i]] <- rq(formula, data=ex.sample.M.h, tau=quant, method="fn")
      
    } else if (length(unique(dataframe[,time.var_names[i]])) == 2){
      
      models_females.w[[i]] <- glm(formula, family=binomial, data=ex.sample.F.w)
      models_females.b[[i]] <- glm(formula, family=binomial, data=ex.sample.F.b)
      models_females.h[[i]] <- glm(formula, family=binomial, data=ex.sample.F.h)
      
      models_males.w[[i]] <- glm(formula, family=binomial, data=ex.sample.M.w)
      models_males.b[[i]] <- glm(formula, family=binomial, data=ex.sample.M.b)
      models_males.h[[i]] <- glm(formula, family=binomial, data=ex.sample.M.h)
      
    } else if (class(dataframe[,time.var_names[i]]) == "factor" & length(unique(dataframe[,time.var_names[i]])) > 2){         
      
      models_females.w[[i]] <- multinom(formula, data=ex.sample.F.w)
      models_females.b[[i]] <- multinom(formula, data=ex.sample.F.b)
      models_females.h[[i]] <- multinom(formula, data=ex.sample.F.h)
      
      models_males.w[[i]] <- multinom(formula, data=ex.sample.M.w)
      models_males.b[[i]] <- multinom(formula, data=ex.sample.M.b)
      models_males.h[[i]] <- multinom(formula, data=ex.sample.M.h)
      
    } else if (time.var_names == "occupation_class"){  
      # occupation model has slightly difference model specification
      models_females.w[[i]] <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                          l.H_child_cat + l.depress +
                                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                        data=ex.sample.F.w)
      models_females.b[[i]] <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                          l.H_child_cat + l.depress +
                                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                        data=ex.sample.F.b)
      models_females.h[[i]] <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                          l.H_child_cat + l.depress +
                                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                        data=ex.sample.F.h)
      
      models_males.w[[i]] <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.M.w)
      models_males.b[[i]] <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.M.b)
      models_males.h[[i]] <- multinom(occupation_class_3cat ~ education + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.M.h)
      
      
    } else {print("Time-varying variable has no unique levels/values!")}
  }
  
  names(models_females.w) <- names(models_females.b) <- names(models_females.h) <- lapply(time.var_names, FUN = paste0,".model.F")
  names(models_males.w)   <- names(models_males.b)   <- names(models_males.h)    <- lapply(time.var_names, FUN = paste0,".model.M")
  
  #make mc results saving arrays
  mc.results.arr.F.w  <- list()
  
  for (i in 1:length(time.var_names)){
    if (length(levels(dataframe[,time.var_names[i]])) > 2){
      matrix  <- matrix(NA,(tunits*length(levels(dataframe[,time.var_names[i]]))),1+mcsize)   
      matrix[,1] <- rep(seq(1,tunits,1),each=length(levels(dataframe[,time.var_names[i]])))
      rownames(matrix) <- rep(levels(dataframe[,time.var_names[i]]),tunits)
      mc.results.arr.F.w[[i]] <- matrix
    } else {
      matrix <- matrix(NA,tunits,1+mcsize)
      matrix[,1] <- seq(1,tunits,1)
      mc.results.arr.F.w[[i]]  <- matrix
    }
  }
  names(mc.results.arr.F.w) <- time.var_names
  
  mc.results.arr.F.h.2 <- mc.results.arr.F.h <- mc.results.arr.F.b.2 <- mc.results.arr.F.b <- mc.results.arr.F.w.2 <- mc.results.arr.F.w
  mc.results.arr.M.h.2 <- mc.results.arr.M.h <- mc.results.arr.M.b.2 <- mc.results.arr.M.b <- mc.results.arr.M.w.2 <- mc.results.arr.M.w <- mc.results.arr.F.w
  

  for(m in 1:mcsize) {
    
    # take individuals at time 1
    # and discard the other observations
    
    ######## start simulation with imputation at age 50 
    ex.sample.F.w <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$ethn=="1. non-Hisp White",]
    ex.sample.F.b <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$ethn=="2. non-Hisp Black",]
    ex.sample.F.h <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$ethn=="3.Hispanic",]
    
    ex.sample.M.w <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male" & df_gform.t1_imp$ethn=="1. non-Hisp White",]
    ex.sample.M.b <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male" & df_gform.t1_imp$ethn=="2. non-Hisp Black",]
    ex.sample.M.h <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male" & df_gform.t1_imp$ethn=="3.Hispanic",]
    
    # give idnr
    ex.sample.F.w$idnr <- 1:length(ex.sample.F.w[,id])
    ex.sample.F.b$idnr <- 1:length(ex.sample.F.b[,id])
    ex.sample.F.h$idnr <- 1:length(ex.sample.F.h[,id])
    
    ex.sample.M.w$idnr <- 1:length(ex.sample.M.w[,id])
    ex.sample.M.b$idnr <- 1:length(ex.sample.M.b[,id])
    ex.sample.M.h$idnr <- 1:length(ex.sample.M.h[,id])
    
    # make copies for intervention and mediation simulations
    ex.sample.F.2.w <- ex.sample.F.w
    ex.sample.F.2.b <- ex.sample.F.b
    ex.sample.F.2.h <- ex.sample.F.h
    
    ex.sample.M.2.w <- ex.sample.M.w
    ex.sample.M.2.b <- ex.sample.M.b
    ex.sample.M.2.h <- ex.sample.M.h
    
    if (intervention=="ALL + PAST INCOME"){
      #### preparation for l.income intervention
      ## take income levels at t1 (age 50, 51) in males 
      ## and calculate mean levels based on lbrf status and occupational group
      ## if we use more grouping variables we need to consider using a regression model
      ## because a prop table will have scarcity issues
      ## and then give females new income levels at t1 based on that prop.table
      
      df_calc_mean.inc <- rbind(ex.sample.F.2.w, ex.sample.F.2.b,ex.sample.F.2.h, 
                                ex.sample.M.2.w,ex.sample.M.2.b,ex.sample.M.2.h)
      df_calc_mean.inc <- df_calc_mean.inc %>% group_by(ethn,LBRF,occupation_class_4cat) %>%
        mutate(ind_earn_adj = mean(ind_earn_adj[RAGENDER=="male"])) %>%
        ungroup()

      ex.sample.F.2.w <-  df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female" & df_calc_mean.inc$ethn=="1. non-Hisp White",]   
      ex.sample.F.2.w <- as.data.frame(ex.sample.F.2.w)
      
      ex.sample.F.2.b <-  df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female" & df_calc_mean.inc$ethn=="2. non-Hisp Black",]   
      ex.sample.F.2.b <- as.data.frame(ex.sample.F.2.b)
      
      ex.sample.F.2.h <-  df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female" & df_calc_mean.inc$ethn=="3.Hispanic",]   
      ex.sample.F.2.h <- as.data.frame(ex.sample.F.2.h)
      
      
    }
    
    # NATURAL COURSE APPROXIMATION LOOP ---------------------------------------
    # start a loop that moves through the follow-up time units
    # this part of the g-formula tries to reproduce the empirical data
    # and is known as the 'natural course'
    for(t in 2:tunits) {


      
      #female
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.F.w <- ex.sample.F.w[ex.sample.F.w$t==(t-1),]
      ex.sample.temp.F.w$t <- ex.sample.temp.F.w$t+1
      ex.sample.temp.F.w$AGEY_E <- ex.sample.temp.F.w$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.w[,col.index.to] <- ex.sample.temp.F.w[,col.index.from]
      
      #hisp
      ex.sample.temp.F.b <- ex.sample.F.b[ex.sample.F.b$t==(t-1),]
      ex.sample.temp.F.b$t <- ex.sample.temp.F.b$t+1
      ex.sample.temp.F.b$AGEY_E <- ex.sample.temp.F.b$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.b[,col.index.to] <- ex.sample.temp.F.b[,col.index.from]
      
      #black
      ex.sample.temp.F.h <- ex.sample.F.h[ex.sample.F.h$t==(t-1),]
      ex.sample.temp.F.h$t <- ex.sample.temp.F.h$t+1
      ex.sample.temp.F.h$AGEY_E <- ex.sample.temp.F.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.h[,col.index.to] <- ex.sample.temp.F.h[,col.index.from]
      
      #males
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.M.w <- ex.sample.M.w[ex.sample.M.w$t==(t-1),]
      ex.sample.temp.M.w$t <- ex.sample.temp.M.w$t+1
      ex.sample.temp.M.w$AGEY_E <- ex.sample.temp.M.w$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.w[,col.index.to] <- ex.sample.temp.M.w[,col.index.from]
      
      #hisp
      ex.sample.temp.M.b <- ex.sample.M.b[ex.sample.M.b$t==(t-1),]
      ex.sample.temp.M.b$t <- ex.sample.temp.M.b$t+1
      ex.sample.temp.M.b$AGEY_E <- ex.sample.temp.M.b$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.b[,col.index.to] <- ex.sample.temp.M.b[,col.index.from]
      
      #black
      ex.sample.temp.M.h <- ex.sample.M.h[ex.sample.M.h$t==(t-1),]
      ex.sample.temp.M.h$t <- ex.sample.temp.M.h$t+1
      ex.sample.temp.M.h$AGEY_E <- ex.sample.temp.M.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.h[,col.index.to] <- ex.sample.temp.M.h[,col.index.from]
      
      
      # predict time varying variables
      # occupation is predicted slightly differently 
      l <- length(time.var_names)-2
      for (i in 1:l){
        ex.sample.temp.F.w[,time.var_names[i]] <- predict_function(model1=models_females.w[[i]], dataset=ex.sample.temp.F.w)
        ex.sample.temp.F.b[,time.var_names[i]] <- predict_function(model1=models_females.b[[i]], dataset=ex.sample.temp.F.b)
        ex.sample.temp.F.h[,time.var_names[i]] <- predict_function(model1=models_females.h[[i]], dataset=ex.sample.temp.F.h)
        
        ex.sample.temp.M.w[,time.var_names[i]] <- predict_function(model1=models_males.w[[i]], dataset=ex.sample.temp.M.w)
        ex.sample.temp.M.b[,time.var_names[i]] <- predict_function(model1=models_males.b[[i]], dataset=ex.sample.temp.M.b)
        ex.sample.temp.M.h[,time.var_names[i]] <- predict_function(model1=models_males.h[[i]], dataset=ex.sample.temp.M.h)
      }
      
      ex.sample.temp.F.w$ind_earn_adj <- predict_function(model1=models_females.w$ind_earn_adj, dataset=ex.sample.temp.F.w,
                                                          quantiles=quant)
      ex.sample.temp.F.b$ind_earn_adj <- predict_function(model1=models_females.b$ind_earn_adj, dataset=ex.sample.temp.F.b,
                                                          quantiles=quant)
      ex.sample.temp.F.h$ind_earn_adj <- predict_function(model1=models_females.h$ind_earn_adj, dataset=ex.sample.temp.F.h,
                                                          quantiles=quant)
      
      
      ex.sample.temp.M.w$ind_earn_adj <- predict_function(model1=models_males.w$ind_earn_adj, dataset=ex.sample.temp.M.w,
                                                          quantiles=quant)
      ex.sample.temp.M.b$ind_earn_adj <- predict_function(model1=models_males.b$ind_earn_adj, dataset=ex.sample.temp.M.b,
                                                          quantiles=quant)
      ex.sample.temp.M.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.h,
                                                          quantiles=quant)
      
      
      ## predict occupation
      ex.sample.temp.F.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.w$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.F.w$LBRF == "5.Retired"|
                                                                     ex.sample.temp.F.w$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.F.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_females.w$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.w))))
      ex.sample.temp.F.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.b$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.F.b$LBRF == "5.Retired"|
                                                                     ex.sample.temp.F.b$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.F.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_females.b$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.b))))
      
      ex.sample.temp.F.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.h$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.F.h$LBRF == "5.Retired"|
                                                                     ex.sample.temp.F.h$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.F.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_females.h$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.h))))
      
      
      ex.sample.temp.M.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.w$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.M.w$LBRF == "5.Retired"|
                                                                     ex.sample.temp.M.w$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.M.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_males.w$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.w))))
      ex.sample.temp.M.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.b$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.M.b$LBRF == "5.Retired"|
                                                                     ex.sample.temp.M.b$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.M.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_males.b$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.b))))
      
      ex.sample.temp.M.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.h$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.M.h$LBRF == "5.Retired"|
                                                                     ex.sample.temp.M.h$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.M.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.h))))
      
      
      #white
      # join this newly created data with the old data
      ex.sample.F.w <- rbind(ex.sample.F.w,ex.sample.temp.F.w)
      rm(ex.sample.temp.F.w)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.F.w <- ex.sample.F.w[order(ex.sample.F.w$idnr,ex.sample.F.w$t),]
      
      #hisp
      ex.sample.F.b <- rbind(ex.sample.F.b,ex.sample.temp.F.b)
      rm(ex.sample.temp.F.b)
      
      ex.sample.F.b <- ex.sample.F.b[order(ex.sample.F.b$idnr,ex.sample.F.b$t),]
      
      #black
      ex.sample.F.h <- rbind(ex.sample.F.h,ex.sample.temp.F.h)
      rm(ex.sample.temp.F.h)
      
      ex.sample.F.h <- ex.sample.F.h[order(ex.sample.F.h$idnr,ex.sample.F.h$t),]
      
      
      #males
      #white
      # join this newly created data with the old data
      ex.sample.M.w <- rbind(ex.sample.M.w,ex.sample.temp.M.w)
      rm(ex.sample.temp.M.w)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.M.w <- ex.sample.M.w[order(ex.sample.M.w$idnr,ex.sample.M.w$t),]
      
      #hisp
      ex.sample.M.b <- rbind(ex.sample.M.b,ex.sample.temp.M.b)
      rm(ex.sample.temp.M.b)
      
      ex.sample.M.b <- ex.sample.M.b[order(ex.sample.M.b$idnr,ex.sample.M.b$t),]
      
      #black
      ex.sample.M.h <- rbind(ex.sample.M.h,ex.sample.temp.M.h)
      rm(ex.sample.temp.M.h)
      
      ex.sample.M.h <- ex.sample.M.h[order(ex.sample.M.h$idnr,ex.sample.M.h$t),]

      
    }
    

# INTERVENTION LOOP -------------------------------------------------------

    ## now let's create a counterfactual dataset for the TOTAL EFFECT
    # we will equalise LBRF distribution among males and females
    #assuming that females have the same LBRF distribution as males
    
    ## for each time unit ## 
    for(t in 2:tunits) {
      
      #female
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.F.2.w <- ex.sample.F.2.w[ex.sample.F.2.w$t==(t-1),]
      ex.sample.temp.F.2.w$t <- ex.sample.temp.F.2.w$t+1
      ex.sample.temp.F.2.w$AGEY_E <- ex.sample.temp.F.2.w$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.2.w[,col.index.to] <- ex.sample.temp.F.2.w[,col.index.from]
      
      #hisp
      ex.sample.temp.F.2.b <- ex.sample.F.2.b[ex.sample.F.2.b$t==(t-1),]
      ex.sample.temp.F.2.b$t <- ex.sample.temp.F.2.b$t+1
      ex.sample.temp.F.2.b$AGEY_E <- ex.sample.temp.F.2.b$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.2.b[,col.index.to] <- ex.sample.temp.F.2.b[,col.index.from]
      
      #black
      ex.sample.temp.F.2.h <- ex.sample.F.2.h[ex.sample.F.2.h$t==(t-1),]
      ex.sample.temp.F.2.h$t <- ex.sample.temp.F.2.h$t+1
      ex.sample.temp.F.2.h$AGEY_E <- ex.sample.temp.F.2.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.2.h[,col.index.to] <- ex.sample.temp.F.2.h[,col.index.from]
      
      
      #male
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.M.2.w <- ex.sample.M.2.w[ex.sample.M.2.w$t==(t-1),]
      ex.sample.temp.M.2.w$t <- ex.sample.temp.M.2.w$t+1
      ex.sample.temp.M.2.w$AGEY_E <- ex.sample.temp.M.2.w$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.2.w[,col.index.to] <- ex.sample.temp.M.2.w[,col.index.from]
      
      #hisp
      ex.sample.temp.M.2.b <- ex.sample.M.2.b[ex.sample.M.2.b$t==(t-1),]
      ex.sample.temp.M.2.b$t <- ex.sample.temp.M.2.b$t+1
      ex.sample.temp.M.2.b$AGEY_E <- ex.sample.temp.M.2.b$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.2.b[,col.index.to] <- ex.sample.temp.M.2.b[,col.index.from]
      
      #black
      ex.sample.temp.M.2.h <- ex.sample.M.2.h[ex.sample.M.2.h$t==(t-1),]
      ex.sample.temp.M.2.h$t <- ex.sample.temp.M.2.h$t+1
      ex.sample.temp.M.2.h$AGEY_E <- ex.sample.temp.M.2.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.2.h[,col.index.to] <- ex.sample.temp.M.2.h[,col.index.from]
      
      
      
      ## predict time varying variables ##
      ## occupation is predicted slightly differently ##
      l <- length(time.var_names)-2
      for (i in 1:l){
        ex.sample.temp.F.2.w[,time.var_names[i]] <- predict_function(model1=models_females.w[[i]], dataset=ex.sample.temp.F.2.w)
        ex.sample.temp.F.2.b[,time.var_names[i]] <- predict_function(model1=models_females.b[[i]], dataset=ex.sample.temp.F.2.b)
        ex.sample.temp.F.2.h[,time.var_names[i]] <- predict_function(model1=models_females.h[[i]], dataset=ex.sample.temp.F.2.h)
        
        ex.sample.temp.M.2.w[,time.var_names[i]] <- predict_function(model1=models_males.w[[i]], dataset=ex.sample.temp.M.2.w)
        ex.sample.temp.M.2.b[,time.var_names[i]] <- predict_function(model1=models_males.b[[i]], dataset=ex.sample.temp.M.2.b)
        ex.sample.temp.M.2.h[,time.var_names[i]] <- predict_function(model1=models_males.h[[i]], dataset=ex.sample.temp.M.2.h)
      }
      
      if (intervention=="LBRF"){
        ## predict income ##
        ex.sample.temp.F.2.w$ind_earn_adj <- predict_function(model1=models_females.w$ind_earn_adj, dataset=ex.sample.temp.F.2.w,
                                                              quantiles=quant)
        ex.sample.temp.F.2.b$ind_earn_adj <- predict_function(model1=models_females.b$ind_earn_adj, dataset=ex.sample.temp.F.2.b,
                                                              quantiles=quant)
        ex.sample.temp.F.2.h$ind_earn_adj <- predict_function(model1=models_females.h$ind_earn_adj, dataset=ex.sample.temp.F.2.h,
                                                              quantiles=quant)
        
        ex.sample.temp.M.2.w$ind_earn_adj <- predict_function(model1=models_males.w$ind_earn_adj, dataset=ex.sample.temp.M.2.w,
                                                              quantiles=quant)
        ex.sample.temp.M.2.b$ind_earn_adj <- predict_function(model1=models_males.b$ind_earn_adj, dataset=ex.sample.temp.M.2.b,
                                                              quantiles=quant)
        ex.sample.temp.M.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.2.h,
                                                              quantiles=quant)
        
        #intervene here so that l.LBRF is changed
        # take the observed distribution of employment in males and randomly assign it to females
        ex.sample.temp.F.2.w$LBRF <- predict_function(models_males.w$LBRF.model.M, ex.sample.temp.F.2.w)
        ex.sample.temp.F.2.b$LBRF <- predict_function(models_males.b$LBRF.model.M, ex.sample.temp.F.2.b)
        ex.sample.temp.F.2.h$LBRF <- predict_function(models_males.h$LBRF.model.M, ex.sample.temp.F.2.h)
        
        
        ## predict occupation
        ex.sample.temp.F.2.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.w$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.w$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.w$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_females.w$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2.w))))
        ex.sample.temp.F.2.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.b$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.b$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.b$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_females.b$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2.b))))
        
        ex.sample.temp.F.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_females.h$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2.h))))
        
        ex.sample.temp.M.2.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.w$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.w$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.w$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.w$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.w))))
        ex.sample.temp.M.2.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.b$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.b$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.b$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.b$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.b))))
        
        ex.sample.temp.M.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.h))))
        
        
        
      } else if (intervention=="LBRF+OCC"){
        
        ## predict income
        ex.sample.temp.F.2.w$ind_earn_adj <- predict_function(model1=models_females.w$ind_earn_adj, dataset=ex.sample.temp.F.2.w,
                                                              quantiles=quant)
        ex.sample.temp.F.2.b$ind_earn_adj <- predict_function(model1=models_females.b$ind_earn_adj, dataset=ex.sample.temp.F.2.b,
                                                              quantiles=quant)
        ex.sample.temp.F.2.h$ind_earn_adj <- predict_function(model1=models_females.h$ind_earn_adj, dataset=ex.sample.temp.F.2.h,
                                                              quantiles=quant)
        
        
        ex.sample.temp.M.2.w$ind_earn_adj <- predict_function(model1=models_males.w$ind_earn_adj, dataset=ex.sample.temp.M.2.w,
                                                              quantiles=quant)
        ex.sample.temp.M.2.b$ind_earn_adj <- predict_function(model1=models_males.b$ind_earn_adj, dataset=ex.sample.temp.M.2.b,
                                                              quantiles=quant)
        ex.sample.temp.M.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.2.h,
                                                              quantiles=quant)
        
        #intervene here so that l.LBRF is changed
        # take the observed distribution of employment in males and randomly assign it to females
        ex.sample.temp.F.2.w$LBRF <- predict_function(models_males.w$LBRF.model.M, ex.sample.temp.F.2.w)
        ex.sample.temp.F.2.b$LBRF <- predict_function(models_males.b$LBRF.model.M, ex.sample.temp.F.2.b)
        ex.sample.temp.F.2.h$LBRF <- predict_function(models_males.h$LBRF.model.M, ex.sample.temp.F.2.h)
        
        
        ## predict occupation ##
        ## intervene on occupation ##
        ex.sample.temp.F.2.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.w$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.w$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.w$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.w$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.w))))
        ex.sample.temp.F.2.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.b$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.b$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.b$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.b$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.b))))
        
        ex.sample.temp.F.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.h))))
        
        
        ex.sample.temp.M.2.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.w$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.w$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.w$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.w$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.w))))
        ex.sample.temp.M.2.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.b$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.b$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.b$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.b$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.b))))
        
        ex.sample.temp.M.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.h))))
        
        
      } else if (intervention=="LBRF+OCC+INCOME" | intervention=="ALL + PAST INCOME"){
        
        ## intervene on income ##
        ex.sample.temp.F.2.w$ind_earn_adj <- predict_function(model1=models_males.w$ind_earn_adj, dataset=ex.sample.temp.F.2.w,
                                                              quantiles=quant)
        ex.sample.temp.F.2.b$ind_earn_adj <- predict_function(model1=models_males.b$ind_earn_adj, dataset=ex.sample.temp.F.2.b,
                                                              quantiles=quant)
        ex.sample.temp.F.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.F.2.h,
                                                              quantiles=quant)
        
        
        ex.sample.temp.M.2.w$ind_earn_adj <- predict_function(model1=models_males.w$ind_earn_adj, dataset=ex.sample.temp.M.2.w,
                                                              quantiles=quant)
        ex.sample.temp.M.2.b$ind_earn_adj <- predict_function(model1=models_males.b$ind_earn_adj, dataset=ex.sample.temp.M.2.b,
                                                              quantiles=quant)
        ex.sample.temp.M.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.2.h,
                                                              quantiles=quant)
        
        # intervene here so that l.LBRF is changed
        # take the observed distribution of employment in males and randomly assign it to females
        ex.sample.temp.F.2.w$LBRF <- predict_function(models_males.w$LBRF.model.M, ex.sample.temp.F.2.w)
        ex.sample.temp.F.2.b$LBRF <- predict_function(models_males.b$LBRF.model.M, ex.sample.temp.F.2.b)
        ex.sample.temp.F.2.h$LBRF <- predict_function(models_males.h$LBRF.model.M, ex.sample.temp.F.2.h)
        
        
        ## predict occupation ##
        ## intervene on occupation ##
        ex.sample.temp.F.2.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.w$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.w$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.w$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.w$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.w))))
        ex.sample.temp.F.2.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.b$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.b$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.b$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.b$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.b))))
        
        ex.sample.temp.F.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.h))))
        
        ex.sample.temp.M.2.w$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.w$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.w$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.w$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.w$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.w$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.w))))
        ex.sample.temp.M.2.b$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.b$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.b$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.b$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.b$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.b$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.b))))
        
        ex.sample.temp.M.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.h))))
        
      } else {return("Intervention is not correctly specified!")}
      
      
      # white
      # join this newly created data with the old data
      ex.sample.F.2.w <- rbind(ex.sample.F.2.w,ex.sample.temp.F.2.w)
      rm(ex.sample.temp.F.2.w)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.F.2.w <- ex.sample.F.2.w[order(ex.sample.F.2.w$idnr,ex.sample.F.2.w$t),]
      
      # hisp
      ex.sample.F.2.b <- rbind(ex.sample.F.2.b,ex.sample.temp.F.2.b)
      rm(ex.sample.temp.F.2.b)
      
      ex.sample.F.2.b <- ex.sample.F.2.b[order(ex.sample.F.2.b$idnr,ex.sample.F.2.b$t),]
      
      # black
      ex.sample.F.2.h <- rbind(ex.sample.F.2.h,ex.sample.temp.F.2.h)
      rm(ex.sample.temp.F.2.h)
      
      ex.sample.F.2.h <- ex.sample.F.2.h[order(ex.sample.F.2.h$idnr,ex.sample.F.2.h$t),]
      
      
      # males
      # white
      # join this newly created data with the old data
      ex.sample.M.2.w <- rbind(ex.sample.M.2.w,ex.sample.temp.M.2.w)
      rm(ex.sample.temp.M.2.w)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.M.2.w <- ex.sample.M.2.w[order(ex.sample.M.2.w$idnr,ex.sample.M.2.w$t),]
      
      # hisp
      ex.sample.M.2.b <- rbind(ex.sample.M.2.b,ex.sample.temp.M.2.b)
      rm(ex.sample.temp.M.2.b)
      
      ex.sample.M.2.b <- ex.sample.M.2.b[order(ex.sample.M.2.b$idnr,ex.sample.M.2.b$t),]
      
      # black
      ex.sample.M.2.h <- rbind(ex.sample.M.2.h,ex.sample.temp.M.2.h)
      rm(ex.sample.temp.M.2.h)
      
      ex.sample.M.2.h <- ex.sample.M.2.h[order(ex.sample.M.2.h$idnr,ex.sample.M.2.h$t),]
      
      
      
    }
    
    # save mc results
    # before we save our results, we will subset != retired
    # we do this because we are not interested in the fact of retirement on this but the effect of the other lbrf variables!
    
    ex.sample.F.w.noRET <- ex.sample.F.w[ex.sample.F.w$LBRF!="5.Retired" & ex.sample.F.w$LBRF!="6.Disabled",]
    ex.sample.F.b.noRET <- ex.sample.F.b[ex.sample.F.b$LBRF!="5.Retired" & ex.sample.F.b$LBRF!="6.Disabled",]
    ex.sample.F.h.noRET <- ex.sample.F.h[ex.sample.F.h$LBRF!="5.Retired" & ex.sample.F.h$LBRF!="6.Disabled",]
    
    ex.sample.M.w.noRET <- ex.sample.M.w[ex.sample.M.w$LBRF!="5.Retired" & ex.sample.M.w$LBRF!="6.Disabled",]
    ex.sample.M.b.noRET <- ex.sample.M.b[ex.sample.M.b$LBRF!="5.Retired" & ex.sample.M.b$LBRF!="6.Disabled",]
    ex.sample.M.h.noRET <- ex.sample.M.h[ex.sample.M.h$LBRF!="5.Retired" & ex.sample.M.h$LBRF!="6.Disabled",]
    
    ex.sample.F.2.w.noRET <- ex.sample.F.2.w[ex.sample.F.2.w$LBRF!="5.Retired" & ex.sample.F.2.w$LBRF!="6.Disabled",]
    ex.sample.F.2.b.noRET <- ex.sample.F.2.b[ex.sample.F.2.b$LBRF!="5.Retired" & ex.sample.F.2.b$LBRF!="6.Disabled",]
    ex.sample.F.2.h.noRET <- ex.sample.F.2.h[ex.sample.F.2.h$LBRF!="5.Retired" & ex.sample.F.2.h$LBRF!="6.Disabled",]
    
    ex.sample.M.2.w.noRET <- ex.sample.M.2.w[ex.sample.M.2.w$LBRF!="5.Retired" & ex.sample.M.2.w$LBRF!="6.Disabled",]
    ex.sample.M.2.b.noRET <- ex.sample.M.2.b[ex.sample.M.2.b$LBRF!="5.Retired" & ex.sample.M.2.b$LBRF!="6.Disabled",]
    ex.sample.M.2.h.noRET <- ex.sample.M.2.h[ex.sample.M.2.h$LBRF!="5.Retired" & ex.sample.M.2.h$LBRF!="6.Disabled",]
    
    
    for (i in 1:length(time.var_names)){
      if(time.var_names[i] == "ind_earn_adj"){
        tab.F.w <- ex.sample.F.w.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.w[[i]][,m+1] <- tab.F.w$mean
        tab.F.b <- ex.sample.F.b.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.b[[i]][,m+1] <- tab.F.b$mean
        tab.F.h <- ex.sample.F.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.h[[i]][,m+1] <- tab.F.h$mean
        
        tab.M.w <- ex.sample.M.w.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.w[[i]][,m+1] <- tab.M.w$mean
        tab.M.b <- ex.sample.M.b.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.b[[i]][,m+1] <- tab.M.b$mean
        tab.M.h <- ex.sample.M.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.h[[i]][,m+1] <- tab.M.h$mean
        
        
        
        tab.F2.w <- ex.sample.F.2.w.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.w.2[[i]][,m+1] <- tab.F2.w$mean
        tab.F2.b <- ex.sample.F.2.b.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.b.2[[i]][,m+1] <- tab.F2.b$mean
        tab.F2.h <- ex.sample.F.2.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.h.2[[i]][,m+1] <- tab.F2.h$mean
        
        
        tab.M2.w <- ex.sample.M.2.w.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.w.2[[i]][,m+1] <- tab.M2.w$mean
        tab.M2.b <- ex.sample.M.2.b.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.b.2[[i]][,m+1] <- tab.M2.b$mean
        tab.M2.h <- ex.sample.M.2.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.h.2[[i]][,m+1] <- tab.M2.h$mean
        
        
      } else if(length(levels(dataframe[,time.var_names[i]]))==2){
        mc.results.arr.F.w[[i]][,m+1] <- prop.table(table(ex.sample.F.w.noRET$t, ex.sample.F.w.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.b[[i]][,m+1] <- prop.table(table(ex.sample.F.b.noRET$t, ex.sample.F.b.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.h[[i]][,m+1] <- prop.table(table(ex.sample.F.h.noRET$t, ex.sample.F.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.M.w[[i]][,m+1] <- prop.table(table(ex.sample.M.w.noRET$t, ex.sample.M.w.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.b[[i]][,m+1] <- prop.table(table(ex.sample.M.b.noRET$t, ex.sample.M.b.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.h[[i]][,m+1] <- prop.table(table(ex.sample.M.h.noRET$t, ex.sample.M.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.F.w.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.w.noRET$t, ex.sample.F.2.w.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.b.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.b.noRET$t, ex.sample.F.2.b.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.h.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.h.noRET$t, ex.sample.F.2.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.M.w.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.w.noRET$t, ex.sample.M.2.w.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.b.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.b.noRET$t, ex.sample.M.2.b.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.h.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.h.noRET$t, ex.sample.M.2.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
      } else {
        mc.results.arr.F.w[[i]][,m+1] <- t(prop.table(table(ex.sample.F.w.noRET$t, ex.sample.F.w.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.b[[i]][,m+1] <- t(prop.table(table(ex.sample.F.b.noRET$t, ex.sample.F.b.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.h[[i]][,m+1] <- t(prop.table(table(ex.sample.F.h.noRET$t, ex.sample.F.h.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.M.w[[i]][,m+1] <- t(prop.table(table(ex.sample.M.w.noRET$t, ex.sample.M.w.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.b[[i]][,m+1] <- t(prop.table(table(ex.sample.M.b.noRET$t, ex.sample.M.b.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.h[[i]][,m+1] <- t(prop.table(table(ex.sample.M.h.noRET$t, ex.sample.M.h.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.F.w.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.w.noRET$t, ex.sample.F.2.w.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.b.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.b.noRET$t, ex.sample.F.2.b.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.h.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.h.noRET$t, ex.sample.F.2.h.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.M.w.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.w.noRET$t, ex.sample.M.2.w.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.b.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.b.noRET$t, ex.sample.M.2.b.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.h.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.h.noRET$t, ex.sample.M.2.h.noRET[,time.var_names[i]]),margin=1))
        
      }
    }
  }
  #need to return all of these in one list ideally!
  results <- list(mc.results.arr.F.w,mc.results.arr.F.b,mc.results.arr.F.h,
                  mc.results.arr.M.w,mc.results.arr.M.b,mc.results.arr.M.h,
                  mc.results.arr.F.w.2,mc.results.arr.F.b.2,mc.results.arr.F.h.2,
                  mc.results.arr.M.w.2,mc.results.arr.M.b.2,mc.results.arr.M.h.2)
  names(results) <- c("mc.results.arr.F.w","mc.results.arr.F.b","mc.results.arr.F.h",
                      "mc.results.arr.M.w","mc.results.arr.M.b","mc.results.arr.M.h",
                      "mc.results.arr.F.w.2","mc.results.arr.F.b.2","mc.results.arr.F.h.2",
                      "mc.results.arr.M.w.2","mc.results.arr.M.b.2","mc.results.arr.M.h.2")
  
  return(results)
}

# LONGITUDINAL G FORMULA FUNCTION FOR SUBGROUP ANALYSIS BY EDUCATION ------------------------

## for more detailed comments within the function, see previous section!

long.gformula_educ <- function(dataframe,dataframe.t1=df_gform.t1_imp,id,formula.outcome, time.var_names, mcsize=50,
                               tunits, quant=c(0.25,0.5,0.75,0.90,0.95),intervention) { 
  
  
  ## sample
  # sample individuals (with replacement) from ex.dat
  re.size <- length(unique(dataframe[,id]))
  ex.sample.F.l <- cluster.resample(dataframe[dataframe$RAGENDER=="female" & dataframe$education == "Lt High-school/GED",], cluster.name = id, 
                                    size = re.size)
  ex.sample.F.m <- cluster.resample(dataframe[dataframe$RAGENDER=="female"& dataframe$education == "High-school graduate",], cluster.name = id, 
                                    size = re.size)
  ex.sample.F.h <- cluster.resample(dataframe[dataframe$RAGENDER=="female"& dataframe$education == "Some college or higher",], cluster.name = id, 
                                    size = re.size)
  
  ex.sample.M.l <- cluster.resample(dataframe[dataframe$RAGENDER=="male" & dataframe$education == "Lt High-school/GED",], cluster.name = id, 
                                    size = re.size)
  ex.sample.M.m <- cluster.resample(dataframe[dataframe$RAGENDER=="male"& dataframe$education == "High-school graduate",], cluster.name = id, 
                                    size = re.size)
  ex.sample.M.h <- cluster.resample(dataframe[dataframe$RAGENDER=="male"& dataframe$education == "Some college or higher",], cluster.name = id, 
                                    size = re.size)
  
  
  ## IMPUTATION CODE - GLM
  # RAFEDUC
  fit_RAFED_F.l <- multinom(RAFEDUC_cat ~ ethn + RAMEDUC_cat + psych_probl_E +
                              marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat + occupation_class_3cat +
                              LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                              l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                              l.occupation_class_3cat + 
                              l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F.l)
  
  fit_RAFED_F.m    <- update(fit_RAFED_F.l,data=ex.sample.F.m)
  fit_RAFED_F.h    <- update(fit_RAFED_F.l,data=ex.sample.F.h)
  
  fit_RAFED_M.l    <- update(fit_RAFED_F.l,data=ex.sample.M.l)
  fit_RAFED_M.m    <- update(fit_RAFED_F.l,data=ex.sample.M.m)
  fit_RAFED_M.h    <- update(fit_RAFED_F.l,data=ex.sample.M.h)
  
  # RAMEDUC
  fit_RAMED_F.l <- multinom(RAMEDUC_cat ~ ethn + RAFEDUC_cat  + psych_probl_E +
                              marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat + occupation_class_3cat +
                              LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                              l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                              l.occupation_class_3cat + 
                              l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F.l)
  
  fit_RAMED_F.m    <- update(fit_RAMED_F.l,data=ex.sample.F.m)
  fit_RAMED_F.h    <- update(fit_RAMED_F.l,data=ex.sample.F.h)
  
  fit_RAMED_M.l    <- update(fit_RAMED_F.l,data=ex.sample.M.l)
  fit_RAMED_M.m    <- update(fit_RAMED_F.l,data=ex.sample.M.m)
  fit_RAMED_M.h    <- update(fit_RAMED_F.l,data=ex.sample.M.h)
  
  #occ
  fit_occ_F.l <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E +
                            marital_stat_cat2 + H_members + H_child_cat + depress + N_cond_cat +
                            LBRF * (ns(AGEY_E,3) + ind_earn_adj) +
                            l.marital_stat_cat2 + l.H_members + l.H_child_cat + l.depress + l.N_cond_cat +
                            l.occupation_class_3cat + 
                            l.LBRF * (ns(AGEY_E,3) + l.ind_earn_adj), data=ex.sample.F.l)
  fit_occ_F.m    <- update(fit_occ_F.l,data=ex.sample.F.m)
  fit_occ_F.h    <- update(fit_occ_F.l,data=ex.sample.F.h)
  
  fit_occ_M.l    <- update(fit_occ_F.l,data=ex.sample.M.l)
  fit_occ_M.m    <- update(fit_occ_F.l,data=ex.sample.M.m)
  fit_occ_M.h    <- update(fit_occ_F.l,data=ex.sample.M.h)
  
  #only for t1 and more
  #for t1 we dont need to predict anymore, we can just use the already imputed dataframe.t1
  #RAFED
  ex.sample.F.l$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.l$RAFEDUC_cat) & ex.sample.F.l$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "Lt High-school/GED"]),
                                                ifelse(is.na(ex.sample.F.l$RAFEDUC_cat) & ex.sample.F.l$t!=1, as.character(predict(fit_RAFED_F.l,type = "class")),
                                                       as.character(ex.sample.F.l$RAFEDUC_cat))))
  
  ex.sample.F.m$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.m$RAFEDUC_cat) & ex.sample.F.m$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "High-school graduate"]),
                                                ifelse(is.na(ex.sample.F.m$RAFEDUC_cat) & ex.sample.F.m$t!=1, as.character(predict(fit_RAFED_F.m,type = "class")),
                                                       as.character(ex.sample.F.m$RAFEDUC_cat))))
  
  ex.sample.F.h$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.h$RAFEDUC_cat) & ex.sample.F.h$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "Some college or higher"]),
                                                ifelse(is.na(ex.sample.F.h$RAFEDUC_cat) & ex.sample.F.h$t!=1, as.character(predict(fit_RAFED_F.h,type = "class")),
                                                       as.character(ex.sample.F.h$RAFEDUC_cat))))
  
  
  ex.sample.M.l$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.l$RAFEDUC_cat) & ex.sample.M.l$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "Lt High-school/GED"]),
                                                ifelse(is.na(ex.sample.M.l$RAFEDUC_cat) & ex.sample.M.l$t!=1, as.character(predict(fit_RAFED_M.l,type = "class")),
                                                       as.character(ex.sample.M.l$RAFEDUC_cat))))
  
  ex.sample.M.m$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.m$RAFEDUC_cat) & ex.sample.M.m$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "High-school graduate"]),
                                                ifelse(is.na(ex.sample.M.m$RAFEDUC_cat) & ex.sample.M.m$t!=1, as.character(predict(fit_RAFED_M.m,type = "class")),
                                                       as.character(ex.sample.M.m$RAFEDUC_cat))))
  
  ex.sample.M.h$RAFEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.h$RAFEDUC_cat) & ex.sample.M.h$t==1,
                                                as.character(dataframe.t1$RAFEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "Some college or higher"]),
                                                ifelse(is.na(ex.sample.M.h$RAFEDUC_cat) & ex.sample.M.h$t!=1, as.character(predict(fit_RAFED_M.h,type = "class")),
                                                       as.character(ex.sample.M.h$RAFEDUC_cat))))
  
  #RAMED
  ex.sample.F.l$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.l$RAMEDUC_cat) & ex.sample.F.l$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "Lt High-school/GED"]),
                                                ifelse(is.na(ex.sample.F.l$RAMEDUC_cat) & ex.sample.F.l$t!=1, as.character(predict(fit_RAMED_F.l,type = "class")),
                                                       as.character(ex.sample.F.l$RAMEDUC_cat))))
  
  ex.sample.F.m$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.m$RAMEDUC_cat) & ex.sample.F.m$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "High-school graduate"]),
                                                ifelse(is.na(ex.sample.F.m$RAMEDUC_cat) & ex.sample.F.m$t!=1, as.character(predict(fit_RAMED_F.m,type = "class")),
                                                       as.character(ex.sample.F.m$RAMEDUC_cat))))
  
  ex.sample.F.h$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.F.h$RAMEDUC_cat) & ex.sample.F.h$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "Some college or higher"]),
                                                ifelse(is.na(ex.sample.F.h$RAMEDUC_cat) & ex.sample.F.h$t!=1, as.character(predict(fit_RAMED_F.h,type = "class")),
                                                       as.character(ex.sample.F.h$RAMEDUC_cat))))
  
  
  ex.sample.M.l$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.l$RAMEDUC_cat) & ex.sample.M.l$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "Lt High-school/GED"]),
                                                ifelse(is.na(ex.sample.M.l$RAMEDUC_cat) & ex.sample.M.l$t!=1, as.character(predict(fit_RAMED_M.l,type = "class")),
                                                       as.character(ex.sample.M.l$RAMEDUC_cat))))
  
  ex.sample.M.m$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.m$RAMEDUC_cat) & ex.sample.M.m$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "High-school graduate"]),
                                                ifelse(is.na(ex.sample.M.m$RAMEDUC_cat) & ex.sample.M.m$t!=1, as.character(predict(fit_RAMED_M.m,type = "class")),
                                                       as.character(ex.sample.M.m$RAMEDUC_cat))))
  
  ex.sample.M.h$RAMEDUC_cat <- as.factor(ifelse(is.na(ex.sample.M.h$RAMEDUC_cat) & ex.sample.M.h$t==1,
                                                as.character(dataframe.t1$RAMEDUC_cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "Some college or higher"]),
                                                ifelse(is.na(ex.sample.M.h$RAMEDUC_cat) & ex.sample.M.h$t!=1, as.character(predict(fit_RAMED_M.h,type = "class")),
                                                       as.character(ex.sample.M.h$RAMEDUC_cat))))
  #occ
  ex.sample.F.l$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F.l$occupation_class_4cat) & ex.sample.F.l$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "Lt High-school/GED"]),
                                                          ifelse(is.na(ex.sample.F.l$occupation_class_4cat) & ex.sample.F.l$t!=1, as.character(predict(fit_occ_F.l,type = "class")),
                                                                 as.character(ex.sample.F.l$occupation_class_4cat))))
  
  ex.sample.F.m$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F.m$occupation_class_4cat) & ex.sample.F.m$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "High-school graduate"]),
                                                          ifelse(is.na(ex.sample.F.m$occupation_class_4cat) & ex.sample.F.m$t!=1, as.character(predict(fit_occ_F.m,type = "class")),
                                                                 as.character(ex.sample.F.m$occupation_class_4cat))))
  
  ex.sample.F.h$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.F.h$occupation_class_4cat) & ex.sample.F.h$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="female" & dataframe.t1$education == "Some college or higher"]),
                                                          ifelse(is.na(ex.sample.F.h$occupation_class_4cat) & ex.sample.F.h$t!=1, as.character(predict(fit_occ_F.h,type = "class")),
                                                                 as.character(ex.sample.F.h$occupation_class_4cat))))
  
  
  ex.sample.M.l$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M.l$occupation_class_4cat) & ex.sample.M.l$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "Lt High-school/GED"]),
                                                          ifelse(is.na(ex.sample.M.l$occupation_class_4cat) & ex.sample.M.l$t!=1, as.character(predict(fit_occ_M.l,type = "class")),
                                                                 as.character(ex.sample.M.l$occupation_class_4cat))))
  
  ex.sample.M.m$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M.m$occupation_class_4cat) & ex.sample.M.m$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "High-school graduate"]),
                                                          ifelse(is.na(ex.sample.M.m$occupation_class_4cat) & ex.sample.M.m$t!=1, as.character(predict(fit_occ_M.m,type = "class")),
                                                                 as.character(ex.sample.M.m$occupation_class_4cat))))
  
  ex.sample.M.h$occupation_class_4cat <- as.factor(ifelse(is.na(ex.sample.M.h$occupation_class_4cat) & ex.sample.M.h$t==1,
                                                          as.character(dataframe.t1$occupation_class_4cat[dataframe.t1$RAGENDER=="male" & dataframe.t1$education == "Some college or higher"]),
                                                          ifelse(is.na(ex.sample.M.h$occupation_class_4cat) & ex.sample.M.h$t!=1, as.character(predict(fit_occ_M.h,type = "class")),
                                                                 as.character(ex.sample.M.h$occupation_class_4cat))))
  
  
  
  # remake 3 class occ and l.occupation
  ex.sample.F.l$occupation_class_3cat <- as.factor(ifelse(ex.sample.F.l$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F.l$occupation_class_4cat)))
  ex.sample.F.m$occupation_class_3cat <- as.factor(ifelse(ex.sample.F.m$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F.m$occupation_class_4cat)))
  ex.sample.F.h$occupation_class_3cat <- as.factor(ifelse(ex.sample.F.h$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.F.h$occupation_class_4cat)))
  
  ex.sample.M.l$occupation_class_3cat <- as.factor(ifelse(ex.sample.M.l$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M.l$occupation_class_4cat)))
  ex.sample.M.m$occupation_class_3cat <- as.factor(ifelse(ex.sample.M.m$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M.m$occupation_class_4cat)))
  ex.sample.M.h$occupation_class_3cat <- as.factor(ifelse(ex.sample.M.h$occupation_class_4cat == "4.doesnt work",NA,as.character(ex.sample.M.h$occupation_class_4cat)))
  
  ex.sample.F.l$l.occupation_class_4cat <- lag1.func(ex.sample.F.l$occupation_class_4cat, ex.sample.F.l$HHIDPN, ex.sample.F.l$wave)
  ex.sample.F.m$l.occupation_class_4cat <- lag1.func(ex.sample.F.m$occupation_class_4cat, ex.sample.F.m$HHIDPN, ex.sample.F.m$wave)
  ex.sample.F.h$l.occupation_class_4cat <- lag1.func(ex.sample.F.h$occupation_class_4cat, ex.sample.F.h$HHIDPN, ex.sample.F.h$wave)
  
  ex.sample.M.l$l.occupation_class_4cat <- lag1.func(ex.sample.M.l$occupation_class_4cat, ex.sample.M.l$HHIDPN, ex.sample.M.l$wave)
  ex.sample.M.m$l.occupation_class_4cat <- lag1.func(ex.sample.M.m$occupation_class_4cat, ex.sample.M.m$HHIDPN, ex.sample.M.m$wave)
  ex.sample.M.h$l.occupation_class_4cat <- lag1.func(ex.sample.M.h$occupation_class_4cat, ex.sample.M.h$HHIDPN, ex.sample.M.h$wave)
  
  ## IMPUTATION CODE - GLM
  
  # (re)fit models to ex.sample
  models_females.l <- models_females.m <- models_females.h <- vector(mode = "list", length = length(time.var_names))
  models_males.l <- models_males.m <- models_males.h <- vector(mode = "list", length = length(time.var_names))
  
  for (i in 1:length(time.var_names)){
    
    formula <- as.formula(sub('outcome',paste0(time.var_names[i]),formula.outcome))
    
    if (class(dataframe[,time.var_names[i]]) == "numeric") {
      
      #multiple quantiles
      
      # we will use quantile regression here because the distribution
      # of income looks pretty bad (zero inflated and outliers)
      # For larger problems it is advantageous to use the Frisch-Newton interior point method "fn"
      models_females.l[[i]] <- rq(formula, data=ex.sample.F.l, tau=quant, method="fn") #replace quant with 0.5 for median regression only
      models_females.m[[i]] <- rq(formula, data=ex.sample.F.m, tau=quant, method="fn")
      models_females.h[[i]] <- rq(formula, data=ex.sample.F.h, tau=quant, method="fn")
      
      models_males.l[[i]] <- rq(formula, data=ex.sample.M.l, tau=quant, method="fn") #replace quant with 0.5 for median regression only
      models_males.m[[i]] <- rq(formula, data=ex.sample.M.m, tau=quant, method="fn")
      models_males.h[[i]] <- rq(formula, data=ex.sample.M.h, tau=quant, method="fn")
      
    } else if (length(unique(dataframe[,time.var_names[i]])) == 2){
      
      models_females.l[[i]] <- glm(formula, family=binomial, data=ex.sample.F.l)
      models_females.m[[i]] <- glm(formula, family=binomial, data=ex.sample.F.m)
      models_females.h[[i]] <- glm(formula, family=binomial, data=ex.sample.F.h)
      
      models_males.l[[i]] <- glm(formula, family=binomial, data=ex.sample.M.l)
      models_males.m[[i]] <- glm(formula, family=binomial, data=ex.sample.M.m)
      models_males.h[[i]] <- glm(formula, family=binomial, data=ex.sample.M.h)
      
    } else if (class(dataframe[,time.var_names[i]]) == "factor" & length(unique(dataframe[,time.var_names[i]])) > 2){         
      
      models_females.l[[i]] <- multinom(formula, data=ex.sample.F.l)
      models_females.m[[i]] <- multinom(formula, data=ex.sample.F.m)
      models_females.h[[i]] <- multinom(formula, data=ex.sample.F.h)
      
      models_males.l[[i]] <- multinom(formula, data=ex.sample.M.l)
      models_males.m[[i]] <- multinom(formula, data=ex.sample.M.m)
      models_males.h[[i]] <- multinom(formula, data=ex.sample.M.h)
      
    } else if (time.var_names == "occupation_class"){  
      # occupation model has slightly difference model specification
      models_females.l[[i]] <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                          l.H_child_cat + l.depress +
                                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                        data=ex.sample.F.l)
      models_females.m[[i]] <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                          l.H_child_cat + l.depress +
                                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                        data=ex.sample.F.m)
      models_females.h[[i]] <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                          l.H_child_cat + l.depress +
                                          LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                        data=ex.sample.F.h)
      
      models_males.l[[i]] <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.M.l)
      models_males.m[[i]] <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.M.m)
      models_males.h[[i]] <- multinom(occupation_class_3cat ~ ethn + RAFEDUC_cat + RAMEDUC_cat + psych_probl_E + l.marital_stat_cat2 + l.H_members +
                                        l.H_child_cat + l.depress +
                                        LBRF * (ns(AGEY_E,3) + ind_earn_adj) + l.N_cond_cat,
                                      data=ex.sample.M.h)
      
      
    } else {print("Time-varying variable has no unique levels/values!")}
  }
  
  names(models_females.l) <- names(models_females.m) <- names(models_females.h) <- lapply(time.var_names, FUN = paste0,".model.F")
  names(models_males.l)   <- names(models_males.m)   <- names(models_males.h)   <- lapply(time.var_names, FUN = paste0,".model.M")
  
  #make mc results saving arrays
  mc.results.arr.F.l  <- list()
  
  for (i in 1:length(time.var_names)){
    if (length(levels(dataframe[,time.var_names[i]])) > 2){
      matrix  <- matrix(NA,(tunits*length(levels(dataframe[,time.var_names[i]]))),1+mcsize)   
      matrix[,1] <- rep(seq(1,tunits,1),each=length(levels(dataframe[,time.var_names[i]])))
      rownames(matrix) <- rep(levels(dataframe[,time.var_names[i]]),tunits)
      mc.results.arr.F.l[[i]] <- matrix
    } else {
      matrix <- matrix(NA,tunits,1+mcsize)
      matrix[,1] <- seq(1,tunits,1)
      mc.results.arr.F.l[[i]]  <- matrix
    }
  }
  names(mc.results.arr.F.l) <- time.var_names
  
  mc.results.arr.F.h.2 <- mc.results.arr.F.h <- mc.results.arr.F.m.2 <- mc.results.arr.F.m <- mc.results.arr.F.l.2 <- mc.results.arr.F.l  
  mc.results.arr.M.h.2 <- mc.results.arr.M.h <- mc.results.arr.M.m.2 <- mc.results.arr.M.m <- mc.results.arr.M.l.2 <- mc.results.arr.M.l <- mc.results.arr.F.l
  
  for(m in 1:mcsize) {
    
    # take individuals at time 1
    # and discard the other observations
    
    ######## start simulation with imputation at age 50 
    ex.sample.F.l <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$education=="Lt High-school/GED",]
    ex.sample.F.m <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$education=="High-school graduate",]
    ex.sample.F.h <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="female" & df_gform.t1_imp$education=="Some college or higher",]
    
    ex.sample.M.l <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male" & df_gform.t1_imp$education=="Lt High-school/GED",]
    ex.sample.M.m <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male" & df_gform.t1_imp$education=="High-school graduate",]
    ex.sample.M.h <- df_gform.t1_imp[df_gform.t1_imp$RAGENDER=="male" & df_gform.t1_imp$education=="Some college or higher",]
    
    # give idnr
    ex.sample.F.l$idnr <- 1:length(ex.sample.F.l[,id])
    ex.sample.F.m$idnr <- 1:length(ex.sample.F.m[,id])
    ex.sample.F.h$idnr <- 1:length(ex.sample.F.h[,id])
    
    ex.sample.M.l$idnr <- 1:length(ex.sample.M.l[,id])
    ex.sample.M.m$idnr <- 1:length(ex.sample.M.m[,id])
    ex.sample.M.h$idnr <- 1:length(ex.sample.M.h[,id])
    
    # make copies for intervention and mediation simulations
    ex.sample.F.2.l <- ex.sample.F.l
    ex.sample.F.2.m <- ex.sample.F.m
    ex.sample.F.2.h <- ex.sample.F.h
    
    ex.sample.M.2.l <- ex.sample.M.l
    ex.sample.M.2.m <- ex.sample.M.m
    ex.sample.M.2.h <- ex.sample.M.h
    
    if (intervention=="ALL + PAST INCOME"){
      #### preparation for l.income intervention
      ## take income levels at t1 (age 50, 51) in males 
      ## and calculate mean levels based on lbrf status and occupational group
      ## if we use more grouping variables we need to consider using a regression model
      ## because a prop table will have scarcity issues
      ## and then give females new income levels at t1 based on that prop.table
      
      df_calc_mean.inc <- rbind(ex.sample.F.2.l, ex.sample.F.2.m,ex.sample.F.2.h,
                                ex.sample.M.2.l,ex.sample.M.2.m,ex.sample.M.2.h)
      df_calc_mean.inc <- df_calc_mean.inc %>% group_by(education,LBRF,occupation_class_4cat) %>%
        mutate(ind_earn_adj = mean(ind_earn_adj[RAGENDER=="male"])) %>%
        ungroup()
      
      ex.sample.F.2.l <- df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female" & df_calc_mean.inc$education=="Lt High-school/GED",]
      ex.sample.F.2.l <- as.data.frame(ex.sample.F.2.l)
      
      ex.sample.F.2.m<- df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female" & df_calc_mean.inc$education=="High-school graduate",]
      ex.sample.F.2.m <- as.data.frame(ex.sample.F.2.m)
      
      ex.sample.F.2.h <- df_calc_mean.inc[df_calc_mean.inc$RAGENDER=="female" & df_calc_mean.inc$education=="Some college or higher",]
      ex.sample.F.2.h <- as.data.frame(ex.sample.F.2.h)
      
      
    }
    
    # NATURAL COURSE APPROXIMATION --------------------------------------------
    # start a loop that moves through the follow-up time units
    # this part of the g-formula tries to reproduce the empirical data
    # and is known as the 'natural course'
    
    ## for each time unit ##
    for(t in 2:tunits) {

      #female
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.F.l <- ex.sample.F.l[ex.sample.F.l$t==(t-1),]
      ex.sample.temp.F.l$t <- ex.sample.temp.F.l$t+1
      ex.sample.temp.F.l$AGEY_E <- ex.sample.temp.F.l$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.l[,col.index.to] <- ex.sample.temp.F.l[,col.index.from]
      
      #hisp
      ex.sample.temp.F.m <- ex.sample.F.m[ex.sample.F.m$t==(t-1),]
      ex.sample.temp.F.m$t <- ex.sample.temp.F.m$t+1
      ex.sample.temp.F.m$AGEY_E <- ex.sample.temp.F.m$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.m[,col.index.to] <- ex.sample.temp.F.m[,col.index.from]
      
      #black
      ex.sample.temp.F.h <- ex.sample.F.h[ex.sample.F.h$t==(t-1),]
      ex.sample.temp.F.h$t <- ex.sample.temp.F.h$t+1
      ex.sample.temp.F.h$AGEY_E <- ex.sample.temp.F.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.h[,col.index.to] <- ex.sample.temp.F.h[,col.index.from]
      
      
      #male
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.M.l <- ex.sample.M.l[ex.sample.M.l$t==(t-1),]
      ex.sample.temp.M.l$t <- ex.sample.temp.M.l$t+1
      ex.sample.temp.M.l$AGEY_E <- ex.sample.temp.M.l$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.l[,col.index.to] <- ex.sample.temp.M.l[,col.index.from]
      
      #hisp
      ex.sample.temp.M.m <- ex.sample.M.m[ex.sample.M.m$t==(t-1),]
      ex.sample.temp.M.m$t <- ex.sample.temp.M.m$t+1
      ex.sample.temp.M.m$AGEY_E <- ex.sample.temp.M.m$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.m[,col.index.to] <- ex.sample.temp.M.m[,col.index.from]
      
      #black
      ex.sample.temp.M.h <- ex.sample.M.h[ex.sample.M.h$t==(t-1),]
      ex.sample.temp.M.h$t <- ex.sample.temp.M.h$t+1
      ex.sample.temp.M.h$AGEY_E <- ex.sample.temp.M.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.h[,col.index.to] <- ex.sample.temp.M.h[,col.index.from]
      
      
      # predict time varying variables
      # occupation is predicted slightly differently 
      l <- length(time.var_names)-2
      for (i in 1:l){
        ex.sample.temp.F.l[,time.var_names[i]] <- predict_function(model1=models_females.l[[i]], dataset=ex.sample.temp.F.l)
        ex.sample.temp.F.m[,time.var_names[i]] <- predict_function(model1=models_females.m[[i]], dataset=ex.sample.temp.F.m)
        ex.sample.temp.F.h[,time.var_names[i]] <- predict_function(model1=models_females.h[[i]], dataset=ex.sample.temp.F.h)
        
        ex.sample.temp.M.l[,time.var_names[i]] <- predict_function(model1=models_males.l[[i]], dataset=ex.sample.temp.M.l)
        ex.sample.temp.M.m[,time.var_names[i]] <- predict_function(model1=models_males.m[[i]], dataset=ex.sample.temp.M.m)
        ex.sample.temp.M.h[,time.var_names[i]] <- predict_function(model1=models_males.h[[i]], dataset=ex.sample.temp.M.h)
      }
      
      ex.sample.temp.F.l$ind_earn_adj <- predict_function(model1=models_females.l$ind_earn_adj, dataset=ex.sample.temp.F.l,
                                                          quantiles=quant)
      ex.sample.temp.F.m$ind_earn_adj <- predict_function(model1=models_females.m$ind_earn_adj, dataset=ex.sample.temp.F.m,
                                                          quantiles=quant)
      ex.sample.temp.F.h$ind_earn_adj <- predict_function(model1=models_females.h$ind_earn_adj, dataset=ex.sample.temp.F.h,
                                                          quantiles=quant)
      
      ex.sample.temp.M.l$ind_earn_adj <- predict_function(model1=models_males.l$ind_earn_adj, dataset=ex.sample.temp.M.l,
                                                          quantiles=quant)
      ex.sample.temp.M.m$ind_earn_adj <- predict_function(model1=models_males.m$ind_earn_adj, dataset=ex.sample.temp.M.m,
                                                          quantiles=quant)
      ex.sample.temp.M.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.h,
                                                          quantiles=quant)
      
      
      ## predict occupation
      ex.sample.temp.F.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.l$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.F.l$LBRF == "5.Retired"|
                                                                     ex.sample.temp.F.l$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.F.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_females.l$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.l))))
      ex.sample.temp.F.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.m$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.F.m$LBRF == "5.Retired"|
                                                                     ex.sample.temp.F.m$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.F.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_females.m$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.m))))
      
      ex.sample.temp.F.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.h$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.F.h$LBRF == "5.Retired"|
                                                                     ex.sample.temp.F.h$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.F.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_females.h$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.h))))
      
      
      ex.sample.temp.M.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.l$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.M.l$LBRF == "5.Retired"|
                                                                     ex.sample.temp.M.l$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.M.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_males.l$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.l))))
      ex.sample.temp.M.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.m$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.M.m$LBRF == "5.Retired"|
                                                                     ex.sample.temp.M.m$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.M.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_males.m$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.m))))
      
      ex.sample.temp.M.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.h$LBRF == "3.Unemployed"|
                                                                     ex.sample.temp.M.h$LBRF == "5.Retired"|
                                                                     ex.sample.temp.M.h$LBRF == "6.Disabled"|
                                                                     ex.sample.temp.M.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                   as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.h))))
      
      
      #white
      # join this newly created data with the old data
      ex.sample.F.l <- rbind(ex.sample.F.l,ex.sample.temp.F.l)
      rm(ex.sample.temp.F.l)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.F.l <- ex.sample.F.l[order(ex.sample.F.l$idnr,ex.sample.F.l$t),]
      
      #hisp
      ex.sample.F.m <- rbind(ex.sample.F.m,ex.sample.temp.F.m)
      rm(ex.sample.temp.F.m)
      
      ex.sample.F.m <- ex.sample.F.m[order(ex.sample.F.m$idnr,ex.sample.F.m$t),]
      
      #black
      ex.sample.F.h <- rbind(ex.sample.F.h,ex.sample.temp.F.h)
      rm(ex.sample.temp.F.h)
      
      ex.sample.F.h <- ex.sample.F.h[order(ex.sample.F.h$idnr,ex.sample.F.h$t),]
      
      
      #males
      #white
      # join this newly created data with the old data
      ex.sample.M.l <- rbind(ex.sample.M.l,ex.sample.temp.M.l)
      rm(ex.sample.temp.M.l)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.M.l <- ex.sample.M.l[order(ex.sample.M.l$idnr,ex.sample.M.l$t),]
      
      #hisp
      ex.sample.M.m <- rbind(ex.sample.M.m,ex.sample.temp.M.m)
      rm(ex.sample.temp.M.m)
      
      ex.sample.M.m <- ex.sample.M.m[order(ex.sample.M.m$idnr,ex.sample.M.m$t),]
      
      #black
      ex.sample.M.h <- rbind(ex.sample.M.h,ex.sample.temp.M.h)
      rm(ex.sample.temp.M.h)
      
      ex.sample.M.h <- ex.sample.M.h[order(ex.sample.M.h$idnr,ex.sample.M.h$t),]
      
      
    }
    
# INTERVENTION LOOP -------------------------------------------------------

    ## now let's create a counterfactual dataset for the TOTAL EFFECT
    # we will equalise LBRF distribution among males and females
    #assuming that females have the same LBRF distribution as males
     
    ## for each time unit ##
    for(t in 2:tunits) {
      
      #female
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.F.2.l <- ex.sample.F.2.l[ex.sample.F.2.l$t==(t-1),]
      ex.sample.temp.F.2.l$t <- ex.sample.temp.F.2.l$t+1
      ex.sample.temp.F.2.l$AGEY_E <- ex.sample.temp.F.2.l$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.2.l[,col.index.to] <- ex.sample.temp.F.2.l[,col.index.from]
      
      #hisp
      ex.sample.temp.F.2.m <- ex.sample.F.2.m[ex.sample.F.2.m$t==(t-1),]
      ex.sample.temp.F.2.m$t <- ex.sample.temp.F.2.m$t+1
      ex.sample.temp.F.2.m$AGEY_E <- ex.sample.temp.F.2.m$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.2.m[,col.index.to] <- ex.sample.temp.F.2.m[,col.index.from]
      
      #black
      ex.sample.temp.F.2.h <- ex.sample.F.2.h[ex.sample.F.2.h$t==(t-1),]
      ex.sample.temp.F.2.h$t <- ex.sample.temp.F.2.h$t+1
      ex.sample.temp.F.2.h$AGEY_E <- ex.sample.temp.F.2.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.F.2.h[,col.index.to] <- ex.sample.temp.F.2.h[,col.index.from]
      
      #male
      #white
      # make a copy of the previous row and then update year and age
      # but don't make a copy if someone was censored in the previous year
      ex.sample.temp.M.2.l <- ex.sample.M.2.l[ex.sample.M.2.l$t==(t-1),]
      ex.sample.temp.M.2.l$t <- ex.sample.temp.M.2.l$t+1
      ex.sample.temp.M.2.l$AGEY_E <- ex.sample.temp.M.2.l$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.2.l[,col.index.to] <- ex.sample.temp.M.2.l[,col.index.from]
      
      #hisp
      ex.sample.temp.M.2.m <- ex.sample.M.2.m[ex.sample.M.2.m$t==(t-1),]
      ex.sample.temp.M.2.m$t <- ex.sample.temp.M.2.m$t+1
      ex.sample.temp.M.2.m$AGEY_E <- ex.sample.temp.M.2.m$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.2.m[,col.index.to] <- ex.sample.temp.M.2.m[,col.index.from]
      
      #black
      ex.sample.temp.M.2.h <- ex.sample.M.2.h[ex.sample.M.2.h$t==(t-1),]
      ex.sample.temp.M.2.h$t <- ex.sample.temp.M.2.h$t+1
      ex.sample.temp.M.2.h$AGEY_E <- ex.sample.temp.M.2.h$AGEY_E+2
      
      # lag values
      ex.sample.temp.M.2.h[,col.index.to] <- ex.sample.temp.M.2.h[,col.index.from]
      
      # predict time varying variables
      # occupation is predicted slightly differently 
      l <- length(time.var_names)-2
      for (i in 1:l){
        ex.sample.temp.F.2.l[,time.var_names[i]] <- predict_function(model1=models_females.l[[i]], dataset=ex.sample.temp.F.2.l)
        ex.sample.temp.F.2.m[,time.var_names[i]] <- predict_function(model1=models_females.m[[i]], dataset=ex.sample.temp.F.2.m)
        ex.sample.temp.F.2.h[,time.var_names[i]] <- predict_function(model1=models_females.h[[i]], dataset=ex.sample.temp.F.2.h)
        
        ex.sample.temp.M.2.l[,time.var_names[i]] <- predict_function(model1=models_males.l[[i]], dataset=ex.sample.temp.M.2.l)
        ex.sample.temp.M.2.m[,time.var_names[i]] <- predict_function(model1=models_males.m[[i]], dataset=ex.sample.temp.M.2.m)
        ex.sample.temp.M.2.h[,time.var_names[i]] <- predict_function(model1=models_males.h[[i]], dataset=ex.sample.temp.M.2.h)
      }
      
      if (intervention=="LBRF"){
        ex.sample.temp.F.2.l$ind_earn_adj <- predict_function(model1=models_females.l$ind_earn_adj, dataset=ex.sample.temp.F.2.l,
                                                              quantiles=quant)
        ex.sample.temp.F.2.m$ind_earn_adj <- predict_function(model1=models_females.m$ind_earn_adj, dataset=ex.sample.temp.F.2.m,
                                                              quantiles=quant)
        ex.sample.temp.F.2.h$ind_earn_adj <- predict_function(model1=models_females.h$ind_earn_adj, dataset=ex.sample.temp.F.2.h,
                                                              quantiles=quant)
        
        ex.sample.temp.M.2.l$ind_earn_adj <- predict_function(model1=models_males.l$ind_earn_adj, dataset=ex.sample.temp.M.2.l,
                                                              quantiles=quant)
        ex.sample.temp.M.2.m$ind_earn_adj <- predict_function(model1=models_males.m$ind_earn_adj, dataset=ex.sample.temp.M.2.m,
                                                              quantiles=quant)
        ex.sample.temp.M.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.2.h,
                                                              quantiles=quant)
        
        #intervene here so that l.LBRF is changed
        # take the observed distribution of employment in males and randomly assign it to females
        ex.sample.temp.F.2.l$LBRF <- predict_function(models_males.l$LBRF.model.M, ex.sample.temp.F.2.l)
        ex.sample.temp.F.2.m$LBRF <- predict_function(models_males.m$LBRF.model.M, ex.sample.temp.F.2.m)
        ex.sample.temp.F.2.h$LBRF <- predict_function(models_males.h$LBRF.model.M, ex.sample.temp.F.2.h)
        
        
        ## predict occupation
        ex.sample.temp.F.2.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.l$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.l$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.l$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_females.l$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2.l))))
        ex.sample.temp.F.2.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.m$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.m$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.m$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_females.m$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2.m))))
        
        ex.sample.temp.F.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_females.h$occupation_class_4cat.model.F, dataset=ex.sample.temp.F.2.h))))
        
        
        ex.sample.temp.M.2.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.l$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.l$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.l$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.l$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.l))))
        ex.sample.temp.M.2.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.m$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.m$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.m$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.m$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.m))))
        
        ex.sample.temp.M.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.h))))
        
        
        
      } else if (intervention=="LBRF+OCC"){
        ex.sample.temp.F.2.l$ind_earn_adj <- predict_function(model1=models_females.l$ind_earn_adj, dataset=ex.sample.temp.F.2.l,
                                                              quantiles=quant)
        ex.sample.temp.F.2.m$ind_earn_adj <- predict_function(model1=models_females.m$ind_earn_adj, dataset=ex.sample.temp.F.2.m,
                                                              quantiles=quant)
        ex.sample.temp.F.2.h$ind_earn_adj <- predict_function(model1=models_females.h$ind_earn_adj, dataset=ex.sample.temp.F.2.h,
                                                              quantiles=quant)
        
        ex.sample.temp.M.2.l$ind_earn_adj <- predict_function(model1=models_males.l$ind_earn_adj, dataset=ex.sample.temp.M.2.l,
                                                              quantiles=quant)
        ex.sample.temp.M.2.m$ind_earn_adj <- predict_function(model1=models_males.m$ind_earn_adj, dataset=ex.sample.temp.M.2.m,
                                                              quantiles=quant)
        ex.sample.temp.M.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.2.h,
                                                              quantiles=quant)
        
        #intervene here so that l.LBRF is changed
        # take the observed distribution of employment in males and randomly assign it to females
        ex.sample.temp.F.2.l$LBRF <- predict_function(models_males.l$LBRF.model.M, ex.sample.temp.F.2.l)
        ex.sample.temp.F.2.m$LBRF <- predict_function(models_males.m$LBRF.model.M, ex.sample.temp.F.2.m)
        ex.sample.temp.F.2.h$LBRF <- predict_function(models_males.h$LBRF.model.M, ex.sample.temp.F.2.h)
        
        
        ## predict occupation
        ex.sample.temp.F.2.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.l$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.l$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.l$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.l$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.l))))
        ex.sample.temp.F.2.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.m$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.m$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.m$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.m$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.m))))
        
        ex.sample.temp.F.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.h))))
        
        
        ex.sample.temp.M.2.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.l$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.l$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.l$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.l$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.l))))
        ex.sample.temp.M.2.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.m$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.m$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.m$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.m$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.m))))
        
        ex.sample.temp.M.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.h))))
        
        
      } else if (intervention=="LBRF+OCC+INCOME" | intervention=="ALL + PAST INCOME"){
        ex.sample.temp.F.2.l$ind_earn_adj <- predict_function(model1=models_males.l$ind_earn_adj, dataset=ex.sample.temp.F.2.l,
                                                              quantiles=quant)
        ex.sample.temp.F.2.m$ind_earn_adj <- predict_function(model1=models_males.m$ind_earn_adj, dataset=ex.sample.temp.F.2.m,
                                                              quantiles=quant)
        ex.sample.temp.F.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.F.2.h,
                                                              quantiles=quant)
        
        ex.sample.temp.M.2.l$ind_earn_adj <- predict_function(model1=models_males.l$ind_earn_adj, dataset=ex.sample.temp.M.2.l,
                                                              quantiles=quant)
        ex.sample.temp.M.2.m$ind_earn_adj <- predict_function(model1=models_males.m$ind_earn_adj, dataset=ex.sample.temp.M.2.m,
                                                              quantiles=quant)
        ex.sample.temp.M.2.h$ind_earn_adj <- predict_function(model1=models_males.h$ind_earn_adj, dataset=ex.sample.temp.M.2.h,
                                                              quantiles=quant)
        
        #intervene here so that l.LBRF is changed
        # take the observed distribution of employment in males and randomly assign it to females
        ex.sample.temp.F.2.l$LBRF <- predict_function(models_males.l$LBRF.model.M, ex.sample.temp.F.2.l)
        ex.sample.temp.F.2.m$LBRF <- predict_function(models_males.m$LBRF.model.M, ex.sample.temp.F.2.m)
        ex.sample.temp.F.2.h$LBRF <- predict_function(models_males.h$LBRF.model.M, ex.sample.temp.F.2.h)
        
        
        ## predict occupation
        ex.sample.temp.F.2.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.l$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.l$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.l$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.l$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.l))))
        ex.sample.temp.F.2.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.m$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.m$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.m$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.m$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.m))))
        
        ex.sample.temp.F.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.F.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.F.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.F.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.F.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.F.2.h))))
        
        
        ex.sample.temp.M.2.l$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.l$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.l$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.l$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.l$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.l$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.l))))
        ex.sample.temp.M.2.m$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.m$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.m$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.m$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.m$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.m$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.m))))
        
        ex.sample.temp.M.2.h$occupation_class_4cat <- as.factor(ifelse(ex.sample.temp.M.2.h$LBRF == "3.Unemployed"|
                                                                         ex.sample.temp.M.2.h$LBRF == "5.Retired"|
                                                                         ex.sample.temp.M.2.h$LBRF == "6.Disabled"|
                                                                         ex.sample.temp.M.2.h$LBRF == "7.Not in LbrF", "4.doesnt work",
                                                                       as.character(predict_function(model1=models_males.h$occupation_class_4cat.model.M, dataset=ex.sample.temp.M.2.h))))
        
        
      } else {return("Intervention is not correctly specified!")}
      
      
      #white
      # join this newly created data with the old data
      ex.sample.F.2.l <- rbind(ex.sample.F.2.l,ex.sample.temp.F.2.l)
      rm(ex.sample.temp.F.2.l)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.F.2.l <- ex.sample.F.2.l[order(ex.sample.F.2.l$idnr,ex.sample.F.2.l$t),]
      
      #hisp
      ex.sample.F.2.m <- rbind(ex.sample.F.2.m,ex.sample.temp.F.2.m)
      rm(ex.sample.temp.F.2.m)
      
      ex.sample.F.2.m <- ex.sample.F.2.m[order(ex.sample.F.2.m$idnr,ex.sample.F.2.m$t),]
      
      #black
      ex.sample.F.2.h <- rbind(ex.sample.F.2.h,ex.sample.temp.F.2.h)
      rm(ex.sample.temp.F.2.h)
      
      ex.sample.F.2.h <- ex.sample.F.2.h[order(ex.sample.F.2.h$idnr,ex.sample.F.2.h$t),]
      
      
      #males
      #white
      # join this newly created data with the old data
      ex.sample.M.2.l <- rbind(ex.sample.M.2.l,ex.sample.temp.M.2.l)
      rm(ex.sample.temp.M.2.l)
      
      # order by ID variable - not necessary, but make sure it doesnt mess anything up!!!
      ex.sample.M.2.l <- ex.sample.M.2.l[order(ex.sample.M.2.l$idnr,ex.sample.M.2.l$t),]
      
      #hisp
      ex.sample.M.2.m <- rbind(ex.sample.M.2.m,ex.sample.temp.M.2.m)
      rm(ex.sample.temp.M.2.m)
      
      ex.sample.M.2.m <- ex.sample.M.2.m[order(ex.sample.M.2.m$idnr,ex.sample.M.2.m$t),]
      
      #black
      ex.sample.M.2.h <- rbind(ex.sample.M.2.h,ex.sample.temp.M.2.h)
      rm(ex.sample.temp.M.2.h)
      
      ex.sample.M.2.h <- ex.sample.M.2.h[order(ex.sample.M.2.h$idnr,ex.sample.M.2.h$t),]
      
      
    }
    
    # save mc results
    # before we save our results, we will subset != retired
    # we do this because we are not interested in the fact of retirement on this but the effect of the other
    # lbrf variables!
    
    ex.sample.F.l.noRET <- ex.sample.F.l[ex.sample.F.l$LBRF!="5.Retired" & ex.sample.F.l$LBRF!="6.Disabled",]
    ex.sample.F.m.noRET <- ex.sample.F.m[ex.sample.F.m$LBRF!="5.Retired" & ex.sample.F.m$LBRF!="6.Disabled",]
    ex.sample.F.h.noRET <- ex.sample.F.h[ex.sample.F.h$LBRF!="5.Retired" & ex.sample.F.h$LBRF!="6.Disabled",]
    
    ex.sample.M.l.noRET <- ex.sample.M.l[ex.sample.M.l$LBRF!="5.Retired" & ex.sample.M.l$LBRF!="6.Disabled",]
    ex.sample.M.m.noRET <- ex.sample.M.m[ex.sample.M.m$LBRF!="5.Retired" & ex.sample.M.m$LBRF!="6.Disabled",]
    ex.sample.M.h.noRET <- ex.sample.M.h[ex.sample.M.h$LBRF!="5.Retired" & ex.sample.M.h$LBRF!="6.Disabled",]
    
    ex.sample.F.2.l.noRET <- ex.sample.F.2.l[ex.sample.F.2.l$LBRF!="5.Retired" & ex.sample.F.2.l$LBRF!="6.Disabled",]
    ex.sample.F.2.m.noRET <- ex.sample.F.2.m[ex.sample.F.2.m$LBRF!="5.Retired" & ex.sample.F.2.m$LBRF!="6.Disabled",]
    ex.sample.F.2.h.noRET <- ex.sample.F.2.h[ex.sample.F.2.h$LBRF!="5.Retired" & ex.sample.F.2.h$LBRF!="6.Disabled",]
    
    ex.sample.M.2.l.noRET <- ex.sample.M.2.l[ex.sample.M.2.l$LBRF!="5.Retired" & ex.sample.M.2.l$LBRF!="6.Disabled",]
    ex.sample.M.2.m.noRET <- ex.sample.M.2.m[ex.sample.M.2.m$LBRF!="5.Retired" & ex.sample.M.2.m$LBRF!="6.Disabled",]
    ex.sample.M.2.h.noRET <- ex.sample.M.2.h[ex.sample.M.2.h$LBRF!="5.Retired" & ex.sample.M.2.h$LBRF!="6.Disabled",]
    
    ## aggregate by t ##
    for (i in 1:length(time.var_names)){
      if(time.var_names[i] == "ind_earn_adj"){
        tab.F.l <- ex.sample.F.l.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.l[[i]][,m+1] <- tab.F.l$mean
        tab.F.m <- ex.sample.F.m.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.m[[i]][,m+1] <- tab.F.m$mean
        tab.F.h <- ex.sample.F.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.h[[i]][,m+1] <- tab.F.h$mean
        
        tab.M.l <- ex.sample.M.l.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.l[[i]][,m+1] <- tab.M.l$mean
        tab.M.m <- ex.sample.M.m.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.m[[i]][,m+1] <- tab.M.m$mean
        tab.M.h <- ex.sample.M.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.h[[i]][,m+1] <- tab.M.h$mean
        
        
        tab.F2.l <- ex.sample.F.2.l.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.l.2[[i]][,m+1] <- tab.F2.l$mean
        tab.F2.m <- ex.sample.F.2.m.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.m.2[[i]][,m+1] <- tab.F2.m$mean
        tab.F2.h <- ex.sample.F.2.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.F.h.2[[i]][,m+1] <- tab.F2.h$mean
        
        tab.M2.l <- ex.sample.M.2.l.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.l.2[[i]][,m+1] <- tab.M2.l$mean
        tab.M2.m <- ex.sample.M.2.m.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.m.2[[i]][,m+1] <- tab.M2.m$mean
        tab.M2.h <- ex.sample.M.2.h.noRET %>% group_by(t) %>% summarise(mean = mean(ind_earn_adj,na.rm=T))
        mc.results.arr.M.h.2[[i]][,m+1] <- tab.M2.h$mean
        
      } else if(length(levels(dataframe[,time.var_names[i]]))==2){
        mc.results.arr.F.l[[i]][,m+1] <- prop.table(table(ex.sample.F.l.noRET$t, ex.sample.F.l.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.m[[i]][,m+1] <- prop.table(table(ex.sample.F.m.noRET$t, ex.sample.F.m.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.h[[i]][,m+1] <- prop.table(table(ex.sample.F.h.noRET$t, ex.sample.F.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.M.l[[i]][,m+1] <- prop.table(table(ex.sample.M.l.noRET$t, ex.sample.M.l.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.m[[i]][,m+1] <- prop.table(table(ex.sample.M.m.noRET$t, ex.sample.M.m.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.h[[i]][,m+1] <- prop.table(table(ex.sample.M.h.noRET$t, ex.sample.M.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.F.l.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.l.noRET$t, ex.sample.F.2.l.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.m.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.m.noRET$t, ex.sample.F.2.m.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.F.h.2[[i]][,m+1] <- prop.table(table(ex.sample.F.2.h.noRET$t, ex.sample.F.2.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
        mc.results.arr.M.l.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.l.noRET$t, ex.sample.M.2.l.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.m.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.m.noRET$t, ex.sample.M.2.m.noRET[,time.var_names[i]]),margin=1)[,2]
        mc.results.arr.M.h.2[[i]][,m+1] <- prop.table(table(ex.sample.M.2.h.noRET$t, ex.sample.M.2.h.noRET[,time.var_names[i]]),margin=1)[,2]
        
      } else {
        mc.results.arr.F.l[[i]][,m+1] <- t(prop.table(table(ex.sample.F.l.noRET$t, ex.sample.F.l.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.m[[i]][,m+1] <- t(prop.table(table(ex.sample.F.m.noRET$t, ex.sample.F.m.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.h[[i]][,m+1] <- t(prop.table(table(ex.sample.F.h.noRET$t, ex.sample.F.h.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.M.l[[i]][,m+1] <- t(prop.table(table(ex.sample.M.l.noRET$t, ex.sample.M.l.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.m[[i]][,m+1] <- t(prop.table(table(ex.sample.M.m.noRET$t, ex.sample.M.m.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.h[[i]][,m+1] <- t(prop.table(table(ex.sample.M.h.noRET$t, ex.sample.M.h.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.F.l.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.l.noRET$t, ex.sample.F.2.l.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.m.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.m.noRET$t, ex.sample.F.2.m.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.F.h.2[[i]][,m+1] <- t(prop.table(table(ex.sample.F.2.h.noRET$t, ex.sample.F.2.h.noRET[,time.var_names[i]]),margin=1))
        
        mc.results.arr.M.l.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.l.noRET$t, ex.sample.M.2.l.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.m.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.m.noRET$t, ex.sample.M.2.m.noRET[,time.var_names[i]]),margin=1))
        mc.results.arr.M.h.2[[i]][,m+1] <- t(prop.table(table(ex.sample.M.2.h.noRET$t, ex.sample.M.2.h.noRET[,time.var_names[i]]),margin=1))
        
      }
    }
  }
  #need to return all of these in one list ideally!
  results <- list(mc.results.arr.F.l,mc.results.arr.F.m,mc.results.arr.F.h,
                  mc.results.arr.M.l,mc.results.arr.M.m,mc.results.arr.M.h,
                  mc.results.arr.F.l.2,mc.results.arr.F.m.2,mc.results.arr.F.h.2,
                  mc.results.arr.M.l.2,mc.results.arr.M.m.2,mc.results.arr.M.h.2)
  names(results) <- c("mc.results.arr.F.l","mc.results.arr.F.m","mc.results.arr.F.h",
                      "mc.results.arr.M.l","mc.results.arr.M.m","mc.results.arr.M.h",
                      "mc.results.arr.F.l.2","mc.results.arr.F.m.2","mc.results.arr.F.h.2",
                      "mc.results.arr.M.l.2","mc.results.arr.M.m.2","mc.results.arr.M.h.2")
  
  return(results)
}
