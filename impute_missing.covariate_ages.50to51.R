library(mice) # multiple imputation with chained equations

load("df_inflation.adj.RData")
df_gform <- df_inflation.adj

## IMPUTE MISSING DATA FOR THE STTARTING POINT OF THE SIMULATION (AGES 50-51)
# impute outside of bootstrap with enough iterations and m ----------------

df_gform.t1 <- df_gform[df_gform$t==1,]
k <- length(names(df_gform.t1))
predictor_matrix <- matrix(0,nrow = k, ncol = k)
colnames(predictor_matrix) <- rownames(predictor_matrix) <- names(df_gform)

#Now we specific which variables we want to impute, by creating a vector of same length of the number of total
#variables, and putting in the imputation method for the columns (variables) we want to impute
imp_meth <- rep("",k)

imp_meth[c(3,4,8)] <- "logreg"  # binomial
imp_meth[c(6,9,12,17)] <- "polyreg"# multinomial
imp_meth[c(7,11,13,14,15,16)] <- "polr"# multinomial,ordered
imp_meth[c(5,31)] <- "norm"  # numeric
covariate.index <- c(3:9,11:17,31)

#Lastly, we specify which variables we want to use as predictors for our variables to be imputed. Note that you want to use
#the other variables to be imputed as a predictors as well, this is the basis of mice.
#using variables at t (and t-1) is richer than our estimation model - will go with that

covarnames <- colnames(df_gform.t1[,covariate.index])

for(i in 1:length(names(df_gform.t1)[covariate.index])) {
  
  varname <- names(df_gform.t1)[covariate.index[i]]
  
  predictor_matrix[varname,covarnames[-i]] <- 1
  # -i because a variable shouldn't be a predictor of itself
  # this is the only way in which the current approach is a bit less congenial potentially
  
}

set.seed(1234)

impute <- mice(df_gform.t1, method = imp_meth, predictorMatrix = predictor_matrix, m = 100)
#summary(impute)
#impute$method
#impute$loggedEvents # occ class - doesnt work omitted,
# imputations are still pretty accurate though (on average)

df_gform.t1_imp <- complete(impute)

#save 
save(df_gform.t1_imp, file="df_gform.t1_imp.RData")
