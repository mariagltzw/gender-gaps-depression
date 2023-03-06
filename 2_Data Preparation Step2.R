
#########################################
########### DATA HANDLING P.2 FOR ####### 
########### GENDER GAPS PAPER    ######## 
########### GUELTZOW ET AL       ######## 
#########################################

library(dplyr)
## contains a lag function that we need later ##
source("functions.R")

# EXCLUDE NON ELIGIBLE RESPONDENTS AND RECODE SOME VARIABLES ----------------------------------------------------------------------
## load previously created dataset ##
load("df_resp_long.RData")

## include everyone below the age of 50 and above the age of 80 ##
## only keep respondents that respondents and are alive in specific wave ## 
df_elig_resp <- df_resp_long[df_resp_long$AGEY_E >= 50 & df_resp_long$AGEY_E <=80 & df_resp_long$response_stat=="1.Resp,alive",]

## exclude respondents with missing values on the depressive symptoms measure ##
df_elig_resp <- df_elig_resp[!is.na(df_elig_resp$CESD),]

## create variable for elevated depressive symptoms cut-offs, 1 = elevated depressive symptoms, 0 = no elevated depressive symptoms ## 
df_elig_resp$depress <- ifelse(df_elig_resp$CESD >= "3",1, 
                               ifelse(df_elig_resp$CESD < "3",0, NA))

#### recode baseline confounders####
## recode gender ##
df_elig_resp$RAGENDER <- factor(df_elig_resp$RAGENDER, levels=c("2.Female","1.Male"), labels=c("female","male"))

## combine race and ethnicity variables into one variable ##
## so we have non-hispanic white, non-hispanic black, black, other ##
df_elig_resp$ethn <- ifelse(df_elig_resp$RARACEM == "1.White/Caucasian" & df_elig_resp$RAHISPAN == "0.Not Hispanic" | 
                              df_elig_resp$RARACEM == "1.White/Caucasian" & is.na(df_elig_resp$RAHISPAN), "1. non-Hisp White",
                            ifelse(df_elig_resp$RARACEM == "1.White/Caucasian" & df_elig_resp$RAHISPAN == "1.Hispanic", "3.Hispanic",
                                   ifelse(df_elig_resp$RARACEM == "2.Black/African American" & df_elig_resp$RAHISPAN == "0.Not Hispanic" |
                                            df_elig_resp$RARACEM == "2.Black/African American" & is.na(df_elig_resp$RAHISPAN), "2. non-Hisp Black",
                                          ifelse(df_elig_resp$RARACEM == "2.Black/African American" & df_elig_resp$RAHISPAN == "1.Hispanic", "3.Hispanic",
                                                 ifelse(df_elig_resp$RARACEM == "3.Other" & df_elig_resp$RAHISPAN == "1.Hispanic", "3.Hispanic", 
                                                        ifelse(df_elig_resp$RARACEM == "3.Other" & df_elig_resp$RAHISPAN == "0.Not Hispanic"|
                                                                 df_elig_resp$RARACEM == "3.Other" & is.na(df_elig_resp$RAHISPAN), "4.Other",NA))))))
df_elig_resp$ethn <- as.factor(df_elig_resp$ethn)

## recode mother's and father's education into low middle and high
## low: 0-8, highschool/middle: 9-12, college/high: 13-16, 17+
df_elig_resp$RAFEDUC <- factor(df_elig_resp$RAFEDUC, labels = c(0,1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,7.5,8,8.5,9))
df_elig_resp$RAFEDUC <- as.numeric(as.character(df_elig_resp$RAFEDUC))

df_elig_resp$RAFEDUC_cat <- as.factor(ifelse(df_elig_resp$RAFEDUC <9, "1.low",
                                   ifelse(df_elig_resp$RAFEDUC >=9 & df_elig_resp$RAFEDUC<13, "2.middle",
                                          ifelse(df_elig_resp$RAFEDUC >=13, "3.high", NA))))

#table(df_elig_resp$RAFEDUC_cat)

df_elig_resp$RAMEDUC <- factor(df_elig_resp$RAMEDUC, labels = c(0,1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,7.5,8,8.5,9))
df_elig_resp$RAMEDUC <- as.numeric(as.character(df_elig_resp$RAMEDUC))

df_elig_resp$RAMEDUC_cat <- as.factor(ifelse(df_elig_resp$RAMEDUC <9, "1.low",
                                   ifelse(df_elig_resp$RAMEDUC >=9 & df_elig_resp$RAMEDUC<13, "2.middle",
                                          ifelse(df_elig_resp$RAMEDUC >=13, "3.high", NA))))
#table(df_elig_resp$RAMEDUC_cat)

## categorize education (of respondent) into 3 levels ##
## less than hig school/GED; high school diploma; and some college or higher
#table(df_elig_resp$RAEDUC)

df_elig_resp$education <- ifelse(df_elig_resp$RAEDUC == "1.Lt High-school" | df_elig_resp$RAEDUC == "2.GED", "Lt High-school/GED",
                                 ifelse(df_elig_resp$RAEDUC == "3.High-school graduate", "High-school graduate",
                                        ifelse(df_elig_resp$RAEDUC == "4.Some college" | df_elig_resp$RAEDUC == "5.College and above", "Some college or higher",
                                               NA)))
df_elig_resp$education <- as.factor(df_elig_resp$education)



#### recode time-varing covariates ####

## recode occupation codes ##
## there are a lot of missings but i think that might also be that variable is set to NA for retired or homemakers / anyone with no ind earn ##
## special missing codes in SAS data (codebook) show that lots have missing because they are not working ##
## let's create these by outselves ##
df_elig_resp$occupation_code_2 <- ifelse(df_elig_resp$LBRF == "3.Unemployed"|df_elig_resp$LBRF == "5.Retired"|
                                           df_elig_resp$LBRF == "6.Disabled"|df_elig_resp$LBRF == "7.Not in LbrF","18",
                                 df_elig_resp$occupation_code)
df_elig_resp$occupation_code_2 <- factor(df_elig_resp$occupation_code_2, levels = c("1","2","3","4","5","6",
                                                                "7","8","9",
                                                                "10","11","12","13","14","15","16","17","18"),
                               labels = c("01.Managerial specialty oper",   "02.Prof specialty opr/tech sup",
                                          "03.Sales",                       "04.Clerical/admin supp"        ,
                                          "05.Svc:prv hhld/clean/bldg svc", "06.Svc:protection"             ,
                                          "07.Svc:food prep" ,              "08.Health svc"                 ,
                                          "09.Personal svc" ,               "10.Farming/forestry/fishing"   ,
                                          "11.Mechanics/repair"  ,          "12.Constr trade/extractors"    ,
                                          "13.Precision production"  ,      "14.Operators: machine"         ,
                                          "15.Operators: transport, etc",   "16.Operators: handlers, etc"   ,
                                          "17.Member of Armed Forces"  ,"18.doesnt work"))

## classify occupation into blue collar (manual), white collar (desk work), pink collar (service related) ##
## assigned army to blue collar even though it could be any of the 3 groups, however they make up only 16 observations so it should not make a big difference ##
df_elig_resp$occupation_class_4cat <- as.factor(ifelse(df_elig_resp$occupation_code_2 == "01.Managerial specialty oper"|df_elig_resp$occupation_code_2 == "02.Prof specialty opr/tech sup"|
                                                         df_elig_resp$occupation_code_2 == "03.Sales"| df_elig_resp$occupation_code_2 == "04.Clerical/admin supp", "1.White collar/desk job",
                                             ifelse(df_elig_resp$occupation_code_2 == "05.Svc:prv hhld/clean/bldg svc"|df_elig_resp$occupation_code_2 == "06.Svc:protection"|
                                                      df_elig_resp$occupation_code_2 == "07.Svc:food prep"| df_elig_resp$occupation_code_2 == "08.Health svc"| df_elig_resp$occupation_code_2 =="09.Personal svc", "2.Pink collar/Service-related job",
                                                    ifelse(df_elig_resp$occupation_code_2 == "10.Farming/forestry/fishing"|df_elig_resp$occupation_code_2 == "11.Mechanics/repair"|
                                                             df_elig_resp$occupation_code_2 == "12.Constr trade/extractors"| df_elig_resp$occupation_code_2 == "13.Precision production"| df_elig_resp$occupation_code_2 =="14.Operators: machine" |
                                                             df_elig_resp$occupation_code_2 =="14.Operators: machine"| df_elig_resp$occupation_code_2 =="15.Operators: transport, etc"|
                                                             df_elig_resp$occupation_code_2 =="16.Operators: handlers, etc"| df_elig_resp$occupation_code_2 =="17.Member of Armed Forces", "3.Blue-collar/Manual job",
                                                           ifelse(is.na(df_elig_resp$occupation_code_2),NA, "4.doesnt work" )))))
## create one that has NAs if they are not working ##
## we will need this for the analysis later ## 
df_elig_resp$occupation_class_3cat <- as.factor(ifelse(df_elig_resp$occupation_code == "01.Managerial specialty oper"|df_elig_resp$occupation_code == "02.Prof specialty opr/tech sup"|
                                                         df_elig_resp$occupation_code == "03.Sales"| df_elig_resp$occupation_code == "04.Clerical/admin supp", "1.White collar/desk job",
                                             ifelse(df_elig_resp$occupation_code == "05.Svc:prv hhld/clean/bldg svc"|df_elig_resp$occupation_code == "06.Svc:protection"|
                                                      df_elig_resp$occupation_code == "07.Svc:food prep"| df_elig_resp$occupation_code == "08.Health svc"| df_elig_resp$occupation_code =="09.Personal svc", "2.Pink collar/Service-related job",
                                                    ifelse(df_elig_resp$occupation_code == "10.Farming/forestry/fishing"|df_elig_resp$occupation_code == "11.Mechanics/repair"|
                                                             df_elig_resp$occupation_code == "12.Constr trade/extractors"| df_elig_resp$occupation_code == "13.Precision production"| df_elig_resp$occupation_code =="14.Operators: machine" |
                                                             df_elig_resp$occupation_code =="14.Operators: machine"| df_elig_resp$occupation_code =="15.Operators: transport, etc"|
                                                             df_elig_resp$occupation_code =="16.Operators: handlers, etc"| df_elig_resp$occupation_code =="17.Member of Armed Forces", "3.Blue-collar/Manual job",
                                                           NA))))

## create a number of chronic conditions variable ## 
## made from whether has condition y/n, if dispute take the one that they reported this wave ## 
## means that we combined dip prev record and no cond and disp prev record (DK if cond) and preld prob:prev had/no new, not alzheimer's disease, TIA/possible stroke in no condition
## we exclude PSYCH because we assume that it will capture the same information as our outcome measure ##
HIBP  <- ifelse(df_elig_resp$HIBP == "1.Yes" |df_elig_resp$HIBP == "3.Disp prev record and has cond", 1,0)
ARTHR <- ifelse(df_elig_resp$ARTHR == "1.Yes" |df_elig_resp$ARTHR == "3.Disp prev record and has cond", 1, 0)
CANCR <- ifelse(df_elig_resp$CANCR == "1.Yes" |df_elig_resp$CANCR == "3.Disp prev record and has cond", 1, 0)
LUNG  <- ifelse(df_elig_resp$LUNG == "1.Yes" |df_elig_resp$LUNG == "3.Disp prev record and has cond", 1, 0)
STROK <- ifelse(df_elig_resp$STROK == "1.Yes" |df_elig_resp$STROK == "3.Disp prev record and has cond", 1, 0)
PSYCH <- ifelse(df_elig_resp$PSYCH == "1.Yes" |df_elig_resp$PSYCH == "3.Disp prev record and has cond", 1, 0)
HEART <- ifelse(df_elig_resp$HEART == "1.Yes" |df_elig_resp$HEART == "3.Disp prev record and has cond", 1, 0)
DIAB  <- ifelse(df_elig_resp$DIAB == "1.Yes" |df_elig_resp$DIAB == "3.Disp prev record and has cond", 1, 0)

## take sum of conditions
df_elig_resp$N_cond <- HIBP + ARTHR + CANCR + LUNG + STROK  + HEART + DIAB 

## transform N_cond to categorical variable to make the analysis easier
## categorise to get somewhat similar group sizes
df_elig_resp$N_cond_cat <- as.factor(ifelse(df_elig_resp$N_cond == 0,"0.none",
                                  ifelse(df_elig_resp$N_cond == 1,"1",
                                         ifelse(df_elig_resp$N_cond == 2,"2",
                                                ifelse(df_elig_resp$N_cond == 3,"3",
                                                       ifelse(df_elig_resp$N_cond >= 4,"4+",NA))))))
df_elig_resp$N_cond_cat <- factor(df_elig_resp$N_cond_cat, levels=c("0.none","1","2","3","4+"))

## categorise household members into 1, 2, 3, >3 ##
df_elig_resp$H_members <- cut(df_elig_resp$H_members, breaks=c(0,1,2,3,19), labels = c("1","2","3","4.>3"))



df_elig_resp$H_child_cat <- as.factor(ifelse(df_elig_resp$H_child == 0,"0.no children",
                                   ifelse(df_elig_resp$H_child == 1,"1. child",
                                          ifelse(df_elig_resp$H_child == 2,"2. children",
                                                 ifelse(df_elig_resp$H_child > 2,"3.more than 2 children",NA)))))
df_elig_resp$H_child_cat <- factor(df_elig_resp$H_child_cat, levels=c("0.no children","1. child","2. children","3.more than 2 children"))

## collapse marital status into smaller categories
df_elig_resp$marital_stat_cat2 <- as.factor(ifelse(df_elig_resp$marital_stat_nopartner == "1.Married" | df_elig_resp$marital_stat_nopartner == "2.Married,spouse absent","1.Married",
                                         ifelse(df_elig_resp$marital_stat_nopartner == "4.Separated"|
                                                  df_elig_resp$marital_stat_nopartner == "5.Divorced"|
                                                  df_elig_resp$marital_stat_nopartner == "6.Separated/divorced","3.Separated/divorced",
                                                ifelse(df_elig_resp$marital_stat_nopartner == "7.Widowed","4.Widowed",
                                                       ifelse(df_elig_resp$marital_stat_nopartner == "8.Never married"|
                                                                df_elig_resp$marital_stat_nopartner == "9.Unknown unmar","5.Not married",NA)))))



## select variables that I will use for the analysis ##
df_elig_resp_comp <- df_elig_resp[,c("HHIDPN","depress", "wave","RAGENDER", "AGEY_E" , "ethn" , "education",
             "psych_probl_E"  , "LBRF" , "ind_earn" , "N_cond_cat" ,
             "marital_stat_cat2", "H_members" , "H_child_cat")] 
#omit NAs
df_elig_resp_comp <- na.omit(df_elig_resp_comp)

#N respondents
length(unique(df_elig_resp_comp$HHIDPN))

## merge parent's education and occupation code but allow for missing data ##
## we will impute this in a later step ##
df_gform_imp.occ <- merge(df_elig_resp_comp, df_elig_resp[,c("HHIDPN","wave","RAFEDUC_cat" , "RAMEDUC_cat","occupation_class_4cat","occupation_class_3cat")], by=c("HHIDPN","wave"))

## make sure everything has correct class ##
str(df_gform_imp.occ)

# CREATE LAGGED VARIABLES FOR OUTCOME AND TIME-VARYING COVARIATES ----------------------------------------------------------------------
## necessary for our analysis ##
## make sure it's ordered correctly ##
df_gform_imp.occ <- df_gform_imp.occ[
  order( df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave ),
]

df_gform_imp.occ$depress <- as.factor(df_gform_imp.occ$depress)

df_gform_imp.occ$l.depress <- lag1.func(df_gform_imp.occ$depress,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)
df_gform_imp.occ$l.LBRF <- lag1.func(df_gform_imp.occ$LBRF,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)

df_gform_imp.occ$l.ind_earn   <- lag1.func(df_gform_imp.occ$ind_earn,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)
df_gform_imp.occ$l.N_cond_cat <- lag1.func(df_gform_imp.occ$N_cond_cat,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)

df_gform_imp.occ$l.marital_stat_cat <- lag1.func(df_gform_imp.occ$marital_stat_cat,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)
df_gform_imp.occ$l.marital_stat_cat2 <- lag1.func(df_gform_imp.occ$marital_stat_cat2,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)
df_gform_imp.occ$l.H_members <- lag1.func(df_gform_imp.occ$H_members,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)
df_gform_imp.occ$l.H_child_cat <- lag1.func(df_gform_imp.occ$H_child_cat,df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)

## make a tunit variable that combines people age 50+51, 52+53, and so on ##
## we do this because we need to model it in 2 year intervals because of HRS ##
df_gform_imp.occ$t <- cut(df_gform_imp.occ$AGEY_E, include.lowest = T,breaks=c(50,seq(51,81,2)),
            labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))
df_gform_imp.occ$t <- as.numeric(as.character(df_gform_imp.occ$t))

df_gform_imp.occ$l.occupation_class_4cat <- lag1.func(df_gform_imp.occ$occupation_class_4cat, df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)

#lagged occupation class3 for imputation within gform
df_gform_imp.occ$l.occupation_class_3cat <- lag1.func(df_gform_imp.occ$occupation_class_3cat, df_gform_imp.occ$HHIDPN, df_gform_imp.occ$wave)

# ADJUST INCOME FOR INFLATION ----------------------------------------------------------------------
## get inflation in relation to the middle of our time period, which is 2002 ##
## obtained from https://www.bls.gov/data/inflation_calculator.html ##
## taking estimate for June of each year ##
inflation_by_year <- data.frame(wave = unique(df_gform_imp.occ$wave),
                                period =             c("94","96","98","00","02","04","06", "08","10","12","14","16","18"),
                                purch.power_rel.to.02 =c(1215.54,1148.05,1103.68,1043.50,1000,948.34,886.64, 822.16,825.36,783.95,754.79,746.42,713.92))
inflation_by_year$inflation <- inflation_by_year$purch.power_rel.to.02/1000


## now that we have the percentage we just need to calculate the new inflation adjusted income ##
## ref is 2002 so for example 1000 in 2002 would be wort 1100 in 98 ##
df_inflation.adj <- left_join(df_gform_imp.occ, inflation_by_year[,c("wave","inflation")], by="wave")

df_inflation.adj$ind_earn_adj <- df_inflation.adj$ind_earn * df_inflation.adj$inflation

#lag the inflation adjusted income
df_inflation.adj$l.ind_earn_adj <- lag1.func(df_inflation.adj$ind_earn_adj, df_inflation.adj$HHIDPN, df_inflation.adj$wave)

## save new df ##
save(file="df_inflation.adj.RData",df_inflation.adj)

