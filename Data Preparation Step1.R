
#########################################
########### DATA HANDLING FOR ########### 
########### GENDER GAPS PAPER ########### 
########### GUELTZOW ET AL    ########### 
#########################################


# > sessionInfo()
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows Server x64 (build 17763)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252
# [4] LC_NUMERIC=C                    LC_TIME=German_Germany.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.2.1 tools_4.2.1   


#### load packages ####
library("foreign")
library(reshape2)
library(dplyr)
library(tidyselect)
library("lattice")


# EXTRACT RELEVANT VARIABLES FROM HRS DATASET ----------------------------------------------------------------------

## load HRS RAND 1922 to 2018 v1 file
df <- read.spss("randhrs1992_2018v1.sav", to.data.frame = TRUE)


## only keep respondents data
## we are not interested in spouse data 
df_resp   <- df %>%
  dplyr::select(vars_select(names(df), !starts_with('S')))

## find all variables that are constant across time 
idnames   <- (df[,grep("[1,2,3,4,5,6,7,8,9,10,11,12,13,14]", names(df),
                          value=TRUE, invert=TRUE)])

## create a new dataset with variables of interest
base      <- df_resp[,c("HHIDPN","HHID","PN","RABYEAR","RAGENDER", "RAHISPAN", "RARACEM",
                        "RAFEDUC","RAMEDUC", #parents education
                        "RAEDUC", "RAEDYRS", "RAEDEGRM")] #respondents education
## whether proxy interview
proxy <- df_resp[,grep("([P][R][O][X][Y])", names(df_resp), value=TRUE)]
proxy <- proxy[,1:14]

## response status
response_stat <- df_resp[,grep("[I][W][S][T][A][T]", names(df_resp), value=TRUE)] 

## age
age       <- df_resp[,grep("([A][G][E][Y])", names(df_resp), value=TRUE)]
age <- age[,-c(43,44,45)]

#### outcome ####
## depressive symptoms measured with CESD
depr      <- df_resp[,grep("[C][E][S][D]", names(df_resp), value=TRUE)] #R1CESD & R13CESDM is missing
## was not measured in wave 1, add in wave 1 as NA for consistency
depr$R1CESD <- as.numeric(NA)
depr        <- depr[,c(28,1:27)]

#### labor market characteristics ####
## employment variables
LBRF        <- df_resp[,grep("[L][B][R][F]$", names(df_resp), value=TRUE)] #has whether in LBRF and LBRF status

## occupation code
## there are newer occupation codes but this variable is most complete
occupation_code <- df_resp[,grep("[J][C][O][C][C]$", names(df_resp), value=TRUE)] 

## individual earnings: sum of Respondent's wage/salary income, bonuses/overtime pay/commissions/tips,
## 2nd job or military reserve earnings, and professional practice or trade income. 
ind_earn <- df_resp[,grep("[I][E][A][R][N]", names(df_resp), value=TRUE)]


#### time varying covariates ####
## marital status
marital_stat <- df_resp[,grep("([M][S][T][A][T]$)", names(df_resp), value=TRUE)]
marital_stat$REMSTAT <- NULL
marital_stat_nopartner <- df_resp[,grep("([M][S][T][A][T][H])", names(df_resp), value=TRUE)]
marital_stat_nopartner <- marital_stat_nopartner[,1:14]

## Number of chronic conditions, i.e. did get diagnosis this wave 
## alzheimers, diabetes, and sleep is available from wave 10 onwards so i will not include them
HIBP      <- df_resp[,grep("[H][I][B][P]$|[H][I][B][P][Q]", names(df_resp), value=TRUE)] 
HEART      <- df_resp[,grep("[H][E][A][R][T]$|[H][E][A][R][T][Q]", names(df_resp), value=TRUE)] 
ARTHR      <- df_resp[,grep("[A][R][T][H][R]$|[A][R][T][H][R][Q]", names(df_resp), value=TRUE)] 
CANCR      <- df_resp[,grep("[C][A][N][C][R]$|[C][A][N][C][R][Q]", names(df_resp), value=TRUE)] 
DIAB      <- df_resp[,grep("[D][I][A][B]$|[D][I][A][B][Q]", names(df_resp), value=TRUE)] 
LUNG      <- df_resp[,grep("[L][U][N][G]$|[L][U][N][G][Q]", names(df_resp), value=TRUE)] 
STROK      <- df_resp[,grep("[S][T][R][O][K]$|[S][T][R][O][K][Q]", names(df_resp), value=TRUE)] 
PSYCH      <- df_resp[,grep("[P][S][Y][C][H]$|[P][S][Y][C][H][Q]", names(df_resp), value=TRUE)] 
CONDS      <- df_resp[,grep("[C][O][N][D][S]", names(df_resp), value=TRUE)] 
CONDS      <- CONDS[,grep("[C][O][N][D][S][F]", names(CONDS), value=TRUE, invert=TRUE)] 
CONDS$R1CONDSM <- as.numeric(NA)
CONDS$R1CONDSP <- as.numeric(NA)
CONDS$R1CONDS  <- as.numeric(NA)
CONDS      <- CONDS[,c(40,1:13,41,14:26,42,27:39)]
HEALTH <- cbind(HIBP, HEART, ARTHR, CANCR, DIAB, LUNG, STROK, PSYCH,CONDS )
rm(HIBP);rm(HEART);rm(ARTHR);rm(CANCR);rm(DIAB);rm(LUNG);rm(STROK);rm(PSYCH);rm(CONDS)


## ever reported psychological problems? 
psych_probl  <- df_resp[,grep("[P][S][Y][C][H][Q]|[P][S][Y][C][H][E]|[P][S][Y][C][H][S]", names(df_resp), value=TRUE)] 
psych_probl$R1PSYCHS <- as.numeric(NA) 
psych_probl <- psych_probl[,c(1:28,42,29:41)]

## number of people in household
H_members <-  df_resp[,grep("[H][H][R][E][S]$", names(df_resp), value=TRUE)] 

## number of living children (either biological or from partner)
H_child  <- df_resp[,grep("[C][H][I][L][D]$", names(df_resp), value=TRUE)] 



new_dat      <- cbind(base, proxy, response_stat, age, depr, LBRF, occupation_code, ind_earn,
                     marital_stat, marital_stat_nopartner, HEALTH, psych_probl, H_members, H_child) 

rm(proxy);rm(response_stat);rm(age);rm(depr);rm(LBRF);rm(occupation_code);rm(ind_earn);rm(marital_stat);rm(marital_stat_nopartner);rm(HEALTH);rm(psych_probl);rm(H_members);rm(H_child)

## transform data from wide to long ##
df_resp_long    <-  reshape(new_dat,
                direction = "long",
                varying = do.call(list,lapply(seq(13,ncol(new_dat),by=14), #first 13 columns are baseline covariates
                                              function(i) {
                                              (names(new_dat)[i:(i+13)])
                                               })),
                v.names = c("proxy","response_stat","AGEY_B","AGEY_E","AGEY_M","CESD", "CESDM",
                            "LBRF","INLBRF", "occupation_code", "ind_earn",
                            "marital_stat", "marital_stat_nopartner", "HIBPQ", "HIBP", "HEARTQ", "HEART", "ARTHRQ", "ARTHR", "CANCRQ",
                            "CANCR", "DIABQ", "DIAB", "LUNGQ", "LUNG", "STROKQ", "STROK", "PSYCHQ", "PSYCH","CONDSM", "CONDSP", "CONDS", 
                            "psych_probl_Q", "psych_probl_E", "psych_probl_S", "H_members","H_child"),
                idvar = "HHIDPN",
                timevar = "wave",
                times = 1:14)

## save this dataset so it is easier to work with ##
save(df_resp_long, file="df_resp_long.RData")
