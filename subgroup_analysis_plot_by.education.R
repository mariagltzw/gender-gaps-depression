library(dplyr)
library(ggplot2)


# CALCULATE SUMMARY MEASURES ----------------------------------------------
bs=499 #specify bootstrap iterations

#### INTERVENTION A ####
#load output
load("output_educ.subgroup_LBRF.int.Rdata")
results_depr.f  <- as.data.frame(rbind(results_output$natural.course_females_low$depress,
                                       results_output$natural.course_females_middle$depress,
                                       results_output$natural.course_females_high$depress,
                                       results_output$TE_int_females_low$depress, 
                                       results_output$TE_int_females_middle$depress,
                                       results_output$TE_int_females_high$depress))

results_depr.m  <- as.data.frame(rbind(results_output$natural.course_males_low$depress,
                                       results_output$natural.course_males_middle$depress,
                                       results_output$natural.course_males_high$depress,
                                      results_output$TE_int_males_low$depress,
                                      results_output$TE_int_males_middle$depress,
                                      results_output$TE_int_males_high$depress))

# calculate gender depression gap in intervention (cf) and natural course (nc)
abs_gap_cf <- results_depr.f[results_depr.f$V3==1,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]
abs_gap_nc <- results_depr.f[results_depr.f$V3==0,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]

# absolute difference
diff_LBRF <- abs_gap_cf-abs_gap_nc

# aggregate across bootstrap iteration and obtain 95%CIs
diff_LBRF$diff <- apply(diff_LBRF, 1, mean)*100
diff_LBRF$lb <- apply(diff_LBRF, 1, quantile, probs = c(0.025))*100
diff_LBRF$ub <- apply(diff_LBRF, 1, quantile, probs = c(0.975))*100

diff_LBRF$abs_gap_cf <- rowMeans(abs_gap_cf)
diff_LBRF$abs_gap_nc <-rowMeans(abs_gap_nc)

# specify which intervention we introduced
diff_LBRF$effect <- "LbrF"
diff_LBRF$educ <- rep(c("low","middle","high"),each=16)
diff_LBRF$t <- rep(seq(50,80,2),3)

# calculate contribution and obtain 95%CIs
contr <- (1-(abs_gap_cf/abs_gap_nc))

diff_LBRF$contr <- apply(contr, 1, median)*100
diff_LBRF$contr_lb <- apply(contr, 1, quantile, probs = c(0.025))*100
diff_LBRF$contr_ub <- apply(contr, 1, quantile, probs = c(0.975))*100


#### INTERVENTION B ####
load("output_educ.subgroup_LBRF.occ.int.Rdata")

results_depr.f  <- as.data.frame(rbind(results_output$natural.course_females_low$depress,
                                       results_output$natural.course_females_middle$depress,
                                       results_output$natural.course_females_high$depress,
                                       results_output$TE_int_females_low$depress,
                                       results_output$TE_int_females_middle$depress,
                                       results_output$TE_int_females_high$depress))

results_depr.m  <- as.data.frame(rbind(results_output$natural.course_males_low$depress,
                                       results_output$natural.course_males_middle$depress,
                                       results_output$natural.course_males_high$depress,
                                       results_output$TE_int_males_low$depress,
                                       results_output$TE_int_males_middle$depress,
                                       results_output$TE_int_males_high$depress))

abs_gap_cf <- results_depr.f[results_depr.f$V3==1,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]
abs_gap_nc <- results_depr.f[results_depr.f$V3==0,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]


diff_LBRF_occ <- abs_gap_cf-abs_gap_nc

diff_LBRF_occ$diff <- apply(diff_LBRF_occ, 1, mean)*100
diff_LBRF_occ$lb <- apply(diff_LBRF_occ, 1, quantile, probs = c(0.025))*100
diff_LBRF_occ$ub <- apply(diff_LBRF_occ, 1, quantile, probs = c(0.975))*100

diff_LBRF_occ$abs_gap_cf <- rowMeans(abs_gap_cf)
diff_LBRF_occ$abs_gap_nc <-rowMeans(abs_gap_nc)

diff_LBRF_occ$effect <- "LbrF + occupation"
diff_LBRF_occ$educ <- rep(c("low","middle","high"),each=16)
diff_LBRF_occ$t <- rep(seq(50,80,2),3)

contr <- (1-(abs_gap_cf/abs_gap_nc))

diff_LBRF_occ$contr <- apply(contr, 1, median)*100
diff_LBRF_occ$contr_lb <- apply(contr, 1, quantile, probs = c(0.025))*100
diff_LBRF_occ$contr_ub <- apply(contr, 1, quantile, probs = c(0.975))*100


#### INTERVENTION C ####
load("output_educ.subgroup_LBRF.occ.inc.int_.Rdata")

results_depr.f  <- as.data.frame(rbind(results_output$natural.course_females_low$depress,
                                       results_output$natural.course_females_middle$depress,
                                       results_output$natural.course_females_high$depress,
                                       results_output$TE_int_females_low$depress, 
                                       results_output$TE_int_females_middle$depress,
                                       results_output$TE_int_females_high$depress))

results_depr.m  <- as.data.frame(rbind(results_output$natural.course_males_low$depress,
                                       results_output$natural.course_males_middle$depress,
                                       results_output$natural.course_males_high$depress,
                                       results_output$TE_int_males_low$depress, 
                                       results_output$TE_int_males_middle$depress,
                                       results_output$TE_int_males_high$depress))

abs_gap_cf <- results_depr.f[results_depr.f$V3==1,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]
abs_gap_nc <- results_depr.f[results_depr.f$V3==0,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]


diff_LBRF_occ_inc <- abs_gap_cf-abs_gap_nc

diff_LBRF_occ_inc$diff <- apply(diff_LBRF_occ_inc, 1, mean)*100
diff_LBRF_occ_inc$lb <- apply(diff_LBRF_occ_inc, 1, quantile, probs = c(0.025))*100
diff_LBRF_occ_inc$ub <- apply(diff_LBRF_occ_inc, 1, quantile, probs = c(0.975))*100

diff_LBRF_occ_inc$abs_gap_cf <- rowMeans(abs_gap_cf)
diff_LBRF_occ_inc$abs_gap_nc <-rowMeans(abs_gap_nc)

diff_LBRF_occ_inc$effect <- "LbrF + occupation + \nincome"
diff_LBRF_occ_inc$educ <- rep(c("low","middle","high"),each=16)
diff_LBRF_occ_inc$t <- rep(seq(50,80,2),3)

contr <- (1-(abs_gap_cf/abs_gap_nc))

diff_LBRF_occ_inc$contr <- apply(contr, 1, median)*100
diff_LBRF_occ_inc$contr_lb <- apply(contr, 1, quantile, probs = c(0.025))*100
diff_LBRF_occ_inc$contr_ub <- apply(contr, 1, quantile, probs = c(0.975))*100


#### INTERVENTION D ####
load("output_educ.subgroup_LBRF.occ.inc.lag.inc.int.Rdata")

results_depr.f  <- as.data.frame(rbind(results_output$natural.course_females_low$depress,
                                       results_output$natural.course_females_middle$depress,
                                       results_output$natural.course_females_high$depress,
                                       results_output$TE_int_females_low$depress,
                                       results_output$TE_int_females_middle$depress,
                                       results_output$TE_int_females_high$depress))

results_depr.m  <- as.data.frame(rbind(results_output$natural.course_males_low$depress,
                                       results_output$natural.course_males_middle$depress,
                                       results_output$natural.course_males_high$depress,
                                       results_output$TE_int_males_low$depress, 
                                       results_output$TE_int_males_middle$depress,
                                       results_output$TE_int_males_high$depress))

abs_gap_cf <- results_depr.f[results_depr.f$V3==1,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]
abs_gap_nc <- results_depr.f[results_depr.f$V3==0,5:(bs+4)] - results_depr.m[results_depr.m$V3==0,5:(bs+4)]


diff_LBRF_occ_inc2 <- abs_gap_cf-abs_gap_nc

diff_LBRF_occ_inc2$diff <- apply(diff_LBRF_occ_inc2, 1, mean)*100
diff_LBRF_occ_inc2$lb <- apply(diff_LBRF_occ_inc2, 1, quantile, probs = c(0.025))*100
diff_LBRF_occ_inc2$ub <- apply(diff_LBRF_occ_inc2, 1, quantile, probs = c(0.975))*100

diff_LBRF_occ_inc2$abs_gap_cf <- rowMeans(abs_gap_cf)
diff_LBRF_occ_inc2$abs_gap_nc <-rowMeans(abs_gap_nc)

diff_LBRF_occ_inc2$effect <- "LbrF + occupation + \nincome + \nprevious income"
diff_LBRF_occ_inc2$educ <- rep(c("low","middle","high"),each=16)
diff_LBRF_occ_inc2$t <- rep(seq(50,80,2))

contr <- (1-(abs_gap_cf/abs_gap_nc))

diff_LBRF_occ_inc2$contr <- apply(contr, 1, median)*100
diff_LBRF_occ_inc2$contr_lb <- apply(contr, 1, quantile, probs = c(0.025))*100
diff_LBRF_occ_inc2$contr_ub <- apply(contr, 1, quantile, probs = c(0.975))*100

# PLOT INTERVENTIONS TOGETHER ---------------------------------------------
# absolute difference
diff <- rbind(diff_LBRF,diff_LBRF_occ,diff_LBRF_occ_inc,diff_LBRF_occ_inc2)
diff$educ <- factor(diff$educ, levels=c("low","middle","high"))
diff$Age <- factor(diff$t,labels=c("50-51", "52-53","54-55","56-57","58-59", "60-61","62-63","64-65","66-67",
                                   "68-69","70-71","72-73","74-75","76-77","78-79","80"))

#structure results output
results <- data.frame(Age=diff$Age, diff$t, diff$educ, diff$effect,
                      abs_diff=paste0(round(diff$diff,1),"(",round(diff$lb,1),"-",round(diff$ub,1),")"),
                      contr=paste0(round(diff$contr,1),"(",round(diff$contr_lb,1),"-",round(diff$contr_ub,1),")"))

results <- results[results$diff.t<72,]

# assign intervention names as they are used in the paper
diff$effect <- factor(diff$effect, labels=c("A","B","C","D"))
diff$line_width <- ifelse(diff$lb * diff$ub > 0, "a", "b")

cbpalette <- c("#E69F00", "#CC79A7", "#009E73", "#0072B2") 

plot_abs.diff <- diff %>%
  rename(Intervention=effect)%>%
  filter(t<72) %>%
  ggplot(aes(x=Age,y=diff, color=Intervention, group=Intervention, shape=line_width)) +
  geom_line()+ geom_point()+
  scale_shape_manual(values=c("a"=8, "b"=46))+
  scale_color_manual(values=cbpalette)+
  geom_hline(yintercept=0)+
  theme_minimal()+
  guides(shape="none")+
  theme(legend.position = "bottom",
        legend.text = element_text(size=10),
        axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        #axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(size=11),
        strip.text = element_text(size=11)) +
  facet_wrap(~educ, nrow=3) +
  labs(y="%-point Difference",
       x="Age")

#ggsave("Figure3.by.educ.svg", plot_abs.diff, dpi=320)
