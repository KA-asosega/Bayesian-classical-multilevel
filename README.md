# Bayesian-classical-multilevel
Malaria among U5 analysis code
R -code for Data analysis
maladata=read.csv("…….\\Malaria_Analysis1.csv", header = T)
maladata$Age_Mother<-as.numeric(as.character.Date( maladata$Age_Mother))
maladata$Type_Residence<-as.factor(maladata$Type_Residence)
maladata$Wealth_StatusNew <- as.factor(maladata$Wealth_StatusNew)
maladata$M_net <- as.factor(maladata$M_net)
maladata$NHIS_coverage <-as.factor((maladata$NHIS_coverage))
maladata$Child_Sex <- as.factor(maladata$Child_Sex)
maladata$Child_Slept_net <- as.factor(maladata$Child_Slept_net)
maladata$Anaemia_Stat  <- as.factor(maladata$Anaemia_Stat)
maladata$Age_Group <- as.factor(maladata$Age_Group)
maladata$Ethnicity<-as.factor(maladata$Ethnicity)
str(maladata)

# Accounting for survey design
library(survey)
library(haven)
library(dplyr)
mysurvey<- svydesign(id=maladata$Cluster.Number, data = maladata, strata = maladata$strata, weights = maladata$weight, nest = T)
                     options(survey.lonely.psu = "adjust")
svytable(~ Malaria, mysurvey)
malaria_P= prop.table(svytable(~ Malaria, mysurvey))
binom.test(859,2895, conf.level = 0.95)
round(malaria_P,2)
table(maladata$Malaria)
prop.table(table(maladata$Malaria))
#crosstabulation of malaria with residence type
svyby(~Malaria, by= ~Child_Sex, design = mysurvey, FUN = svymean, vartype = c("se","ci"))
#malaria and age group
str(mysurvey)
svyby(~Malaria, by= ~ Age_Group, design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Ethnicity, design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Type_Residence , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Wealth_StatusNew , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Child_Sex  , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Child_Slept_net  , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Anaemia_Stat  , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ Edu_Stat  , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
svyby(~Malaria, by= ~ M_net    , design = mysurvey, FUN = svymean, vartype = c("se","ci"))
#Chisquare test of association
svychisq(~Malaria +Age_Group, mysurvey)
svychisq(~Malaria +Ethnicity, mysurvey)
svychisq(~Malaria + Type_Residence, mysurvey)
svychisq(~Malaria + Wealth_StatusNew, mysurvey)
svychisq(~Malaria + Child_Sex, mysurvey)
svychisq(~Malaria + Child_Slept_net, mysurvey)
svychisq(~Malaria + Anaemia_Stat, mysurvey)
svychisq(~Malaria + Edu_Stat, mysurvey)
svychisq(~Malaria + M_net, mysurvey)


#Bayesian Multilevel Model
library(brms)
library(tidybayes)
library(tidyverse)
library(devtools)
library(ROCR)
library(posterior)
library(HDInterval)
library(rstan)
library(rstanarm)
library(rstantools)
library(ggplot2)
library(dplyr)
library(merTools)
library(lme4)
library(jtools)
library(performance)
# Bayesian Random Intercept Only Model
m0 = brm(Malaria~ 1 +(1|cluster_NO.),data = maladata, family = bernoulli(link = "logit"))
summary(m0)
prior_summary(m0)
loo(m0)
model_performance(m0)

# Final Bayesian Model
prior1<- c(prior(normal(0,10), class= Intercept),
           prior(cauchy(0,5), class= sd), 
           prior(normal(0,5), class= b))
final_model<- brm( Malaria~ 1+ Age_Group + Wealth_StatusNew + Anaemia_Stat + Child_Slept_net +Ethnicity + BUILT_Population +  LS_Temperature + Log_rain +(1|Cluster.Number),data = maladata, family = bernoulli(link = "logit"), prior = prior1)
summary(final_model)
loo(final_model)

prior_summary(final_model)
# Assessing for multidisciplinary among predictors 
check_collinearity(final_model)

plot(final_model)
# Posterior predictive plot
pp_check(final_model, nsamples = 200 + theme(length.position= c(0.5, 0.5)))

conditional_effects(final_model)
# Trace and Density plots
final_model %>%
  plot(combo=c("hist","trace"), width=c(1,1.5), theme=theme_bw(base_size = 10))
library(shinystan)
launch_shinystan(final_model)

#Classical Method
m1<- glmer(Malaria~ 1 +(1|cluster_NO.),data = maladata, family = binomial )
summary(m1)
summ(m1)
#Final Classical model
Cl_model<- glmer(Malaria~ 1+ Age_Group + Wealth_StatusNew +Anaemia_Stat  + Child_Slept_net +Ethnicity +BUILT_Population +  LS_Temperature + Log_rain +(1|cluster_NO.),data = maladata, family = binomial,nAGQ = 10, control = glmerControl(optimizer = "bobyqa") )
summary(Cl_model)
plot(Cl_model)


#Cluster Level Covariates
maU5=read.csv("C: ….\\Cluster Covariates\\Cluster_level1.csv",header = TRUE)
str(maU5)
library(PerformanceAnalytics)
library(ggcorrplot)
ggcorrplot(cor(maU5),
           hc.order = TRUE,
           type = "lower",
           lab = TRUE)
chart.Correlation(maU5, histogram = TRUE, method = "pearson")
cor(maU5)
ggplot(cor(maU5))

#Descriptive Statistics for Cluster level covariates
library(psych)
describe(maU5)
