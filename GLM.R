######LIBRARIES#####
###lme4:Linear, generalized linear, and nonlinear mixed models
###lme4 provides functions for fitting and analyzing mixed models: linear (lmer), generalized linear (glmer) and nonlinear (nlmer.)

library(lme4)

###car:Companion to Applied Regression
###car includes many functions that require an object created by a modeling function like lm, glm or nls as input.

library(car)

###DHARMa:Residual Diagnostics for HierArchical (Multi-level / Mixed) Regression Models
###DHARMa uses a simulation-based approach to create readily interpretable scaled (quantile) residuals for fitted generalized linear mixed models.

library(DHARMa)

###multcomp: Simultaneous Inference in General Parametric Models
###multcomp includes simultaneous tests and confidence intervals for general linear hypotheses in parametric models, including linear, generalized linear, linear mixed effects, and survival models.

library(multcomp)

###
library(survival)
library(ggplot2)
library(caTools)
library(jtools)
library(Hmisc)
library(plyr)
library(nlme)
library(mgcv)
library(lsmeans)
library(PAmeasures)
library(afex)
library(plotrix)
library(r2glmm)

###

#####DATA#####
###Germination
germ = read.csv(file = "TEMPERATURE_EXPERIMENT_2.csv", header = TRUE)
germ$Temperature <-as.numeric(germ$Temperature)
germ$Cross_Region <-as.factor(germ$Cross_Region)
germ$Mother_Population <-as.factor(germ$Mother_Population)
germ$Father_Population <-as.factor(germ$Father_Population)

g9 = (germ[germ$Cross_Region=="PANXCAR" & germ$Temperature==9,])
g14 = (germ[germ$Cross_Region=="PANXCAR" & germ$Temperature==14,])
g19 = (germ[germ$Cross_Region=="PANXCAR" & germ$Temperature==19,])
g24 = (germ[germ$Cross_Region=="PANXCAR" & germ$Temperature==24,])

summary(g9$Germination)
std.error(g9$Germination)

summary(g14$Germination)
std.error(g14$Germination)

summary(g19$Germination)
std.error(g19$Germination)

summary(g24$Germination)
std.error(g24$Germination)

#germ_pre = germ[germ$Cross=="Conspecific",]

###Survival
surv = germ[germ$Germination=="1",]

s9 =  (surv[surv$Cross_Region=="PANXCAR" & surv$Temperature==9,])
s14 = (surv[surv$Cross_Region=="PANXCAR" & surv$Temperature==14,])
s19 = (surv[surv$Cross_Region=="PANXCAR" & surv$Temperature==19,])
s24 = (surv[surv$Cross_Region=="PANXCAR" & surv$Temperature==24,])

summary(s9$Survival)
std.error(s9$Survival)

summary(s14$Survival)
std.error(s14$Survival)

summary(s19$Survival)
std.error(s19$Survival)

summary(s24$Survival)
std.error(s24$Survival)


#surv_pre = surv[surv$Cross=="Conspecific",]

###Flowering
fl = surv[surv$Survival=="1",]

f9 =  (fl[fl$Cross_Region=="PANXCAR" & fl$Temperature==9,])
f14 = (fl[fl$Cross_Region=="PANXCAR" & fl$Temperature==14,])
f19 = (fl[fl$Cross_Region=="PANXCAR" & fl$Temperature==19,])
f24 = (fl[fl$Cross_Region=="PANXCAR" & fl$Temperature==24,])

summary(f9$Flowering)
std.error(f9$Flowering)

summary(f14$Flowering)
std.error(f14$Flowering)

summary(f19$Flowering)
std.error(f19$Flowering)

summary(f24$Flowering)
std.error(f24$Flowering)



#fl_pre = fl[fl$Cross=="Conspecific",]

###Flowering time
ft = surv[surv$Survival=="1",]
ft = fl[fl$Flowering=="1",]

ft9 =  (ft[ft$Cross_Region=="PANXCAR" & ft$Temperature==9,])
ft14 = (ft[ft$Cross_Region=="PANXCAR" & ft$Temperature==14,])
ft19 = (ft[ft$Cross_Region=="PANXCAR" & ft$Temperature==19,])
ft24 = (ft[ft$Cross_Region=="PANXCAR" & ft$Temperature==24,])

summary(ft9$Flowering_Time)
std.error(ft9$Flowering_Time)

summary(ft14$Flowering_Time)
std.error(ft14$Flowering_Time)

summary(ft19$Flowering_Time)
std.error(ft19$Flowering_Time)

summary(ft24$Flowering_Time)
std.error(ft24$Flowering_Time)

#ft_pre = ft[ft$Cross=="Conspecific",]

###Pollen viability
pv = fl[fl$Flowering=="1",]
pv$Pollen_Viab <- pv$Pollen_Viable/pv$Pollen_Total
pv$Pollen_Viability <- asin(sqrt(abs(pv$Pollen_Viable/pv$Pollen_Total)))


pv9 =  (pv[pv$Cross_Region=="PANXCAR" & pv$Temperature==9,])
pv14 = (pv[pv$Cross_Region=="PANXCAR" & pv$Temperature==14,])
pv19 = (pv[pv$Cross_Region=="PANXCAR" & pv$Temperature==19,])
pv24 = (pv[pv$Cross_Region=="PANXCAR" & pv$Temperature==24,])

summary(pv9$Pollen_Viab)
std.error(pv9$Pollen_Viab)

summary(pv14$Pollen_Viab)
std.error(pv14$Pollen_Viab)

summary(pv19$Pollen_Viab)
std.error(pv19$Pollen_Viab)

summary(pv24$Pollen_Viab)
std.error(pv24$Pollen_Viab)

#pv_pre = pv[pv$Cross=="Conspecific",]

###Seed set
ss = read.csv(file = "Seed_Set.csv", header = TRUE)

ss$Temperature <- as.numeric(ss$Temperature)
ss$Cross_Region <- as.factor(ss$Cross_Region)
ss$Mother_Population <- as.factor(ss$Mother_Population)
ss$Father_Population <- as.factor(ss$Father_Population)

ss9 =  (ss[ss$Cross_Region=="CPXCC" & ss$Temperature==9,])
ss14 = (ss[ss$Cross_Region=="CPXCC" & ss$Temperature==14,])
ss19 = (ss[ss$Cross_Region=="CPXCC" & ss$Temperature==19,])
ss24 = (ss[ss$Cross_Region=="CPXCC" & ss$Temperature==24,])

summary  (ss9$Seed_Set)
std.error(ss9$Seed_Set)

summary  (ss14$Seed_Set)
std.error(ss14$Seed_Set)

summary  (ss19$Seed_Set)
std.error(ss19$Seed_Set)

summary  (ss24$Seed_Set)
std.error(ss24$Seed_Set)

ss$Seed <- sqrt(ss$Seed_Set)

ss_pre <-ss[ss$Cross_Female=="Conspecific",]

ss_het<-ss[ss$Cross_Female=="Heterospecific",]
ss_con1<-ss[ss$Cross_Region=="CCXCC",]
ss_con2<-ss[ss$Cross_Region=="PPXPP",]
ss_pos<-rbind(ss_het,ss_con1,ss_con2)

#####MODELS#####
###Germination
germ_model = glmer(Germination ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Mother_ID) + (1|Father_ID), germ, family=binomial)
Anova(germ_model)
summ(germ_model)

###Survival
surv_model = glmer(Survival ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Mother_ID) + (1|Father_ID), surv, family=binomial)
Anova(surv_model)
summary(surv_model)
library(asbio)
press(surv_model)

r2beta(surv_model, partial = TRUE, method = "sgv", data = NULL)


###Flowering
fl_model = glmer(Flowering ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Mother_ID) + (1|Father_ID), fl, family=binomial)
Anova(fl_model)
summary(fl_model)

r2beta(fl_model, partial = TRUE, method = "sgv", data = NULL)

###Flowering time
ft_model_surv <- survreg(survival.object ~ Temperature*Cross_Region + Mother_Population + Father_Population + frailty(Mother_ID), data=ft, dist="logistic",x=TRUE,y=TRUE)

Anova(ft_model_surv)
summary(ft_model_surv)
r2beta(ft_model_surv, partial = TRUE, method = "sgv", data = NULL)

ft_model_lmer <-lmer(Flowering_Time~Temperature*Cross_Region + Mother_Population + Father_Population + (1|Mother_ID)+ (1|Father_ID), ft)
Anova(ft_model_lmer)
summary(ft_model_lmer)
r2beta(ft_model_lmer, partial = TRUE, method = "sgv", data = NULL)

###Pollen viability
pv_model   <- lmer(Pollen_Viability ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Mother_ID) + (1|Father_ID), pv)
Anova(pv_model)
summary(pv_model)
r2beta(pv_model, partial = TRUE, method = "sgv", data = NULL)

###Seed set
ss_pre_model<- lmer(Seed ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Female_ID) + (1|Male_ID), ss_pre)
Anova(ss_pre_model)
summary(ss_pre_model)
r2beta(ss_pre_model, partial = TRUE, method = "sgv", data = NULL)

ss_pos_model<- lmer(Seed ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Female_ID) + (1|Male_ID), ss_pos)
Anova(ss_pos_model)
summary(ss_pos_model)
r2beta(ss_pos_model, partial = TRUE, method = "sgv", data = NULL)

ss_model   <- lmer(Seed_Set ~ Temperature*Cross_Region + Mother_Population + Father_Population + (1|Female_ID) + (1|Male_ID), ss)
Anova(ss_model)
summary(ss_model)
r2beta(ss_model, partial = TRUE, method = "sgv", data = NULL)

######DIAGNOSTICS#####
###Germination
germ_pre_res <- simulateResiduals(fittedModel = germ_pre_model)
plotResiduals(germ_pre_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nGERMINATION", ylab="Standardized residual")
plotQQunif(germ_pre_res, testUniformity = F)
testDispersion(germ_pre_res)
testZeroInflation(germ_pre_res)

germ_res <- simulateResiduals(fittedModel = germ_model)
plotResiduals(germ_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nGERMINATION", ylab="Standardized residual")
plotQQunif(germ_res, testUniformity = F)
testDispersion(germ_res)
testZeroInflation(germ_res)

###Survival 
surv_pre_res <- simulateResiduals(fittedModel = surv_pre_model)
plotResiduals(surv_pre_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nSURVIVAL", ylab="Standardized residual")
plotQQunif(surv_pre_res, testUniformity = F)
testDispersion(surv_pre_res)
testZeroInflation(surv_pre_res)

surv_res <- simulateResiduals(fittedModel = surv_model)
plotResiduals(surv_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nSURVIVAL", ylab="Standardized residual")
plotQQunif(surv_res, testUniformity = F)
testDispersion(surv_res)
testZeroInflation(surv_res)


###Flowering
fl_pre_res <- simulateResiduals(fittedModel = fl_pre_model)
plotResiduals(fl_pre_res, rank = TRUE, quantreg = T,xlab="Predicted values (rank transformed)\nFLOWERING", ylab="Standardized residual")
plotQQunif(fl_pre_res, testUniformity = F)
testDispersion(fl_pre_res)
testZeroInflation(fl_pre_res)

fl_res <- simulateResiduals(fittedModel = fl_model)
plotResiduals(fl_res, rank = TRUE, quantreg = T,xlab="Predicted values (rank transformed)\nFLOWERING", ylab="Standardized residual")
plotQQunif(fl_res, testUniformity = F)
testDispersion(fl_res)
testZeroInflation(fl_res)

###Flowering time
ft_pre_res <- simulateResiduals(fittedModel = ft_pre_model_lmer)
plotResiduals(ft_pre_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nFLOWERING TIME", ylab="Standardized residual")
plotQQunif(ft_pre_res, testUniformity = F)
testDispersion(ft_pre_res)
testZeroInflation(ft_pre_res)

ft_res <- simulateResiduals(fittedModel = ft_model_lmer)
plotResiduals(ft_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nFLOWERING TIME", ylab="Standardized residual")
plotQQunif(ft_res, testUniformity = F)
testDispersion(ft_res)
testZeroInflation(ft_res)

###Pollen viability
pv_pre_res <- simulateResiduals(fittedModel = pv_pre_model)
plotResiduals(pv_pre_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nPOLLEN VIABILITY", ylab="Standardized residual")
plotQQunif(pv_pre_res, testUniformity = F)
testDispersion(pv_pre_res)
testZeroInflation(pv_pre_res)####Zero-Inflated

pv_res <- simulateResiduals(fittedModel = pv_model)
plotResiduals(pv_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nPOLLEN VIABILITY", ylab="Standardized residual")
plotQQunif(pv_res, testUniformity = F)
testDispersion(pv_res)
testZeroInflation(pv_res)####Zero-Inflated

###Seed set
ss_pre_res <- simulateResiduals(fittedModel = ss_pre_model)
plotResiduals(ss_pre_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nSEED SET", ylab="Standardized residual")
plotQQunif(ss_pre_res, testUniformity = F)
testDispersion(ss_pre_res)
testZeroInflation(ss_pre_res)####Zero-Inflated

ss_pos_res <- simulateResiduals(fittedModel = ss_pos_model)
plotResiduals(ss_pos_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nSEED SET", ylab="Standardized residual")
plotQQunif(ss_pos_res, testUniformity = F)
testDispersion(ss_pos_res)
testZeroInflation(ss_pos_res)

ss_res <- simulateResiduals(fittedModel = ss_model)
plotResiduals(ss_res, rank = TRUE, quantreg = T, xlab="Predicted values (rank transformed)\nSEED SET", ylab="Standardized residual")
plotQQunif(ss_res, testUniformity = F)
testDispersion(ss_res)
testZeroInflation(ss_res)####Zero-Inflated