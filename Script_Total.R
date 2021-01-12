library(CircStats)
library(circular)
library(ggplot2)
library(visreg)
library(nlme)
library(MASS)
library(MuMIn)
library(cowplot)
library(blme)
library(car)
library(lubridate)
library(glmmTMB)
library(dplyr)
library(DHARMa)


#### I°) Rayleigh tests on conceptions, births, and cycle resumptions

#we use the function r.test to compute the Rayleigh test vector length and pvalue associated
#we use the function circ.summary to compute the mean date (in radian, that we further converted in a date)

#A. Conceptions
TAB=read.csv2('TABLE3')
TAB=TAB[which(TAB$Uncertainty_Conc_Days<31),] #we only keep conceptions with less than 31 days of uncertainty
r.test(TAB$Conc_Radian) 
circ.summary(TAB$Conc_Radian) 

#B. Births
TAB=read.csv2('TABLE2')
TAB=TAB[which(TAB$Uncertainty.Days.<31),] #we only keep births with less than 31 days of uncertainty
TAB$DOB=dmy(TAB$DOB)
TAB=TAB[which(TAB$DOB>'2012-06-01'), ] #we only consider births occuring after the 1st of June 2012. 
r.test(TAB$DOB_Radian) 
circ.summary(TAB$DOB_Radian) 

#C. Cycle resumptions
TAB=read.csv2('TABLE5')

TAB=TAB[which(!TAB$TYPE=="Miscarriage Probable"),] #we did not consider the cycle resumption following a female miscarriage
r.test(TAB$CR_Radian) 
circ.summary(TAB$CR_Radian) 



#### II°) Proportions of time spent foraging (Model 1)

TAB=read.csv2('TABLE1')
TAB=TAB[which(TAB$temps_focal>60),]
TAB$id=as.factor(TAB$id)
TAB$YearFocal=as.factor(TAB$YearFocal)

#here, we first select the best phase, as the one minimizing the AIC of the following models;
MOD0=glmmTMB(Temps_Feed ~ sin(Date_Radian) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD1=glmmTMB(Temps_Feed ~ sin(Date_Radian + 1*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD2=glmmTMB(Temps_Feed ~ sin(Date_Radian + 2*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD3=glmmTMB(Temps_Feed ~ sin(Date_Radian + 3*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD4=glmmTMB(Temps_Feed ~ sin(Date_Radian + 4*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD5=glmmTMB(Temps_Feed ~ sin(Date_Radian + 5*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD6=glmmTMB(Temps_Feed ~ sin(Date_Radian + 6*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD7=glmmTMB(Temps_Feed ~ sin(Date_Radian + 7*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD8=glmmTMB(Temps_Feed ~ sin(Date_Radian + 8*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD9=glmmTMB(Temps_Feed ~ sin(Date_Radian + 9*pi/12) + 
               (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
             family='nbinom1',
             data=TAB)
MOD10=glmmTMB(Temps_Feed ~ sin(Date_Radian + 10*pi/12) + 
                (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
              family='nbinom1',
              data=TAB)
MOD11=glmmTMB(Temps_Feed ~ sin(Date_Radian + 11*pi/12) + 
                (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
              family='nbinom1',
              data=TAB)

AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)
# df      AIC
# MOD0   5 358551.2
# MOD1   5 358373.6
# MOD2   5 358248.1
# MOD3   5 358218.1
# MOD4   5 358289.7
# MOD5   5 358433.9
# MOD6   5 358605.1
# MOD7   5 358759.2
# MOD8   5 358862.5
# MOD9   5 358894.3
# MOD10  5 358846.2
# MOD11  5 358723.9

#the Model 3 is selected here. 
BEST_MOD=glmmTMB(Temps_Feed ~ sin(Date_Radian + 3*pi/12) + 
                   (1|id) + (1|YearFocal) + offset(log(temps_focal)), 
                 family='nbinom1',
                 data=TAB)
summary(BEST_MOD)
Anova(BEST_MOD)
confint(BEST_MOD)

#we ran the following lines for each ran model, in order to check the goddness of the fit of the models (these lines are not shown for following models)
vif(MOD)
r.squaredGLMM(BEST_MOD)
simulationOutput <- simulateResiduals(fittedModel = BEST_MOD, plot = F)
plot(simulationOutput)
testOverdispersion(simulationOutput)
testZeroInflation(simulationOutput)


#### III°) Infant mortality probability (Model 2)

TAB=read.csv2('TABLE2')
TAB=TAB[which(TAB$utiliser.pour.analyse.des.deces=='o'),]
TAB=TAB[which(TAB$Uncertainty.Days.<31),]
TAB$DOB=dmy(TAB$DOB)
TAB=TAB[which(TAB$DOB>'2012-06-01'), ]
TAB$Dead_Before_6months_old=as.factor(TAB$Dead_Before_6months_old)
TAB$Mother_id=as.factor(TAB$Mother_id)
TAB$Year_Birth_V2=as.factor(TAB$Year_Birth_V2)
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))
TAB$In_Peak_2Months_Yearly=as.factor(TAB$In_Peak_2Months_Yearly)
TAB$Birth_Lag_Yearly=as.numeric(scale(TAB$Birth_Lag_Yearly))

#MODEL 2A: 
MOD_Peak=bglmer(Dead_Before_6months_old ~ In_Peak_2Months_Yearly + 
                  Relative_Rank +  Parity + sexe +  
                  + (1|Mother_id) + (1|Year_Birth_V2), 
                glmerControl(optCtrl = list(maxfun = 20000)),
                family='binomial',data=TAB)
summary(MOD_Peak)
Anova(MOD_Peak)
confint(MOD_Peak, method='Wald')

#MODEL 2B:
MOD=bglmer(Dead_Before_6months_old ~ Birth_Lag_Yearly + 
             Relative_Rank +  Parity + sexe +  
             + (1|Mother_id) + (1|Year_Birth_V2), 
           glmerControl(optCtrl = list(maxfun = 20000)),
           family='binomial',data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')


#### Iv°) Female miscarriage probability (Model 3)

TAB=read.csv2('TABLE3')
TAB=TAB[which(is.na(TAB$Miscarriage)==FALSE),]
TAB=TAB[which(TAB$Uncertainty_Conc_Days<31),]
TAB$In_Peak_2months_Yearly=as.factor(TAB$In_Peak_2months_Yearly)
TAB$Miscarriage=as.factor(TAB$Miscarriage)
TAB$ID=as.factor(TAB$ID)
TAB$Year_Conc=as.factor(TAB$Year_Conc)
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))
TAB$Lag_Conc_Year=as.numeric(scale(TAB$Lag_Conc_Year))

#MODEL 3A:
MOD=bglmer(Miscarriage ~ In_Peak_2months_Yearly + 
             Relative_Rank + Parity +  
             (1|ID) + (1|Year_Conc), 
           glmerControl(optCtrl = list(maxfun = 20000)),
           family='binomial',data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')

#MODEL 3B:
MOD=glmer(Miscarriage ~ Lag_Conc_Year + 
            Relative_Rank + Parity +  
            (1|ID) + (1|Year_Conc), 
          glmerControl(optCtrl = list(maxfun = 20000)),
          family='binomial',data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')

#### v°) Female interbirth intervals (IBI) (Model 4)

TAB=read.csv2('TABLE4')
TAB=TAB[which(TAB$Dead_Before_6months_old1==0),]
TAB=TAB[which(TAB$Uncertainty.Days.1<31),]
TAB=TAB[which(TAB$Uncertainty.Days.2<31),]
TAB$DOB1=dmy(TAB$DOB1)
TAB=TAB[which(TAB$DOB1>"2012-06-01"), ]
summary(TAB$IBI)
hist(TAB$IBI)
sd(TAB$IBI) 
TAB$Mother_id=as.factor(TAB$Mother_id)
TAB$Year_Birth_V2=as.factor(TAB$Year_Birth_V2)
TAB$In_Peak_2Months_Yearly=as.factor(TAB$In_Peak_2Months_Yearly)
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))
TAB$Birth_Lag_Yearly=as.numeric(scale(TAB$Birth_Lag_Yearly))

#MODEL 4A:
MOD=blmer(IBI ~ In_Peak_2Months_Yearly + 
            Relative_Rank + Parity + sexe1 + 
            + (1|Mother_id) + (1|Year_Birth_V2), 
          data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')

#MODEL 4B: 
MOD=blmer(IBI ~ Birth_Lag_Yearly +
            Relative_Rank + Parity + sexe1 + 
            + (1|Mother_id) + (1|Year_Birth_V2), 
          data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')


#### vI°) Birth timing individual variation (Models 5 & 6)

TAB=read.csv2('TABLE2')
TAB=TAB[which(TAB$Uncertainty.Days.<31),]
TAB$DOB=dmy(TAB$DOB)
TAB=TAB[which(TAB$DOB>"2012-06-01"), ]
TAB=TAB[which(is.na(TAB$In_Peak_2Months_Yearly)==FALSE),]
TAB$Mother_id=as.factor(TAB$Mother_id)

TAB$Year_Birth_V2=as.factor(TAB$Year_Birth_V2)
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))
TAB$Age_Mother_Year=as.numeric(scale(TAB$Age_Mother_Year))

#MODEL 5:
MOD=bglmer(In_Peak_2Months_Yearly ~ Relative_Rank + Age_Mother_Year + Past_Repro + 
             (1|Mother_id) + (1|Year_Birth_V2),family='binomial', data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')

#to compute the effect of the random effect 'Mother_id':
MOD=bglmer(In_Peak_2Months_Yearly ~ Relative_Rank + Age_Mother_Year + Past_Repro + 
            (1|Mother_id) + (1|Year_Birth_V2),family='binomial', data=TAB)
MOD2=bglmer(In_Peak_2Months_Yearly ~ Relative_Rank + Age_Mother_Year + Past_Repro + 
             (1|Year_Birth_V2),family='binomial', data=TAB)
anova(MOD,MOD2)
  

#MODEL 6:
MOD=blmer(Birth_Lag_Yearly ~ Relative_Rank + Age_Mother_Year + Past_Repro + 
            (1|Mother_id) + (1|Year_Birth_V2),data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')

#to compute the effect of the random effect 'Mother_id':
MOD1=blmer(Birth_Lag_Yearly ~ Relative_Rank + Age_Mother_Year + Past_Repro + 
             (1|Mother_id) + (1|Year_Birth_V2),data=TAB)
MOD2=blmer(Birth_Lag_Yearly ~ Relative_Rank + Age_Mother_Year + Past_Repro + 
             (1|Year_Birth_V2),data=TAB)
anova(MOD1,MOD2)


