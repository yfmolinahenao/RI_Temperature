#####Germination#####
data <- read.csv(file="TEMPERATURE_EXPERIMENT_2.csv", stringsAsFactors = FALSE, na.strings = c("NA", ""))

germ_CC9  = data[data$Temperature=="9"  & data$Cross_Region=="CARXCAR", 16]
germ_CP9  = data[data$Temperature=="9"  & data$Cross_Region=="CARXPAN", 16]
germ_PC9  = data[data$Temperature=="9"  & data$Cross_Region=="PANXCAR", 16]
germ_PP9  = data[data$Temperature=="9"  & data$Cross_Region=="PANXPAN", 16]

germ_CC14 = data[data$Temperature=="14" & data$Cross_Region=="CARXCAR", 16]
germ_CP14 = data[data$Temperature=="14" & data$Cross_Region=="CARXPAN", 16]
germ_PC14 = data[data$Temperature=="14" & data$Cross_Region=="PANXCAR", 16]
germ_PP14 = data[data$Temperature=="14" & data$Cross_Region=="PANXPAN", 16]

germ_CC19 = data[data$Temperature=="19" & data$Cross_Region=="CARXCAR", 16]
germ_CP19 = data[data$Temperature=="19" & data$Cross_Region=="CARXPAN", 16]
germ_PC19 = data[data$Temperature=="19" & data$Cross_Region=="PANXCAR", 16]
germ_PP19 = data[data$Temperature=="19" & data$Cross_Region=="PANXPAN", 16]

germ_CC24 = data[data$Temperature=="24" & data$Cross_Region=="CARXCAR", 16]
germ_CP24 = data[data$Temperature=="24" & data$Cross_Region=="CARXPAN", 16]
germ_PC24 = data[data$Temperature=="24" & data$Cross_Region=="PANXCAR", 16]
germ_PP24 = data[data$Temperature=="24" & data$Cross_Region=="PANXPAN", 16]

##Bootstrap CIs##
set.seed(1)
Pre_Igerm_P2C9  = vector()
Pre_Igerm_C2P9  = vector()
Pre_Igerm_P2C14 = vector()
Pre_Igerm_C2P14 = vector()
Pre_Igerm_P2C19 = vector()
Pre_Igerm_C2P19 = vector()
Pre_Igerm_P2C24 = vector()
Pre_Igerm_C2P24 = vector()
Pos_Igerm_P2C9  = vector()
Pos_Igerm_C2P9  = vector()
Pos_Igerm_P2C14 = vector()
Pos_Igerm_C2P14 = vector()
Pos_Igerm_P2C19 = vector()
Pos_Igerm_C2P19 = vector()
Pos_Igerm_P2C24 = vector()
Pos_Igerm_C2P24 = vector()


for (i in 1:10000) {
  Rgerm_CC9 = sum(sample(germ_CC9 , 288, replace = T))
  Rgerm_CP9 = sum(sample(germ_CP9 , 288, replace = T))
  Rgerm_PC9 = sum(sample(germ_PC9 , 288, replace = T))
  Rgerm_PP9 = sum(sample(germ_PP9 , 288, replace = T))
  Rgerm_CC14= sum(sample(germ_CC14, 288, replace = T))
  Rgerm_CP14= sum(sample(germ_CP14, 288, replace = T))
  Rgerm_PC14= sum(sample(germ_PC14, 288, replace = T))
  Rgerm_PP14= sum(sample(germ_PP14, 288, replace = T))
  Rgerm_CC19= sum(sample(germ_CC19, 288, replace = T))
  Rgerm_CP19= sum(sample(germ_CP19, 288, replace = T))
  Rgerm_PC19= sum(sample(germ_PC19, 288, replace = T))
  Rgerm_PP19= sum(sample(germ_PP19, 288, replace = T))
  Rgerm_CC24= sum(sample(germ_CC24, 288, replace = T))
  Rgerm_CP24= sum(sample(germ_CP24, 288, replace = T))
  Rgerm_PC24= sum(sample(germ_PC24, 288, replace = T))
  Rgerm_PP24= sum(sample(germ_PP24, 288, replace = T))
  
  Pre_Igerm_P2C9  = c(Pre_Igerm_P2C9  , 1-2*Rgerm_PP9/(Rgerm_CC9+Rgerm_PP9))
  Pre_Igerm_C2P9  = c(Pre_Igerm_C2P9  , 1-2*Rgerm_CC9/(Rgerm_PP9+Rgerm_CC9))
  Pre_Igerm_P2C14 = c(Pre_Igerm_P2C14 , 1-2*Rgerm_PP14/(Rgerm_CC14+Rgerm_PP14))
  Pre_Igerm_C2P14 = c(Pre_Igerm_C2P14 , 1-2*Rgerm_CC14/(Rgerm_PP14+Rgerm_CC14))
  Pre_Igerm_P2C19 = c(Pre_Igerm_P2C19 , 1-2*Rgerm_PP19/(Rgerm_CC19+Rgerm_PP19))
  Pre_Igerm_C2P19 = c(Pre_Igerm_C2P19 , 1-2*Rgerm_CC19/(Rgerm_PP19+Rgerm_CC19))
  Pre_Igerm_P2C24 = c(Pre_Igerm_P2C24 , 1-2*Rgerm_PP24/(Rgerm_CC24+Rgerm_PP24))
  Pre_Igerm_C2P24 = c(Pre_Igerm_C2P24 , 1-2*Rgerm_CC24/(Rgerm_PP24+Rgerm_CC24))
  Pos_Igerm_P2C9  = c(Pos_Igerm_P2C9  , 1-2*(((Rgerm_CP9+Rgerm_PC9)/2)/  (Rgerm_CC9+((Rgerm_CP9+Rgerm_PC9)/2))))
  Pos_Igerm_C2P9  = c(Pos_Igerm_C2P9  , 1-2*(((Rgerm_CP9+Rgerm_PC9)/2)/  (Rgerm_PP9+((Rgerm_CP9+Rgerm_PC9)/2))))
  Pos_Igerm_P2C14 = c(Pos_Igerm_P2C14 , 1-2*(((Rgerm_CP14+Rgerm_PC14)/2)/(Rgerm_CC14+((Rgerm_CP14+Rgerm_PC14)/2))))
  Pos_Igerm_C2P14 = c(Pos_Igerm_C2P14 , 1-2*(((Rgerm_CP14+Rgerm_PC14)/2)/(Rgerm_PP14+((Rgerm_CP14+Rgerm_PC14)/2))))
  Pos_Igerm_P2C19 = c(Pos_Igerm_P2C19 , 1-2*(((Rgerm_CP19+Rgerm_PC19)/2)/(Rgerm_CC19+((Rgerm_CP19+Rgerm_PC19)/2))))
  Pos_Igerm_C2P19 = c(Pos_Igerm_C2P19 , 1-2*(((Rgerm_CP19+Rgerm_PC19)/2)/(Rgerm_PP19+((Rgerm_CP19+Rgerm_PC19)/2))))
  Pos_Igerm_P2C24 = c(Pos_Igerm_P2C24 , 1-2*(((Rgerm_CP24+Rgerm_PC24)/2)/(Rgerm_CC24+((Rgerm_CP24+Rgerm_PC24)/2))))
  Pos_Igerm_C2P24 = c(Pos_Igerm_C2P24 , 1-2*(((Rgerm_CP24+Rgerm_PC24)/2)/(Rgerm_PP24+((Rgerm_CP24+Rgerm_PC24)/2))))
  
}

#####Survival#####
data = data[data$Germination=="1",]

surv_CC9  = data[data$Temperature=="9"  & data$Cross_Region=="CARXCAR", 17]
surv_CP9  = data[data$Temperature=="9"  & data$Cross_Region=="CARXPAN", 17]
surv_PC9  = data[data$Temperature=="9"  & data$Cross_Region=="PANXCAR", 17]
surv_PP9  = data[data$Temperature=="9"  & data$Cross_Region=="PANXPAN", 17]

surv_CC14 = data[data$Temperature=="14" & data$Cross_Region=="CARXCAR", 17]
surv_CP14 = data[data$Temperature=="14" & data$Cross_Region=="CARXPAN", 17]
surv_PC14 = data[data$Temperature=="14" & data$Cross_Region=="PANXCAR", 17]
surv_PP14 = data[data$Temperature=="14" & data$Cross_Region=="PANXPAN", 17]

surv_CC19 = data[data$Temperature=="19" & data$Cross_Region=="CARXCAR", 17]
surv_CP19 = data[data$Temperature=="19" & data$Cross_Region=="CARXPAN", 17]
surv_PC19 = data[data$Temperature=="19" & data$Cross_Region=="PANXCAR", 17]
surv_PP19 = data[data$Temperature=="19" & data$Cross_Region=="PANXPAN", 17]

surv_CC24 = data[data$Temperature=="24" & data$Cross_Region=="CARXCAR", 17]
surv_CP24 = data[data$Temperature=="24" & data$Cross_Region=="CARXPAN", 17]
surv_PC24 = data[data$Temperature=="24" & data$Cross_Region=="PANXCAR", 17]
surv_PP24 = data[data$Temperature=="24" & data$Cross_Region=="PANXPAN", 17]

##Bootstrap CIs##
set.seed(1)
Pre_Isurv_P2C9  = vector()
Pre_Isurv_C2P9  = vector()
Pre_Isurv_P2C14 = vector()
Pre_Isurv_C2P14 = vector()
Pre_Isurv_P2C19 = vector()
Pre_Isurv_C2P19 = vector()
Pre_Isurv_P2C24 = vector()
Pre_Isurv_C2P24 = vector()
Pos_Isurv_P2C9  = vector()
Pos_Isurv_C2P9  = vector()
Pos_Isurv_P2C14 = vector()
Pos_Isurv_C2P14 = vector()
Pos_Isurv_P2C19 = vector()
Pos_Isurv_C2P19 = vector()
Pos_Isurv_P2C24 = vector()
Pos_Isurv_C2P24 = vector()


for (i in 1:10000) {
  Rsurv_CC9 = sum(sample(surv_CC9 , 100, replace = T))
  Rsurv_CP9 = sum(sample(surv_CP9 , 100, replace = T))
  Rsurv_PC9 = sum(sample(surv_PC9 , 100, replace = T))
  Rsurv_PP9 = sum(sample(surv_PP9 , 100, replace = T))
  Rsurv_CC14= sum(sample(surv_CC14, 100, replace = T))
  Rsurv_CP14= sum(sample(surv_CP14, 100, replace = T))
  Rsurv_PC14= sum(sample(surv_PC14, 100, replace = T))
  Rsurv_PP14= sum(sample(surv_PP14, 100, replace = T))
  Rsurv_CC19= sum(sample(surv_CC19, 100, replace = T))
  Rsurv_CP19= sum(sample(surv_CP19, 100, replace = T))
  Rsurv_PC19= sum(sample(surv_PC19, 100, replace = T))
  Rsurv_PP19= sum(sample(surv_PP19, 100, replace = T))
  Rsurv_CC24= sum(sample(surv_CC24, 100, replace = T))
  Rsurv_CP24= sum(sample(surv_CP24, 100, replace = T))
  Rsurv_PC24= sum(sample(surv_PC24, 100, replace = T))
  Rsurv_PP24= sum(sample(surv_PP24, 100, replace = T))
  
  Pre_Isurv_P2C9  = c(Pre_Isurv_P2C9  , 1-2*Rsurv_PP9/(Rsurv_CC9+Rsurv_PP9))
  Pre_Isurv_C2P9  = c(Pre_Isurv_C2P9  , 1-2*Rsurv_CC9/(Rsurv_PP9+Rsurv_CC9))
  Pre_Isurv_P2C14 = c(Pre_Isurv_P2C14 , 1-2*Rsurv_PP14/(Rsurv_CC14+Rsurv_PP14))
  Pre_Isurv_C2P14 = c(Pre_Isurv_C2P14 , 1-2*Rsurv_CC14/(Rsurv_PP14+Rsurv_CC14))
  Pre_Isurv_P2C19 = c(Pre_Isurv_P2C19 , 1-2*Rsurv_PP19/(Rsurv_CC19+Rsurv_PP19))
  Pre_Isurv_C2P19 = c(Pre_Isurv_C2P19 , 1-2*Rsurv_CC19/(Rsurv_PP19+Rsurv_CC19))
  Pre_Isurv_P2C24 = c(Pre_Isurv_P2C24 , 1-2*Rsurv_PP24/(Rsurv_CC24+Rsurv_PP24))
  Pre_Isurv_C2P24 = c(Pre_Isurv_C2P24 , 1-2*Rsurv_CC24/(Rsurv_PP24+Rsurv_CC24))
  Pos_Isurv_P2C9  = c(Pos_Isurv_P2C9  , 1-2*(((Rsurv_CP9+Rsurv_PC9)/2)/  (Rsurv_CC9+((Rsurv_CP9+Rsurv_PC9)/2))))
  Pos_Isurv_C2P9  = c(Pos_Isurv_C2P9  , 1-2*(((Rsurv_CP9+Rsurv_PC9)/2)/  (Rsurv_PP9+((Rsurv_CP9+Rsurv_PC9)/2))))
  Pos_Isurv_P2C14 = c(Pos_Isurv_P2C14 , 1-2*(((Rsurv_CP14+Rsurv_PC14)/2)/(Rsurv_CC14+((Rsurv_CP14+Rsurv_PC14)/2))))
  Pos_Isurv_C2P14 = c(Pos_Isurv_C2P14 , 1-2*(((Rsurv_CP14+Rsurv_PC14)/2)/(Rsurv_PP14+((Rsurv_CP14+Rsurv_PC14)/2))))
  Pos_Isurv_P2C19 = c(Pos_Isurv_P2C19 , 1-2*(((Rsurv_CP19+Rsurv_PC19)/2)/(Rsurv_CC19+((Rsurv_CP19+Rsurv_PC19)/2))))
  Pos_Isurv_C2P19 = c(Pos_Isurv_C2P19 , 1-2*(((Rsurv_CP19+Rsurv_PC19)/2)/(Rsurv_PP19+((Rsurv_CP19+Rsurv_PC19)/2))))
  Pos_Isurv_P2C24 = c(Pos_Isurv_P2C24 , 1-2*(((Rsurv_CP24+Rsurv_PC24)/2)/(Rsurv_CC24+((Rsurv_CP24+Rsurv_PC24)/2))))
  Pos_Isurv_C2P24 = c(Pos_Isurv_C2P24 , 1-2*(((Rsurv_CP24+Rsurv_PC24)/2)/(Rsurv_PP24+((Rsurv_CP24+Rsurv_PC24)/2))))
  
}

#####Flowering Time#####
data <- data[data$Germination=="1",]
data <- data[data$Survival=="1",]
data <- data[data$Flowering=="1",]


FT_CC9   = data[data$Temperature=="9"  & data$Cross_Region=="CARXCAR", 20]
FT_CP9   = data[data$Temperature=="9"  & data$Cross_Region=="CARXPAN", 20]
FT_PC9   = data[data$Temperature=="9"  & data$Cross_Region=="PANXCAR", 20]
FT_PP9   = data[data$Temperature=="9"  & data$Cross_Region=="PANXPAN", 20]
FT_H9    = data[data$Temperature=="9"  & data$Cross=="Heterospecific", 20]
FT_CON9  = data[data$Temperature=="9"  & data$Cross=="Conspecific"   , 20]

FT_CC14  = data[data$Temperature=="14" & data$Cross_Region=="CARXCAR", 20]
FT_CP14  = data[data$Temperature=="14" & data$Cross_Region=="CARXPAN", 20]
FT_PC14  = data[data$Temperature=="14" & data$Cross_Region=="PANXCAR", 20]
FT_PP14  = data[data$Temperature=="14" & data$Cross_Region=="PANXPAN", 20]
FT_H14   = data[data$Temperature=="14" & data$Cross=="Heterospecific", 20]
FT_CON14 = data[data$Temperature=="14" & data$Cross=="Conspecific"   , 20]

FT_CC19  = data[data$Temperature=="19" & data$Cross_Region=="CARXCAR", 20]
FT_CP19  = data[data$Temperature=="19" & data$Cross_Region=="CARXPAN", 20]
FT_PC19  = data[data$Temperature=="19" & data$Cross_Region=="PANXCAR", 20]
FT_PP19  = data[data$Temperature=="19" & data$Cross_Region=="PANXPAN", 20]
FT_H19   = data[data$Temperature=="19" & data$Cross=="Heterospecific", 20]
FT_CON19 = data[data$Temperature=="19" & data$Cross=="Conspecific"   , 20]

FT_CC24  = data[data$Temperature=="24" & data$Cross_Region=="CARXCAR", 20]
FT_CP24  = data[data$Temperature=="24" & data$Cross_Region=="CARXPAN", 20]
FT_PC24  = data[data$Temperature=="24" & data$Cross_Region=="PANXCAR", 20]
FT_PP24  = data[data$Temperature=="24" & data$Cross_Region=="PANXPAN", 20]
FT_H24   = data[data$Temperature=="24" & data$Cross=="Heterospecific", 20]
FT_CON24 = data[data$Temperature=="24" & data$Cross=="Conspecific"   , 20]

library(overlapping)

##Bootstrap CIs##
set.seed(1)
Pre_IFT_9    = vector()
Pre_IFT_14   = vector()
Pre_IFT_19   = vector()
Pre_IFT_24   = vector()
Pos_IFT_P2C9 = vector()
Pos_IFT_C2P9 = vector()
Pos_IFT_P2C14= vector()
Pos_IFT_C2P14= vector()
Pos_IFT_P2C19= vector()
Pos_IFT_C2P19= vector()
Pos_IFT_P2C24= vector()
Pos_IFT_C2P24= vector()

for (i in 1:1000) {
  
  
  RFT_CC9 =sample(FT_CC9, replace = T)
  RFT_PP9 =sample(FT_PP9, replace = T)
  RFT9 = list(RFT_CC9, RFT_PP9)
  out_9 <- overlap(RFT9, plot=F)
  S_9 = out_9$OV
  
  RFT_CC14=sample(FT_CC14, replace = T)
  RFT_PP14=sample(FT_PP14, replace = T)
  RFT14 = list(RFT_CC14, RFT_PP14)
  out_14 <- overlap(RFT14,plot=F)
  S_14 = out_14$OV
  
  RFT_CC19=sample(FT_CC19,replace = T)
  RFT_PP19=sample(FT_PP19, replace = T)
  RFT19 = list(RFT_CC19, RFT_PP19)
  out_19 <- overlap(RFT19,plot=F)
  S_19 = out_19$OV
  
  
  RFT_CC24=sample(FT_CC24,replace = T)
  RFT_PP24=sample(FT_PP24,replace = T)
  RFT24 = list(RFT_CC24, RFT_PP24)
  out_24 <- overlap(RFT24,plot=F)
  S_24 = out_24$OV
  
  RFT_CC9 =sample(FT_CC9, replace = T)
  RFT_H9 =sample(FT_H9, replace = T)
  RFTC9 = list(RFT_CC9, RFT_H9)
  out_C9 <- overlap(RFTC9, plot=F)
  S_C9 = out_C9$OV
  
  RFT_H9 =sample(FT_H9, replace = T)
  RFT_PP9 =sample(FT_PP9, replace = T)
  RFTP9 = list(RFT_H9, RFT_PP9)
  out_P9 <- overlap(RFTP9, plot=F)
  S_P9 = out_P9$OV
  
  RFT_CC14 =sample(FT_CC14, replace = T)
  RFT_H14 =sample(FT_H14, replace = T)
  RFTC14 = list(RFT_CC14, RFT_H14)
  out_C14 <- overlap(RFTC14, plot=F)
  S_C14 = out_C14$OV
  
  RFT_H14 =sample(FT_H14, replace = T)
  RFT_PP14 =sample(FT_PP14, replace = T)
  RFTP14 = list(RFT_H14, RFT_PP14)
  out_P14 <- overlap(RFTP14, plot=F)
  S_P14 = out_P14$OV
  
  RFT_CC19 =sample(FT_CC19, replace = T)
  RFT_H19 =sample(FT_H19, replace = T)
  RFTC19 = list(RFT_CC19, RFT_H19)
  out_C19 <- overlap(RFTC19, plot=F)
  S_C19 = out_C19$OV
  
  RFT_H19 =sample(FT_H19, replace = T)
  RFT_PP19 =sample(FT_PP19, replace = T)
  RFTP19 = list(RFT_H19, RFT_PP19)
  out_P19 <- overlap(RFTP19, plot=F)
  S_P19 = out_P19$OV
  
  RFT_CC24 =sample(FT_CC24, replace = T)
  RFT_H24 =sample(FT_H24, replace = T)
  RFTC24 = list(RFT_CC24, RFT_H24)
  out_C24 <- overlap(RFTC24, plot=F)
  S_C24 = out_C24$OV
  
  RFT_H24 =sample(FT_H24, replace = T)
  RFT_PP24 =sample(FT_PP24, replace = T)
  RFTP24 = list(RFT_H24, RFT_PP24)
  out_P24 <- overlap(RFTP24, plot=F)
  S_P24 = out_P24$OV
  
  Pre_IFT_9     = c(Pre_IFT_9     , (1 -  S_9))  
  Pre_IFT_14    = c(Pre_IFT_14    , (1 - S_14))
  Pre_IFT_19    = c(Pre_IFT_19    , (1 - S_19))
  Pre_IFT_24    = c(Pre_IFT_24    , (1 - S_24))
  Pos_IFT_P2C9  = c(Pos_IFT_P2C9  , (1 -  S_C9)) 
  Pos_IFT_C2P9  = c(Pos_IFT_C2P9  , (1 -  S_P9))
  Pos_IFT_P2C14 = c(Pos_IFT_P2C14 , (1 - S_C14))
  Pos_IFT_C2P14 = c(Pos_IFT_C2P14 , (1 - S_P14))
  Pos_IFT_P2C19 = c(Pos_IFT_P2C19 , (1 - S_C19))
  Pos_IFT_C2P19 = c(Pos_IFT_C2P19 , (1 - S_P19))
  Pos_IFT_P2C24 = c(Pos_IFT_P2C24 , (1 - S_C24))
  Pos_IFT_C2P24 = c(Pos_IFT_C2P24 , (1 - S_P24))
  
}


#####Pollen Viability#####

data$pv = data$Pollen_Viable/data$Pollen_Total

pv_CC9  = data[data$Temperature=="9"  & data$Cross_Region=="CARXCAR", 24]
pv_CP9  = data[data$Temperature=="9"  & data$Cross_Region=="CARXPAN", 24]
pv_PC9  = data[data$Temperature=="9"  & data$Cross_Region=="PANXCAR", 24]
pv_PP9  = data[data$Temperature=="9"  & data$Cross_Region=="PANXPAN", 24]

pv_CC14 = data[data$Temperature=="14" & data$Cross_Region=="CARXCAR", 24]
pv_CP14 = data[data$Temperature=="14" & data$Cross_Region=="CARXPAN", 24]
pv_PC14 = data[data$Temperature=="14" & data$Cross_Region=="PANXCAR", 24]
pv_PP14 = data[data$Temperature=="14" & data$Cross_Region=="PANXPAN", 24]

pv_CC19 = data[data$Temperature=="19" & data$Cross_Region=="CARXCAR", 24]
pv_CP19 = data[data$Temperature=="19" & data$Cross_Region=="CARXPAN", 24]
pv_PC19 = data[data$Temperature=="19" & data$Cross_Region=="PANXCAR", 24]
pv_PP19 = data[data$Temperature=="19" & data$Cross_Region=="PANXPAN", 24]

pv_CC24 = data[data$Temperature=="24" & data$Cross_Region=="CARXCAR", 24]
pv_CP24 = data[data$Temperature=="24" & data$Cross_Region=="CARXPAN", 24]
pv_PC24 = data[data$Temperature=="24" & data$Cross_Region=="PANXCAR", 24]
pv_PP24 = data[data$Temperature=="24" & data$Cross_Region=="PANXPAN", 24]


##Bootstrap CIs##
set.seed(1)
Pre_Ipv_P2C9  = vector()
Pre_Ipv_C2P9  = vector()
Pre_Ipv_P2C14 = vector()
Pre_Ipv_C2P14 = vector()
Pre_Ipv_P2C19 = vector()
Pre_Ipv_C2P19 = vector()
Pre_Ipv_P2C24 = vector()
Pre_Ipv_C2P24 = vector()
Pos_Ipv_P2C9  = vector()
Pos_Ipv_C2P9  = vector()
Pos_Ipv_P2C14 = vector()
Pos_Ipv_C2P14 = vector()
Pos_Ipv_P2C19 = vector()
Pos_Ipv_C2P19 = vector()
Pos_Ipv_P2C24 = vector()
Pos_Ipv_C2P24 = vector()

for (i in 1:10000) {
  Rpv_CC9 = sum(sample(pv_CC9 , 100, replace = T))
  Rpv_CP9 = sum(sample(pv_CP9 , 100, replace = T))
  Rpv_PC9 = sum(sample(pv_PC9 , 100, replace = T))
  Rpv_PP9 = sum(sample(pv_PP9 , 100, replace = T))
  Rpv_CC14= sum(sample(pv_CC14, 100, replace = T))
  Rpv_CP14= sum(sample(pv_CP14, 100, replace = T))
  Rpv_PC14= sum(sample(pv_PC14, 100, replace = T))
  Rpv_PP14= sum(sample(pv_PP14, 100, replace = T))
  Rpv_CC19= sum(sample(pv_CC19, 100, replace = T))
  Rpv_CP19= sum(sample(pv_CP19, 100, replace = T))
  Rpv_PC19= sum(sample(pv_PC19, 100, replace = T))
  Rpv_PP19= sum(sample(pv_PP19, 100, replace = T))
  Rpv_CC24= sum(sample(pv_CC24, 100, replace = T))
  Rpv_CP24= sum(sample(pv_CP24, 100, replace = T))
  Rpv_PC24= sum(sample(pv_PC24, 100, replace = T))
  Rpv_PP24= sum(sample(pv_PP24, 100, replace = T))
  
  Pre_Ipv_P2C9  = c(Pre_Ipv_P2C9  , 1-2*Rpv_PP9/(Rpv_CC9+Rpv_PP9))
  Pre_Ipv_C2P9  = c(Pre_Ipv_C2P9  , 1-2*Rpv_CC9/(Rpv_PP9+Rpv_CC9))
  Pre_Ipv_P2C14 = c(Pre_Ipv_P2C14 , 1-2*Rpv_PP14/(Rpv_CC14+Rpv_PP14))
  Pre_Ipv_C2P14 = c(Pre_Ipv_C2P14 , 1-2*Rpv_CC14/(Rpv_PP14+Rpv_CC14))
  Pre_Ipv_P2C19 = c(Pre_Ipv_P2C19 , 1-2*Rpv_PP19/(Rpv_CC19+Rpv_PP19))
  Pre_Ipv_C2P19 = c(Pre_Ipv_C2P19 , 1-2*Rpv_CC19/(Rpv_PP19+Rpv_CC19))
  Pre_Ipv_P2C24 = c(Pre_Ipv_P2C24 , 1-2*Rpv_PP24/(Rpv_CC24+Rpv_PP24))
  Pre_Ipv_C2P24 = c(Pre_Ipv_C2P24 , 1-2*Rpv_CC24/(Rpv_PP24+Rpv_CC24))
  Pos_Ipv_P2C9  = c(Pos_Ipv_P2C9  , 1-2*(((Rpv_CP9+Rpv_PC9)/2)/  (Rpv_CC9+((Rpv_CP9+Rpv_PC9)/2))))
  Pos_Ipv_C2P9  = c(Pos_Ipv_C2P9  , 1-2*(((Rpv_CP9+Rpv_PC9)/2)/  (Rpv_PP9+((Rpv_CP9+Rpv_PC9)/2))))
  Pos_Ipv_P2C14 = c(Pos_Ipv_P2C14 , 1-2*(((Rpv_CP14+Rpv_PC14)/2)/(Rpv_CC14+((Rpv_CP14+Rpv_PC14)/2))))
  Pos_Ipv_C2P14 = c(Pos_Ipv_C2P14 , 1-2*(((Rpv_CP14+Rpv_PC14)/2)/(Rpv_PP14+((Rpv_CP14+Rpv_PC14)/2))))
  Pos_Ipv_P2C19 = c(Pos_Ipv_P2C19 , 1-2*(((Rpv_CP19+Rpv_PC19)/2)/(Rpv_CC19+((Rpv_CP19+Rpv_PC19)/2))))
  Pos_Ipv_C2P19 = c(Pos_Ipv_C2P19 , 1-2*(((Rpv_CP19+Rpv_PC19)/2)/(Rpv_PP19+((Rpv_CP19+Rpv_PC19)/2))))
  Pos_Ipv_P2C24 = c(Pos_Ipv_P2C24 , 1-2*(((Rpv_CP24+Rpv_PC24)/2)/(Rpv_CC24+((Rpv_CP24+Rpv_PC24)/2))))
  Pos_Ipv_C2P24 = c(Pos_Ipv_C2P24 , 1-2*(((Rpv_CP24+Rpv_PC24)/2)/(Rpv_PP24+((Rpv_CP24+Rpv_PC24)/2))))

}
#####Seed Set#####

data_ss <- read.csv(file="Seed_Set.csv", stringsAsFactors = FALSE, na.strings = c("NA", ""))

ss_CC9 = data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="CCXCC", 13]
ss_CP9 = data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="CCXPP", 13]
ss_PC9 = data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="PPXCC", 13]
ss_PP9 = data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="PPXPP", 13]
ss_CPC9= data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="CPXCC", 13]
ss_CPP9= data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="CPXPP", 13]
ss_PCC9= data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="PCXCC", 13]
ss_PCP9= data_ss[data_ss$Temperature=="9" & data_ss$Cross_Region=="PCXPP", 13]

ss_CC14 = data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="CCXCC", 13]
ss_CP14 = data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="CCXPP", 13]
ss_PC14 = data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="PPXCC", 13]
ss_PP14 = data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="PPXPP", 13]
ss_CPC14= data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="CPXCC", 13]
ss_CPP14= data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="CPXPP", 13]
ss_PCC14= data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="PCXCC", 13]
ss_PCP14= data_ss[data_ss$Temperature=="14" & data_ss$Cross_Region=="PCXPP", 13]

ss_CC19 = data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="CCXCC", 13]
ss_CP19 = data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="CCXPP", 13]
ss_PC19 = data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="PPXCC", 13]
ss_PP19 = data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="PPXPP", 13]
ss_CPC19= data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="CPXCC", 13]
ss_CPP19= data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="CPXPP", 13]
ss_PCC19= data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="PCXCC", 13]
ss_PCP19= data_ss[data_ss$Temperature=="19" & data_ss$Cross_Region=="PCXPP", 13]

ss_CC24 = data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="CCXCC", 13]
ss_CP24 = data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="CCXPP", 13]
ss_PC24 = data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="PPXCC", 13]
ss_PP24 = data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="PPXPP", 13]
ss_CPC24= data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="CPXCC", 13]
ss_CPP24= data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="CPXPP", 13]
ss_PCC24= data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="PCXCC", 13]
ss_PCP24= data_ss[data_ss$Temperature=="24" & data_ss$Cross_Region=="PCXPP", 13]


##Bootstrap CIs##
set.seed(1)
Pre_Iss_P2C9  = vector()
Pre_Iss_C2P9  = vector()
Pre_Iss_P2C14 = vector()
Pre_Iss_C2P14 = vector()
Pre_Iss_P2C19 = vector()
Pre_Iss_C2P19 = vector()
Pre_Iss_P2C24 = vector()
Pre_Iss_C2P24 = vector()
Pos_Iss_P2C9  = vector()
Pos_Iss_C2P9  = vector()
Pos_Iss_P2C14 = vector()
Pos_Iss_C2P14 = vector()
Pos_Iss_P2C19 = vector()
Pos_Iss_C2P19 = vector()
Pos_Iss_P2C24 = vector()
Pos_Iss_C2P24 = vector()


for (i in 1:10000) {
  Rss_CC9= sum(sample(ss_CC9,   24, replace = T))
  Rss_PP9= sum(sample(ss_PP9,   24, replace = T))
  Rss_CP9= sum(sample(ss_CP9,   24, replace = T))
  Rss_PC9= sum(sample(ss_PC9,   24, replace = T))
  
  Rss_CC14= sum(sample(ss_CC14, 24, replace = T))
  Rss_PP14= sum(sample(ss_PP14, 24, replace = T))
  Rss_CP14= sum(sample(ss_CP14, 24, replace = T))
  Rss_PC14= sum(sample(ss_PC14, 24, replace = T))
  
  Rss_CC19= sum(sample(ss_CC19, 24, replace = T))
  Rss_PP19= sum(sample(ss_PP19, 24, replace = T))
  Rss_CP19= sum(sample(ss_CP19, 24, replace = T))
  Rss_PC19= sum(sample(ss_PC19, 24, replace = T))
  
  Rss_CC24= sum(sample(ss_CC24, 24, replace = T))
  Rss_PP24= sum(sample(ss_PP24, 24, replace = T))
  Rss_CP24= sum(sample(ss_CP24, 24, replace = T))
  Rss_PC24= sum(sample(ss_PC24, 24, replace = T))
  
  
  Pre_Iss_P2C9  =  c(Pre_Iss_P2C9 , 1-2*Rss_CP9/(Rss_CP9+Rss_CC9))
  Pre_Iss_C2P9  =  c(Pre_Iss_C2P9 , 1-2*Rss_PC9/(Rss_PC9+Rss_PP9))
  Pre_Iss_P2C14 =  c(Pre_Iss_P2C14, 1-2*Rss_CP14/(Rss_CP14+Rss_CC14))
  Pre_Iss_C2P14 =  c(Pre_Iss_C2P14, 1-2*Rss_PC14/(Rss_PC14+Rss_PP14))
  Pre_Iss_P2C19 =  c(Pre_Iss_P2C19, 1-2*Rss_CP19/(Rss_CP19+Rss_CC19))
  Pre_Iss_C2P19 =  c(Pre_Iss_C2P19, 1-2*Rss_PC19/(Rss_PC19+Rss_PP19))
  Pre_Iss_P2C24 =  c(Pre_Iss_P2C24, 1-2*Rss_CP24/(Rss_CP24+Rss_CC24))
  Pre_Iss_C2P24 =  c(Pre_Iss_C2P24, 1-2*Rss_PC24/(Rss_PC24+Rss_PP24))
  Pos_Iss_P2C9  =  c(Pos_Iss_P2C9 ,  1-2*((Rss_CP9+Rss_PC9)/2)/(((Rss_CP9+Rss_PC9)/2)+Rss_CC9))
  Pos_Iss_C2P9  =  c(Pos_Iss_C2P9 ,  1-2*((Rss_CP9+Rss_PC9)/2)/(((Rss_CP9+Rss_PC9)/2)+Rss_PP9))
  Pos_Iss_P2C14 =  c(Pos_Iss_P2C14,  1-2*((Rss_CP14+Rss_PC14)/2)/(((Rss_CP14+Rss_PC14)/2)+Rss_CC14))
  Pos_Iss_C2P14 =  c(Pos_Iss_C2P14,  1-2*((Rss_CP14+Rss_PC14)/2)/(((Rss_CP14+Rss_PC14)/2)+Rss_PP14))
  Pos_Iss_P2C19 =  c(Pos_Iss_P2C19,  1-2*((Rss_CP19+Rss_PC19)/2)/(((Rss_CP19+Rss_PC19)/2)+Rss_CC19))
  Pos_Iss_C2P19 =  c(Pos_Iss_C2P19,  1-2*((Rss_CP19+Rss_PC19)/2)/(((Rss_CP19+Rss_PC19)/2)+Rss_PP19))
  Pos_Iss_P2C24 =  c(Pos_Iss_P2C24,  1-2*((Rss_CP24+Rss_PC24)/2)/(((Rss_CP24+Rss_PC24)/2)+Rss_CC24))
  Pos_Iss_C2P24 =  c(Pos_Iss_C2P24,  1-2*((Rss_CP24+Rss_PC24)/2)/(((Rss_CP24+Rss_PC24)/2)+Rss_PP24))
}

#####TOTAL#####

##Bootstrap CIs##
set.seed(1)

APre_Isurv_P2C9  = vector()
APre_Isurv_C2P9  = vector()
APre_Isurv_P2C14 = vector()
APre_Isurv_C2P14 = vector()
APre_Isurv_P2C19 = vector()
APre_Isurv_C2P19 = vector()
APre_Isurv_P2C24 = vector()
APre_Isurv_C2P24 = vector()

APre_IFT_P2C9 = vector()
APre_IFT_C2P9 = vector()
APre_IFT_P2C14= vector()
APre_IFT_C2P14= vector()
APre_IFT_P2C19= vector()
APre_IFT_C2P19= vector()
APre_IFT_P2C24= vector()
APre_IFT_C2P24= vector()


APre_Ipv_P2C9  = vector()
APre_Ipv_C2P9  = vector()
APre_Ipv_P2C14 = vector()
APre_Ipv_C2P14 = vector()
APre_Ipv_P2C19 = vector()
APre_Ipv_C2P19 = vector()
APre_Ipv_P2C24 = vector()
APre_Ipv_C2P24 = vector()
 
APre_Iss_P2C9  = vector()
APre_Iss_C2P9  = vector()
APre_Iss_P2C14 = vector()
APre_Iss_C2P14 = vector()
APre_Iss_P2C19 = vector()
APre_Iss_C2P19 = vector()
APre_Iss_P2C24 = vector()
APre_Iss_C2P24 = vector()

APos_Igerm_P2C9  = vector()
APos_Igerm_C2P9  = vector()
APos_Igerm_P2C14 = vector()
APos_Igerm_C2P14 = vector()
APos_Igerm_P2C19 = vector()
APos_Igerm_C2P19 = vector()
APos_Igerm_P2C24 = vector()
APos_Igerm_C2P24 = vector()

APos_Isurv_P2C9  = vector()
APos_Isurv_C2P9  = vector()
APos_Isurv_P2C14 = vector()
APos_Isurv_C2P14 = vector()
APos_Isurv_P2C19 = vector()
APos_Isurv_C2P19 = vector()
APos_Isurv_P2C24 = vector()
APos_Isurv_C2P24 = vector()

APos_IFT_P2C9  = vector()
APos_IFT_C2P9  = vector()
APos_IFT_P2C14 = vector()
APos_IFT_C2P14 = vector()
APos_IFT_P2C19 = vector()
APos_IFT_C2P19 = vector()
APos_IFT_P2C24 = vector()
APos_IFT_C2P24 = vector()

APos_Ipv_P2C9  = vector()
APos_Ipv_C2P9  = vector()
APos_Ipv_P2C14 = vector()
APos_Ipv_C2P14 = vector()
APos_Ipv_P2C19 = vector()
APos_Ipv_C2P19 = vector()
APos_Ipv_P2C24 = vector()
APos_Ipv_C2P24 = vector()

APos_Iss_P2C9  = vector()
APos_Iss_C2P9  = vector()
APos_Iss_P2C14 = vector()
APos_Iss_C2P14 = vector()
APos_Iss_P2C19 = vector()
APos_Iss_C2P19 = vector()
APos_Iss_P2C24 = vector()
APos_Iss_C2P24 = vector()

Itotal_P2C9  = vector()
Itotal_C2P9  = vector()
Itotal_P2C14 = vector()
Itotal_C2P14 = vector()
Itotal_P2C19 = vector()
Itotal_C2P19 = vector()
Itotal_P2C24 = vector()
Itotal_C2P24 = vector()


for (i in 1:10000) {
  
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  RPre_Isurv_P2C9  = (sample(Pre_Isurv_P2C9  , 1, replace = T))
  RPre_Isurv_C2P9  = (sample(Pre_Isurv_C2P9  , 1, replace = T))
  RPre_Isurv_P2C14 = (sample(Pre_Isurv_P2C14 , 1, replace = T))
  RPre_Isurv_C2P14 = (sample(Pre_Isurv_C2P14 , 1, replace = T))
  RPre_Isurv_P2C19 = (sample(Pre_Isurv_P2C19 , 1, replace = T))
  RPre_Isurv_C2P19 = (sample(Pre_Isurv_C2P19 , 1, replace = T))
  RPre_Isurv_P2C24 = (sample(Pre_Isurv_P2C24 , 1, replace = T))
  RPre_Isurv_C2P24 = (sample(Pre_Isurv_C2P24 , 1, replace = T))
  
  APre_Isurv_P2C9  = c(APre_Isurv_P2C9 , (1-RPre_Igerm_P2C9 )*RPre_Isurv_P2C9 )
  APre_Isurv_C2P9  = c(APre_Isurv_C2P9 , (1-RPre_Igerm_C2P9 )*RPre_Isurv_C2P9 )
  APre_Isurv_P2C14 = c(APre_Isurv_P2C14, (1-RPre_Igerm_P2C14)*RPre_Isurv_P2C14)
  APre_Isurv_C2P14 = c(APre_Isurv_C2P14, (1-RPre_Igerm_C2P14)*RPre_Isurv_C2P14)
  APre_Isurv_P2C19 = c(APre_Isurv_P2C19, (1-RPre_Igerm_P2C19)*RPre_Isurv_P2C19)
  APre_Isurv_C2P19 = c(APre_Isurv_C2P19, (1-RPre_Igerm_C2P19)*RPre_Isurv_C2P19)
  APre_Isurv_P2C24 = c(APre_Isurv_P2C24, (1-RPre_Igerm_P2C24)*RPre_Isurv_P2C24)
  APre_Isurv_C2P24 = c(APre_Isurv_C2P24, (1-RPre_Igerm_C2P24)*RPre_Isurv_C2P24)
  
}
  
for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  RPre_IFT_9   = (sample(Pre_IFT_9  , 1, replace=T))
  RPre_IFT_14  = (sample(Pre_IFT_14 , 1, replace=T))
  RPre_IFT_19  = (sample(Pre_IFT_19 , 1, replace=T))
  RPre_IFT_24  = (sample(Pre_IFT_24 , 1, replace=T))
  
  APre_IFT_P2C9 =c(APre_IFT_P2C9 , (1-(aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPre_IFT_9)
  APre_IFT_C2P9 =c(APre_IFT_C2P9 , (1-(aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPre_IFT_9)
  APre_IFT_P2C14=c(APre_IFT_P2C14, (1-(aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPre_IFT_14)
  APre_IFT_C2P14=c(APre_IFT_C2P14, (1-(aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPre_IFT_14)
  APre_IFT_P2C19=c(APre_IFT_P2C19, (1-(aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPre_IFT_19)
  APre_IFT_C2P19=c(APre_IFT_C2P19, (1-(aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPre_IFT_19)
  APre_IFT_P2C24=c(APre_IFT_P2C24, (1-(aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPre_IFT_24)
  APre_IFT_C2P24=c(APre_IFT_C2P24, (1-(aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPre_IFT_24)
  
}

for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  RPre_Ipv_P2C9 = (sample(Pre_Ipv_P2C9 ,1, replace= T))
  RPre_Ipv_C2P9 = (sample(Pre_Ipv_C2P9 ,1, replace= T))
  RPre_Ipv_P2C14= (sample(Pre_Ipv_P2C14,1, replace= T))
  RPre_Ipv_C2P14= (sample(Pre_Ipv_C2P14,1, replace= T))
  RPre_Ipv_P2C19= (sample(Pre_Ipv_P2C19,1, replace= T))
  RPre_Ipv_C2P19= (sample(Pre_Ipv_C2P19,1, replace= T))
  RPre_Ipv_P2C24= (sample(Pre_Ipv_P2C24,1, replace= T))
  RPre_Ipv_C2P24= (sample(Pre_Ipv_C2P24,1, replace= T))
  
  APre_Ipv_P2C9 = c(APre_Ipv_P2C9 , (1-(aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPre_Ipv_P2C9 )
  APre_Ipv_C2P9 = c(APre_Ipv_C2P9 , (1-(aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPre_Ipv_C2P9 )
  APre_Ipv_P2C14= c(APre_Ipv_P2C14, (1-(aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPre_Ipv_P2C14)
  APre_Ipv_C2P14= c(APre_Ipv_C2P14, (1-(aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPre_Ipv_C2P14)
  APre_Ipv_P2C19= c(APre_Ipv_P2C19, (1-(aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPre_Ipv_P2C19)
  APre_Ipv_C2P19= c(APre_Ipv_C2P19, (1-(aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPre_Ipv_C2P19)
  APre_Ipv_P2C24= c(APre_Ipv_P2C24, (1-(aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPre_Ipv_P2C24)
  APre_Ipv_C2P24= c(APre_Ipv_C2P24, (1-(aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPre_Ipv_C2P24)
}
  
  
for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  RPre_Iss_P2C9 = (sample(Pre_Iss_P2C9 ,1, replace= T)) 
  RPre_Iss_C2P9 = (sample(Pre_Iss_C2P9 ,1, replace= T)) 
  RPre_Iss_P2C14= (sample(Pre_Iss_P2C14,1, replace= T)) 
  RPre_Iss_C2P14= (sample(Pre_Iss_C2P14,1, replace= T)) 
  RPre_Iss_P2C19= (sample(Pre_Iss_P2C19,1, replace= T)) 
  RPre_Iss_C2P19= (sample(Pre_Iss_C2P19,1, replace= T)) 
  RPre_Iss_P2C24= (sample(Pre_Iss_P2C24,1, replace= T))
  RPre_Iss_C2P24= (sample(Pre_Iss_C2P24,1, replace= T))
  
  APre_Iss_P2C9 =c(APre_Iss_P2C9 , (1-(aPre_Ipv_P2C9 +aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPre_Iss_P2C9 )
  APre_Iss_C2P9 =c(APre_Iss_C2P9 , (1-(aPre_Ipv_C2P9 +aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPre_Iss_C2P9 )
  APre_Iss_P2C14=c(APre_Iss_P2C14, (1-(aPre_Ipv_P2C14+aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPre_Iss_P2C14)
  APre_Iss_C2P14=c(APre_Iss_C2P14, (1-(aPre_Ipv_C2P14+aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPre_Iss_C2P14)
  APre_Iss_P2C19=c(APre_Iss_P2C19, (1-(aPre_Ipv_P2C19+aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPre_Iss_P2C19)
  APre_Iss_C2P19=c(APre_Iss_C2P19, (1-(aPre_Ipv_C2P19+aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPre_Iss_C2P19)
  APre_Iss_P2C24=c(APre_Iss_P2C24, (1-(aPre_Ipv_P2C24+aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPre_Iss_P2C24)
  APre_Iss_C2P24=c(APre_Iss_C2P24, (1-(aPre_Ipv_C2P24+aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPre_Iss_C2P24)
  
}

for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T)) 
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T)) 
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T)) 
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T)) 
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T)) 
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T)) 
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  
  RPos_Igerm_P2C9 = (sample(Pos_Igerm_P2C9 , 1, replace = T))
  RPos_Igerm_C2P9 = (sample(Pos_Igerm_C2P9 , 1, replace = T))
  RPos_Igerm_P2C14= (sample(Pos_Igerm_P2C14, 1, replace = T))
  RPos_Igerm_C2P14= (sample(Pos_Igerm_C2P14, 1, replace = T))
  RPos_Igerm_P2C19= (sample(Pos_Igerm_P2C19, 1, replace = T))
  RPos_Igerm_C2P19= (sample(Pos_Igerm_C2P19, 1, replace = T))
  RPos_Igerm_P2C24= (sample(Pos_Igerm_P2C24, 1, replace = T))
  RPos_Igerm_C2P24= (sample(Pos_Igerm_C2P24, 1, replace = T))
  
  APos_Igerm_P2C9 =c(APos_Igerm_P2C9 ,(1-(aPre_Iss_P2C9 +aPre_Ipv_P2C9 +aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Igerm_P2C9 )
  APos_Igerm_C2P9 =c(APos_Igerm_C2P9 ,(1-(aPre_Iss_C2P9 +aPre_Ipv_C2P9 +aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Igerm_C2P9 )
  APos_Igerm_P2C14=c(APos_Igerm_P2C14,(1-(aPre_Iss_P2C14+aPre_Ipv_P2C14+aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Igerm_P2C14)
  APos_Igerm_C2P14=c(APos_Igerm_C2P14,(1-(aPre_Iss_C2P14+aPre_Ipv_C2P14+aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Igerm_C2P14)
  APos_Igerm_P2C19=c(APos_Igerm_P2C19,(1-(aPre_Iss_P2C19+aPre_Ipv_P2C19+aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Igerm_P2C19)
  APos_Igerm_C2P19=c(APos_Igerm_C2P19,(1-(aPre_Iss_C2P19+aPre_Ipv_C2P19+aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Igerm_C2P19)
  APos_Igerm_P2C24=c(APos_Igerm_P2C24,(1-(aPre_Iss_P2C24+aPre_Ipv_P2C24+aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Igerm_P2C24)
  APos_Igerm_C2P24=c(APos_Igerm_C2P24,(1-(aPre_Iss_C2P24+aPre_Ipv_C2P24+aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Igerm_C2P24)
}  
  
for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T)) 
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T)) 
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T)) 
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T)) 
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T)) 
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T)) 
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  
  aPos_Igerm_P2C9 = (sample(APos_Igerm_P2C9 , 1, replace = T))
  aPos_Igerm_C2P9 = (sample(APos_Igerm_C2P9 , 1, replace = T))
  aPos_Igerm_P2C14= (sample(APos_Igerm_P2C14, 1, replace = T))
  aPos_Igerm_C2P14= (sample(APos_Igerm_C2P14, 1, replace = T))
  aPos_Igerm_P2C19= (sample(APos_Igerm_P2C19, 1, replace = T))
  aPos_Igerm_C2P19= (sample(APos_Igerm_C2P19, 1, replace = T))
  aPos_Igerm_P2C24= (sample(APos_Igerm_P2C24, 1, replace = T))
  aPos_Igerm_C2P24= (sample(APos_Igerm_C2P24, 1, replace = T))
  
  RPos_Isurv_P2C9  = (sample(Pos_Isurv_P2C9  , 1, replace = T))
  RPos_Isurv_C2P9  = (sample(Pos_Isurv_C2P9  , 1, replace = T))
  RPos_Isurv_P2C14 = (sample(Pos_Isurv_P2C14 , 1, replace = T))
  RPos_Isurv_C2P14 = (sample(Pos_Isurv_C2P14 , 1, replace = T))
  RPos_Isurv_P2C19 = (sample(Pos_Isurv_P2C19 , 1, replace = T))
  RPos_Isurv_C2P19 = (sample(Pos_Isurv_C2P19 , 1, replace = T))
  RPos_Isurv_P2C24 = (sample(Pos_Isurv_P2C24 , 1, replace = T))
  RPos_Isurv_C2P24 = (sample(Pos_Isurv_C2P24 , 1, replace = T))
  
  APos_Isurv_P2C9  = c(APos_Isurv_P2C9 , (1-(aPos_Igerm_P2C9 +aPre_Iss_P2C9 +aPre_Ipv_P2C9 +aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Isurv_P2C9 )
  APos_Isurv_C2P9  = c(APos_Isurv_C2P9 , (1-(aPos_Igerm_C2P9 +aPre_Iss_C2P9 +aPre_Ipv_C2P9 +aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Isurv_C2P9 )
  APos_Isurv_P2C14 = c(APos_Isurv_P2C14, (1-(aPos_Igerm_P2C14+aPre_Iss_P2C14+aPre_Ipv_P2C14+aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Isurv_P2C14)
  APos_Isurv_C2P14 = c(APos_Isurv_C2P14, (1-(aPos_Igerm_C2P14+aPre_Iss_C2P14+aPre_Ipv_C2P14+aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Isurv_C2P14)
  APos_Isurv_P2C19 = c(APos_Isurv_P2C19, (1-(aPos_Igerm_P2C19+aPre_Iss_P2C19+aPre_Ipv_P2C19+aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Isurv_P2C19)
  APos_Isurv_C2P19 = c(APos_Isurv_C2P19, (1-(aPos_Igerm_C2P19+aPre_Iss_C2P19+aPre_Ipv_C2P19+aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Isurv_C2P19)
  APos_Isurv_P2C24 = c(APos_Isurv_P2C24, (1-(aPos_Igerm_P2C24+aPre_Iss_P2C24+aPre_Ipv_P2C24+aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Isurv_P2C24)
  APos_Isurv_C2P24 = c(APos_Isurv_C2P24, (1-(aPos_Igerm_C2P24+aPre_Iss_C2P24+aPre_Ipv_C2P24+aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Isurv_C2P24)
}

for (i in 1:10000) {  
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T)) 
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T)) 
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T)) 
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T)) 
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T)) 
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T)) 
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  
  aPos_Igerm_P2C9 = (sample(APos_Igerm_P2C9 , 1, replace = T))
  aPos_Igerm_C2P9 = (sample(APos_Igerm_C2P9 , 1, replace = T))
  aPos_Igerm_P2C14= (sample(APos_Igerm_P2C14, 1, replace = T))
  aPos_Igerm_C2P14= (sample(APos_Igerm_C2P14, 1, replace = T))
  aPos_Igerm_P2C19= (sample(APos_Igerm_P2C19, 1, replace = T))
  aPos_Igerm_C2P19= (sample(APos_Igerm_C2P19, 1, replace = T))
  aPos_Igerm_P2C24= (sample(APos_Igerm_P2C24, 1, replace = T))
  aPos_Igerm_C2P24= (sample(APos_Igerm_C2P24, 1, replace = T))
  
  aPos_Isurv_P2C9  = (sample(APos_Isurv_P2C9  , 1, replace = T))
  aPos_Isurv_C2P9  = (sample(APos_Isurv_C2P9  , 1, replace = T))
  aPos_Isurv_P2C14 = (sample(APos_Isurv_P2C14 , 1, replace = T))
  aPos_Isurv_C2P14 = (sample(APos_Isurv_C2P14 , 1, replace = T))
  aPos_Isurv_P2C19 = (sample(APos_Isurv_P2C19 , 1, replace = T))
  aPos_Isurv_C2P19 = (sample(APos_Isurv_C2P19 , 1, replace = T))
  aPos_Isurv_P2C24 = (sample(APos_Isurv_P2C24 , 1, replace = T))
  aPos_Isurv_C2P24 = (sample(APos_Isurv_C2P24 , 1, replace = T))
  
  RPos_IFT_P2C9 = (sample(Pos_IFT_P2C9 , 1, replace=T))
  RPos_IFT_C2P9 = (sample(Pos_IFT_C2P9 , 1, replace=T))
  RPos_IFT_P2C14= (sample(Pos_IFT_P2C14, 1, replace=T))
  RPos_IFT_C2P14= (sample(Pos_IFT_C2P14, 1, replace=T))
  RPos_IFT_P2C19= (sample(Pos_IFT_P2C19, 1, replace=T))
  RPos_IFT_C2P19= (sample(Pos_IFT_C2P19, 1, replace=T))
  RPos_IFT_P2C24= (sample(Pos_IFT_P2C24, 1, replace=T))
  RPos_IFT_C2P24= (sample(Pos_IFT_C2P24, 1, replace=T))
  
  APos_IFT_P2C9 =c(APos_IFT_P2C9 , (1-(aPos_Isurv_P2C9 +aPos_Igerm_P2C9 +aPre_Iss_P2C9 +aPre_Ipv_P2C9 +aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_IFT_P2C9 )
  APos_IFT_C2P9 =c(APos_IFT_C2P9 , (1-(aPos_Isurv_C2P9 +aPos_Igerm_C2P9 +aPre_Iss_C2P9 +aPre_Ipv_C2P9 +aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_IFT_C2P9 )
  APos_IFT_P2C14=c(APos_IFT_P2C14, (1-(aPos_Isurv_P2C14+aPos_Igerm_P2C14+aPre_Iss_P2C14+aPre_Ipv_P2C14+aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_IFT_P2C14)
  APos_IFT_C2P14=c(APos_IFT_C2P14, (1-(aPos_Isurv_C2P14+aPos_Igerm_C2P14+aPre_Iss_C2P14+aPre_Ipv_C2P14+aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_IFT_C2P14)
  APos_IFT_P2C19=c(APos_IFT_P2C19, (1-(aPos_Isurv_P2C19+aPos_Igerm_P2C19+aPre_Iss_P2C19+aPre_Ipv_P2C19+aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_IFT_P2C19)
  APos_IFT_C2P19=c(APos_IFT_C2P19, (1-(aPos_Isurv_C2P19+aPos_Igerm_C2P19+aPre_Iss_C2P19+aPre_Ipv_C2P19+aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_IFT_C2P19)
  APos_IFT_P2C24=c(APos_IFT_P2C24, (1-(aPos_Isurv_P2C24+aPos_Igerm_P2C24+aPre_Iss_P2C24+aPre_Ipv_P2C24+aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_IFT_P2C24)
  APos_IFT_C2P24=c(APos_IFT_C2P24, (1-(aPos_Isurv_C2P24+aPos_Igerm_C2P24+aPre_Iss_C2P24+aPre_Ipv_C2P24+aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_IFT_C2P24)
}

for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T)) 
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T)) 
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T)) 
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T)) 
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T)) 
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T)) 
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  
  aPos_Igerm_P2C9 = (sample(APos_Igerm_P2C9 , 1, replace = T))
  aPos_Igerm_C2P9 = (sample(APos_Igerm_C2P9 , 1, replace = T))
  aPos_Igerm_P2C14= (sample(APos_Igerm_P2C14, 1, replace = T))
  aPos_Igerm_C2P14= (sample(APos_Igerm_C2P14, 1, replace = T))
  aPos_Igerm_P2C19= (sample(APos_Igerm_P2C19, 1, replace = T))
  aPos_Igerm_C2P19= (sample(APos_Igerm_C2P19, 1, replace = T))
  aPos_Igerm_P2C24= (sample(APos_Igerm_P2C24, 1, replace = T))
  aPos_Igerm_C2P24= (sample(APos_Igerm_C2P24, 1, replace = T))
  
  aPos_Isurv_P2C9  = (sample(APos_Isurv_P2C9  , 1, replace = T))
  aPos_Isurv_C2P9  = (sample(APos_Isurv_C2P9  , 1, replace = T))
  aPos_Isurv_P2C14 = (sample(APos_Isurv_P2C14 , 1, replace = T))
  aPos_Isurv_C2P14 = (sample(APos_Isurv_C2P14 , 1, replace = T))
  aPos_Isurv_P2C19 = (sample(APos_Isurv_P2C19 , 1, replace = T))
  aPos_Isurv_C2P19 = (sample(APos_Isurv_C2P19 , 1, replace = T))
  aPos_Isurv_P2C24 = (sample(APos_Isurv_P2C24 , 1, replace = T))
  aPos_Isurv_C2P24 = (sample(APos_Isurv_C2P24 , 1, replace = T))
  
  aPos_IFT_P2C9 = (sample(APos_IFT_P2C9 , 1, replace=T))
  aPos_IFT_C2P9 = (sample(APos_IFT_C2P9 , 1, replace=T))
  aPos_IFT_P2C14= (sample(APos_IFT_P2C14, 1, replace=T))
  aPos_IFT_C2P14= (sample(APos_IFT_C2P14, 1, replace=T))
  aPos_IFT_P2C19= (sample(APos_IFT_P2C19, 1, replace=T))
  aPos_IFT_C2P19= (sample(APos_IFT_C2P19, 1, replace=T))
  aPos_IFT_P2C24= (sample(APos_IFT_P2C24, 1, replace=T))
  aPos_IFT_C2P24= (sample(APos_IFT_C2P24, 1, replace=T))
  
  RPos_Ipv_P2C9 = (sample(Pos_Ipv_P2C9 ,1, replace= T))
  RPos_Ipv_C2P9 = (sample(Pos_Ipv_C2P9 ,1, replace= T))
  RPos_Ipv_P2C14= (sample(Pos_Ipv_P2C14,1, replace= T))
  RPos_Ipv_C2P14= (sample(Pos_Ipv_C2P14,1, replace= T))
  RPos_Ipv_P2C19= (sample(Pos_Ipv_P2C19,1, replace= T))
  RPos_Ipv_C2P19= (sample(Pos_Ipv_C2P19,1, replace= T))
  RPos_Ipv_P2C24= (sample(Pos_Ipv_P2C24,1, replace= T))
  RPos_Ipv_C2P24= (sample(Pos_Ipv_C2P24,1, replace= T))
  
  APos_Ipv_P2C9 = c(APos_Ipv_P2C9 , (1-(aPos_IFT_P2C9 +aPos_Isurv_P2C9 +aPos_Igerm_P2C9 +aPre_Iss_P2C9 +aPre_Ipv_P2C9 +aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Ipv_P2C9 )
  APos_Ipv_C2P9 = c(APos_Ipv_C2P9 , (1-(aPos_IFT_C2P9 +aPos_Isurv_C2P9 +aPos_Igerm_C2P9 +aPre_Iss_C2P9 +aPre_Ipv_C2P9 +aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Ipv_C2P9 )
  APos_Ipv_P2C14= c(APos_Ipv_P2C14, (1-(aPos_IFT_P2C14+aPos_Isurv_P2C14+aPos_Igerm_P2C14+aPre_Iss_P2C14+aPre_Ipv_P2C14+aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Ipv_P2C14)
  APos_Ipv_C2P14= c(APos_Ipv_C2P14, (1-(aPos_IFT_C2P14+aPos_Isurv_C2P14+aPos_Igerm_C2P14+aPre_Iss_C2P14+aPre_Ipv_C2P14+aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Ipv_C2P14)
  APos_Ipv_P2C19= c(APos_Ipv_P2C19, (1-(aPos_IFT_P2C19+aPos_Isurv_P2C19+aPos_Igerm_P2C19+aPre_Iss_P2C19+aPre_Ipv_P2C19+aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Ipv_P2C19)
  APos_Ipv_C2P19= c(APos_Ipv_C2P19, (1-(aPos_IFT_C2P19+aPos_Isurv_C2P19+aPos_Igerm_C2P19+aPre_Iss_C2P19+aPre_Ipv_C2P19+aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Ipv_C2P19)
  APos_Ipv_P2C24= c(APos_Ipv_P2C24, (1-(aPos_IFT_P2C24+aPos_Isurv_P2C24+aPos_Igerm_P2C24+aPre_Iss_P2C24+aPre_Ipv_P2C24+aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Ipv_P2C24)
  APos_Ipv_C2P24= c(APos_Ipv_C2P24, (1-(aPos_IFT_C2P24+aPos_Isurv_C2P24+aPos_Igerm_C2P24+aPre_Iss_C2P24+aPre_Ipv_C2P24+aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Ipv_C2P24)
}
  
for (i in 1:10000) { 
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T)) 
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T)) 
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T)) 
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T)) 
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T)) 
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T)) 
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  
  aPos_Igerm_P2C9 = (sample(APos_Igerm_P2C9 , 1, replace = T))
  aPos_Igerm_C2P9 = (sample(APos_Igerm_C2P9 , 1, replace = T))
  aPos_Igerm_P2C14= (sample(APos_Igerm_P2C14, 1, replace = T))
  aPos_Igerm_C2P14= (sample(APos_Igerm_C2P14, 1, replace = T))
  aPos_Igerm_P2C19= (sample(APos_Igerm_P2C19, 1, replace = T))
  aPos_Igerm_C2P19= (sample(APos_Igerm_C2P19, 1, replace = T))
  aPos_Igerm_P2C24= (sample(APos_Igerm_P2C24, 1, replace = T))
  aPos_Igerm_C2P24= (sample(APos_Igerm_C2P24, 1, replace = T))
  
  aPos_Isurv_P2C9  = (sample(APos_Isurv_P2C9  , 1, replace = T))
  aPos_Isurv_C2P9  = (sample(APos_Isurv_C2P9  , 1, replace = T))
  aPos_Isurv_P2C14 = (sample(APos_Isurv_P2C14 , 1, replace = T))
  aPos_Isurv_C2P14 = (sample(APos_Isurv_C2P14 , 1, replace = T))
  aPos_Isurv_P2C19 = (sample(APos_Isurv_P2C19 , 1, replace = T))
  aPos_Isurv_C2P19 = (sample(APos_Isurv_C2P19 , 1, replace = T))
  aPos_Isurv_P2C24 = (sample(APos_Isurv_P2C24 , 1, replace = T))
  aPos_Isurv_C2P24 = (sample(APos_Isurv_C2P24 , 1, replace = T))
  
  aPos_IFT_P2C9 = (sample(APos_IFT_P2C9 , 1, replace=T))
  aPos_IFT_C2P9 = (sample(APos_IFT_C2P9 , 1, replace=T))
  aPos_IFT_P2C14= (sample(APos_IFT_P2C14, 1, replace=T))
  aPos_IFT_C2P14= (sample(APos_IFT_C2P14, 1, replace=T))
  aPos_IFT_P2C19= (sample(APos_IFT_P2C19, 1, replace=T))
  aPos_IFT_C2P19= (sample(APos_IFT_C2P19, 1, replace=T))
  aPos_IFT_P2C24= (sample(APos_IFT_P2C24, 1, replace=T))
  aPos_IFT_C2P24= (sample(APos_IFT_C2P24, 1, replace=T))
  
  aPos_Ipv_P2C9 = (sample(APos_Ipv_P2C9 ,1, replace= T))
  aPos_Ipv_C2P9 = (sample(APos_Ipv_C2P9 ,1, replace= T))
  aPos_Ipv_P2C14= (sample(APos_Ipv_P2C14,1, replace= T))
  aPos_Ipv_C2P14= (sample(APos_Ipv_C2P14,1, replace= T))
  aPos_Ipv_P2C19= (sample(APos_Ipv_P2C19,1, replace= T))
  aPos_Ipv_C2P19= (sample(APos_Ipv_C2P19,1, replace= T))
  aPos_Ipv_P2C24= (sample(APos_Ipv_P2C24,1, replace= T))
  aPos_Ipv_C2P24= (sample(APos_Ipv_C2P24,1, replace= T))
  
  RPos_Iss_P2C9 = (sample(Pos_Iss_P2C9 ,1, replace= T)) 
  RPos_Iss_C2P9 = (sample(Pos_Iss_C2P9 ,1, replace= T)) 
  RPos_Iss_P2C14= (sample(Pos_Iss_P2C14,1, replace= T)) 
  RPos_Iss_C2P14= (sample(Pos_Iss_C2P14,1, replace= T)) 
  RPos_Iss_P2C19= (sample(Pos_Iss_P2C19,1, replace= T)) 
  RPos_Iss_C2P19= (sample(Pos_Iss_C2P19,1, replace= T)) 
  RPos_Iss_P2C24= (sample(Pos_Iss_P2C24,1, replace= T))
  RPos_Iss_C2P24= (sample(Pos_Iss_C2P24,1, replace= T))
  
  APos_Iss_P2C9 =c(APos_Iss_P2C9 , (1-(aPos_Ipv_P2C9 +aPos_IFT_P2C9 +aPos_Isurv_P2C9 +aPos_Igerm_P2C9 +aPre_Iss_P2C9 +aPre_Ipv_P2C9 +aPre_IFT_P2C9 +aPre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Iss_P2C9 )
  APos_Iss_C2P9 =c(APos_Iss_C2P9 , (1-(aPos_Ipv_C2P9 +aPos_IFT_C2P9 +aPos_Isurv_C2P9 +aPos_Igerm_C2P9 +aPre_Iss_C2P9 +aPre_Ipv_C2P9 +aPre_IFT_C2P9 +aPre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Iss_C2P9 )
  APos_Iss_P2C14=c(APos_Iss_P2C14, (1-(aPos_Ipv_P2C14+aPos_IFT_P2C14+aPos_Isurv_P2C14+aPos_Igerm_P2C14+aPre_Iss_P2C14+aPre_Ipv_P2C14+aPre_IFT_P2C14+aPre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Iss_P2C14)
  APos_Iss_C2P14=c(APos_Iss_C2P14, (1-(aPos_Ipv_C2P14+aPos_IFT_C2P14+aPos_Isurv_C2P14+aPos_Igerm_C2P14+aPre_Iss_C2P14+aPre_Ipv_C2P14+aPre_IFT_C2P14+aPre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Iss_C2P14)
  APos_Iss_P2C19=c(APos_Iss_P2C19, (1-(aPos_Ipv_P2C19+aPos_IFT_P2C19+aPos_Isurv_P2C19+aPos_Igerm_P2C19+aPre_Iss_P2C19+aPre_Ipv_P2C19+aPre_IFT_P2C19+aPre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Iss_P2C19)
  APos_Iss_C2P19=c(APos_Iss_C2P19, (1-(aPos_Ipv_C2P19+aPos_IFT_C2P19+aPos_Isurv_C2P19+aPos_Igerm_C2P19+aPre_Iss_C2P19+aPre_Ipv_C2P19+aPre_IFT_C2P19+aPre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Iss_C2P19)
  APos_Iss_P2C24=c(APos_Iss_P2C24, (1-(aPos_Ipv_P2C24+aPos_IFT_P2C24+aPos_Isurv_P2C24+aPos_Igerm_P2C24+aPre_Iss_P2C24+aPre_Ipv_P2C24+aPre_IFT_P2C24+aPre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Iss_P2C24)
  APos_Iss_C2P24=c(APos_Iss_C2P24, (1-(aPos_Ipv_C2P24+aPos_IFT_C2P24+aPos_Isurv_C2P24+aPos_Igerm_C2P24+aPre_Iss_C2P24+aPre_Ipv_C2P24+aPre_IFT_C2P24+aPre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Iss_C2P24)
}

set.seed(1)

aItotal_P2C9 =vector()
aItotal_C2P9 =vector()
aItotal_P2C14=vector()
aItotal_C2P14=vector()
aItotal_P2C19=vector()
aItotal_C2P19=vector()
aItotal_P2C24=vector()
aItotal_C2P24=vector()

for (i in 1:1000000) {
  
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  aPre_Isurv_P2C9  = (sample(APre_Isurv_P2C9  , 1, replace = T))
  aPre_Isurv_C2P9  = (sample(APre_Isurv_C2P9  , 1, replace = T))
  aPre_Isurv_P2C14 = (sample(APre_Isurv_P2C14 , 1, replace = T))
  aPre_Isurv_C2P14 = (sample(APre_Isurv_C2P14 , 1, replace = T))
  aPre_Isurv_P2C19 = (sample(APre_Isurv_P2C19 , 1, replace = T))
  aPre_Isurv_C2P19 = (sample(APre_Isurv_C2P19 , 1, replace = T))
  aPre_Isurv_P2C24 = (sample(APre_Isurv_P2C24 , 1, replace = T))
  aPre_Isurv_C2P24 = (sample(APre_Isurv_C2P24 , 1, replace = T))
  
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 , 1, replace=T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 , 1, replace=T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14, 1, replace=T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14, 1, replace=T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19, 1, replace=T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19, 1, replace=T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24, 1, replace=T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24, 1, replace=T))
  
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T)) 
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T)) 
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T)) 
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T)) 
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T)) 
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T)) 
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  
  aPos_Igerm_P2C9 = (sample(APos_Igerm_P2C9 , 1, replace = T))
  aPos_Igerm_C2P9 = (sample(APos_Igerm_C2P9 , 1, replace = T))
  aPos_Igerm_P2C14= (sample(APos_Igerm_P2C14, 1, replace = T))
  aPos_Igerm_C2P14= (sample(APos_Igerm_C2P14, 1, replace = T))
  aPos_Igerm_P2C19= (sample(APos_Igerm_P2C19, 1, replace = T))
  aPos_Igerm_C2P19= (sample(APos_Igerm_C2P19, 1, replace = T))
  aPos_Igerm_P2C24= (sample(APos_Igerm_P2C24, 1, replace = T))
  aPos_Igerm_C2P24= (sample(APos_Igerm_C2P24, 1, replace = T))
  
  aPos_Isurv_P2C9  = (sample(APos_Isurv_P2C9  , 1, replace = T))
  aPos_Isurv_C2P9  = (sample(APos_Isurv_C2P9  , 1, replace = T))
  aPos_Isurv_P2C14 = (sample(APos_Isurv_P2C14 , 1, replace = T))
  aPos_Isurv_C2P14 = (sample(APos_Isurv_C2P14 , 1, replace = T))
  aPos_Isurv_P2C19 = (sample(APos_Isurv_P2C19 , 1, replace = T))
  aPos_Isurv_C2P19 = (sample(APos_Isurv_C2P19 , 1, replace = T))
  aPos_Isurv_P2C24 = (sample(APos_Isurv_P2C24 , 1, replace = T))
  aPos_Isurv_C2P24 = (sample(APos_Isurv_C2P24 , 1, replace = T))
  
  aPos_IFT_P2C9 = (sample(APos_IFT_P2C9 , 1, replace=T))
  aPos_IFT_C2P9 = (sample(APos_IFT_C2P9 , 1, replace=T))
  aPos_IFT_P2C14= (sample(APos_IFT_P2C14, 1, replace=T))
  aPos_IFT_C2P14= (sample(APos_IFT_C2P14, 1, replace=T))
  aPos_IFT_P2C19= (sample(APos_IFT_P2C19, 1, replace=T))
  aPos_IFT_C2P19= (sample(APos_IFT_C2P19, 1, replace=T))
  aPos_IFT_P2C24= (sample(APos_IFT_P2C24, 1, replace=T))
  aPos_IFT_C2P24= (sample(APos_IFT_C2P24, 1, replace=T))
  
  aPos_Ipv_P2C9 = (sample(APos_Ipv_P2C9 ,1, replace= T))
  aPos_Ipv_C2P9 = (sample(APos_Ipv_C2P9 ,1, replace= T))
  aPos_Ipv_P2C14= (sample(APos_Ipv_P2C14,1, replace= T))
  aPos_Ipv_C2P14= (sample(APos_Ipv_C2P14,1, replace= T))
  aPos_Ipv_P2C19= (sample(APos_Ipv_P2C19,1, replace= T))
  aPos_Ipv_C2P19= (sample(APos_Ipv_C2P19,1, replace= T))
  aPos_Ipv_P2C24= (sample(APos_Ipv_P2C24,1, replace= T))
  aPos_Ipv_C2P24= (sample(APos_Ipv_C2P24,1, replace= T))
  
  aPos_Iss_P2C9 = (sample(APos_Iss_P2C9 ,1, replace= T)) 
  aPos_Iss_C2P9 = (sample(APos_Iss_C2P9 ,1, replace= T)) 
  aPos_Iss_P2C14= (sample(APos_Iss_P2C14,1, replace= T)) 
  aPos_Iss_C2P14= (sample(APos_Iss_C2P14,1, replace= T)) 
  aPos_Iss_P2C19= (sample(APos_Iss_P2C19,1, replace= T)) 
  aPos_Iss_C2P19= (sample(APos_Iss_C2P19,1, replace= T)) 
  aPos_Iss_P2C24= (sample(APos_Iss_P2C24,1, replace= T))
  aPos_Iss_C2P24= (sample(APos_Iss_C2P24,1, replace= T))
  
  aItotal_P2C9 = c(aItotal_P2C9 ,(RPre_Igerm_P2C9 +aPre_Isurv_P2C9 +aPre_IFT_P2C9 +aPre_Ipv_P2C9 +aPre_Iss_P2C9 +aPos_Igerm_P2C9 +aPos_Isurv_P2C9 +aPos_IFT_P2C9 +aPos_Ipv_P2C9 +aPos_Iss_P2C9 ))
  aItotal_C2P9 = c(aItotal_C2P9 ,(RPre_Igerm_C2P9 +aPre_Isurv_C2P9 +aPre_IFT_C2P9 +aPre_Ipv_C2P9 +aPre_Iss_C2P9 +aPos_Igerm_C2P9 +aPos_Isurv_C2P9 +aPos_IFT_C2P9 +aPos_Ipv_C2P9 +aPos_Iss_C2P9 ))
  aItotal_P2C14= c(aItotal_P2C14,(RPre_Igerm_P2C14+aPre_Isurv_P2C14+aPre_IFT_P2C14+aPre_Ipv_P2C14+aPre_Iss_P2C14+aPos_Igerm_P2C14+aPos_Isurv_P2C14+aPos_IFT_P2C14+aPos_Ipv_P2C14+aPos_Iss_P2C14))
  aItotal_C2P14= c(aItotal_C2P14,(RPre_Igerm_C2P14+aPre_Isurv_C2P14+aPre_IFT_C2P14+aPre_Ipv_C2P14+aPre_Iss_C2P14+aPos_Igerm_C2P14+aPos_Isurv_C2P14+aPos_IFT_C2P14+aPos_Ipv_C2P14+aPos_Iss_C2P14))
  aItotal_P2C19= c(aItotal_P2C19,(RPre_Igerm_P2C19+aPre_Isurv_P2C19+aPre_IFT_P2C19+aPre_Ipv_P2C19+aPre_Iss_P2C19+aPos_Igerm_P2C19+aPos_Isurv_P2C19+aPos_IFT_P2C19+aPos_Ipv_P2C19+aPos_Iss_P2C19))
  aItotal_C2P19= c(aItotal_C2P19,(RPre_Igerm_C2P19+aPre_Isurv_C2P19+aPre_IFT_C2P19+aPre_Ipv_C2P19+aPre_Iss_C2P19+aPos_Igerm_C2P19+aPos_Isurv_C2P19+aPos_IFT_C2P19+aPos_Ipv_C2P19+aPos_Iss_C2P19))
  aItotal_P2C24= c(aItotal_P2C24,(RPre_Igerm_P2C24+aPre_Isurv_P2C24+aPre_IFT_P2C24+aPre_Ipv_P2C24+aPre_Iss_P2C24+aPos_Igerm_P2C24+aPos_Isurv_P2C24+aPos_IFT_P2C24+aPos_Ipv_P2C24+aPos_Iss_P2C24))
  aItotal_C2P24= c(aItotal_C2P24,(RPre_Igerm_C2P24+aPre_Isurv_C2P24+aPre_IFT_C2P24+aPre_Ipv_C2P24+aPre_Iss_C2P24+aPos_Igerm_C2P24+aPos_Isurv_C2P24+aPos_IFT_C2P24+aPos_Ipv_C2P24+aPos_Iss_C2P24))
  
}

summary(sort(aItotal_P2C9 )[c(25000, 975000)])
summary(sort(aItotal_C2P9 )[c(25000, 975000)])
summary(sort(aItotal_P2C14)[c(25000, 975000)])
summary(sort(aItotal_C2P14)[c(25000, 975000)])
summary(sort(aItotal_P2C19)[c(25000, 975000)])
summary(sort(aItotal_C2P19)[c(25000, 975000)])
summary(sort(aItotal_P2C24)[c(25000, 975000)])
summary(sort(aItotal_C2P24)[c(25000, 975000)])




for (i in 1:10000) {
  RPre_Igerm_P2C9 = (sample(RPre_Igerm_P2C9 ,1, replace= T))
  aPre_Isurv_P2C9 = (sample(APre_Isurv_P2C9 ,1, replace= T))
  aPre_IFT_P2C9 = (sample(APre_IFT_P2C9 ,1, replace= T))
  aPre_Ipv_P2C9 = (sample(APre_Ipv_P2C9 ,1, replace= T))
  aPre_Iss_P2C9 = (sample(APre_Iss_P2C9 ,1, replace= T))
  aPos_Igerm_P2C9 = (sample(APos_Igerm_P2C9 ,1, replace= T))
  aPos_Isurv_P2C9 = (sample(APos_Isurv_P2C9 ,1, replace= T))
  aPos_IFT_P2C9 = (sample(APos_IFT_P2C9 ,1, replace= T))
  aPos_Ipv_P2C9 = (sample(APos_Ipv_P2C9 ,1, replace= T))
  aPos_Iss_P2C9 = (sample(APos_Iss_P2C9 ,1, replace= T))
  RPre_Igerm_C2P9 = (sample(RPre_Igerm_C2P9 ,1, replace= T))
  aPre_Isurv_C2P9 = (sample(APre_Isurv_C2P9 ,1, replace= T))
  aPre_IFT_C2P9 = (sample(APre_IFT_C2P9 ,1, replace= T))
  aPre_Ipv_C2P9 = (sample(APre_Ipv_C2P9 ,1, replace= T))
  aPre_Iss_C2P9 = (sample(APre_Iss_C2P9 ,1, replace= T))
  aPos_Igerm_C2P9 = (sample(APos_Igerm_C2P9 ,1, replace= T))
  aPos_Isurv_C2P9 = (sample(APos_Isurv_C2P9 ,1, replace= T))
  aPos_IFT_C2P9 = (sample(APos_IFT_C2P9 ,1, replace= T))
  aPos_Ipv_C2P9 = (sample(APos_Ipv_C2P9 ,1, replace= T))
  aPos_Iss_C2P9 = (sample(APos_Iss_C2P9 ,1, replace= T))
  RPre_Igerm_P2C14= (sample(RPre_Igerm_P2C14,1, replace= T))
  aPre_Isurv_P2C14= (sample(APre_Isurv_P2C14,1, replace= T))
  aPre_IFT_P2C14= (sample(APre_IFT_P2C14,1, replace= T))
  aPre_Ipv_P2C14= (sample(APre_Ipv_P2C14,1, replace= T))
  aPre_Iss_P2C14= (sample(APre_Iss_P2C14,1, replace= T))
  aPos_Igerm_P2C14= (sample(APos_Igerm_P2C14,1, replace= T))
  aPos_Isurv_P2C14= (sample(APos_Isurv_P2C14,1, replace= T))
  aPos_IFT_P2C14= (sample(APos_IFT_P2C14,1, replace= T))
  aPos_Ipv_P2C14= (sample(APos_Ipv_P2C14,1, replace= T))
  aPos_Iss_P2C14= (sample(APos_Iss_P2C14,1, replace= T))
  RPre_Igerm_C2P14= (sample(RPre_Igerm_C2P14,1, replace= T))
  aPre_Isurv_C2P14= (sample(APre_Isurv_C2P14,1, replace= T))
  aPre_IFT_C2P14= (sample(APre_IFT_C2P14,1, replace= T))
  aPre_Ipv_C2P14= (sample(APre_Ipv_C2P14,1, replace= T))
  aPre_Iss_C2P14= (sample(APre_Iss_C2P14,1, replace= T))
  aPos_Igerm_C2P14= (sample(APos_Igerm_C2P14,1, replace= T))
  aPos_Isurv_C2P14= (sample(APos_Isurv_C2P14,1, replace= T))
  aPos_IFT_C2P14= (sample(APos_IFT_C2P14,1, replace= T))
  aPos_Ipv_C2P14= (sample(APos_Ipv_C2P14,1, replace= T))
  aPos_Iss_C2P14= (sample(APos_Iss_C2P14,1, replace= T))
  RPre_Igerm_P2C19= (sample(RPre_Igerm_P2C19,1, replace= T))
  aPre_Isurv_P2C19= (sample(APre_Isurv_P2C19,1, replace= T))
  aPre_IFT_P2C19= (sample(APre_IFT_P2C19,1, replace= T))
  aPre_Ipv_P2C19= (sample(APre_Ipv_P2C19,1, replace= T))
  aPre_Iss_P2C19= (sample(APre_Iss_P2C19,1, replace= T))
  aPos_Igerm_P2C19= (sample(APos_Igerm_P2C19,1, replace= T))
  aPos_Isurv_P2C19= (sample(APos_Isurv_P2C19,1, replace= T))
  aPos_IFT_P2C19= (sample(APos_IFT_P2C19,1, replace= T))
  aPos_Ipv_P2C19= (sample(APos_Ipv_P2C19,1, replace= T))
  aPos_Iss_P2C19= (sample(APos_Iss_P2C19,1, replace= T))
  RPre_Igerm_C2P19= (sample(RPre_Igerm_C2P19,1, replace= T))
  aPre_Isurv_C2P19= (sample(APre_Isurv_C2P19,1, replace= T))
  aPre_IFT_C2P19= (sample(APre_IFT_C2P19,1, replace= T))
  aPre_Ipv_C2P19= (sample(APre_Ipv_C2P19,1, replace= T))
  aPre_Iss_C2P19= (sample(APre_Iss_C2P19,1, replace= T))
  aPos_Igerm_C2P19= (sample(APos_Igerm_C2P19,1, replace= T))
  aPos_Isurv_C2P19= (sample(APos_Isurv_C2P19,1, replace= T))
  aPos_IFT_C2P19= (sample(APos_IFT_C2P19,1, replace= T))
  aPos_Ipv_C2P19= (sample(APos_Ipv_C2P19,1, replace= T))
  aPos_Iss_C2P19= (sample(APos_Iss_C2P19,1, replace= T))
  RPre_Igerm_P2C24= (sample(RPre_Igerm_P2C24,1, replace= T))
  aPre_Isurv_P2C24= (sample(APre_Isurv_P2C24,1, replace= T))
  aPre_IFT_P2C24= (sample(APre_IFT_P2C24,1, replace= T))
  aPre_Ipv_P2C24= (sample(APre_Ipv_P2C24,1, replace= T))
  aPre_Iss_P2C24= (sample(APre_Iss_P2C24,1, replace= T))
  aPos_Igerm_P2C24= (sample(APos_Igerm_P2C24,1, replace= T))
  aPos_Isurv_P2C24= (sample(APos_Isurv_P2C24,1, replace= T))
  aPos_IFT_P2C24= (sample(APos_IFT_P2C24,1, replace= T))
  aPos_Ipv_P2C24= (sample(APos_Ipv_P2C24,1, replace= T))
  aPos_Iss_P2C24= (sample(APos_Iss_P2C24,1, replace= T))
  RPre_Igerm_C2P24= (sample(RPre_Igerm_C2P24,1, replace= T))
  aPre_Isurv_C2P24= (sample(APre_Isurv_C2P24,1, replace= T))
  aPre_IFT_C2P24= (sample(APre_IFT_C2P24,1, replace= T))
  aPre_Ipv_C2P24= (sample(APre_Ipv_C2P24,1, replace= T))
  aPre_Iss_C2P24= (sample(APre_Iss_C2P24,1, replace= T))
  aPos_Igerm_C2P24= (sample(APos_Igerm_C2P24,1, replace= T))
  aPos_Isurv_C2P24= (sample(APos_Isurv_C2P24,1, replace= T))
  aPos_IFT_C2P24= (sample(APos_IFT_C2P24,1, replace= T))
  aPos_Ipv_C2P24= (sample(APos_Ipv_C2P24,1, replace= T))
  aPos_Iss_C2P24= (sample(APos_Iss_C2P24,1, replace= T))
  
  Itotal_P2C9 = c(Itotal_P2C9 ,(RPre_Igerm_P2C9 +aPre_Isurv_P2C9 +aPre_IFT_P2C9 +aPre_Ipv_P2C9 +aPre_Iss_P2C9 +aPos_Igerm_P2C9 +aPos_Isurv_P2C9 +aPos_IFT_P2C9 +aPos_Ipv_P2C9 +aPos_Iss_P2C9 ))
  Itotal_C2P9 = c(Itotal_C2P9 ,(RPre_Igerm_C2P9 +aPre_Isurv_C2P9 +aPre_IFT_C2P9 +aPre_Ipv_C2P9 +aPre_Iss_C2P9 +aPos_Igerm_C2P9 +aPos_Isurv_C2P9 +aPos_IFT_C2P9 +aPos_Ipv_C2P9 +aPos_Iss_C2P9 ))
  Itotal_P2C14= c(Itotal_P2C14,(RPre_Igerm_P2C14+aPre_Isurv_P2C14+aPre_IFT_P2C14+aPre_Ipv_P2C14+aPre_Iss_P2C14+aPos_Igerm_P2C14+aPos_Isurv_P2C14+aPos_IFT_P2C14+aPos_Ipv_P2C14+aPos_Iss_P2C14))
  Itotal_C2P14= c(Itotal_C2P14,(RPre_Igerm_C2P14+aPre_Isurv_C2P14+aPre_IFT_C2P14+aPre_Ipv_C2P14+aPre_Iss_C2P14+aPos_Igerm_C2P14+aPos_Isurv_C2P14+aPos_IFT_C2P14+aPos_Ipv_C2P14+aPos_Iss_C2P14))
  Itotal_P2C19= c(Itotal_P2C19,(RPre_Igerm_P2C19+aPre_Isurv_P2C19+aPre_IFT_P2C19+aPre_Ipv_P2C19+aPre_Iss_P2C19+aPos_Igerm_P2C19+aPos_Isurv_P2C19+aPos_IFT_P2C19+aPos_Ipv_P2C19+aPos_Iss_P2C19))
  Itotal_C2P19= c(Itotal_C2P19,(RPre_Igerm_C2P19+aPre_Isurv_C2P19+aPre_IFT_C2P19+aPre_Ipv_C2P19+aPre_Iss_C2P19+aPos_Igerm_C2P19+aPos_Isurv_C2P19+aPos_IFT_C2P19+aPos_Ipv_C2P19+aPos_Iss_C2P19))
  Itotal_P2C24= c(Itotal_P2C24,(RPre_Igerm_P2C24+aPre_Isurv_P2C24+aPre_IFT_P2C24+aPre_Ipv_P2C24+aPre_Iss_P2C24+aPos_Igerm_P2C24+aPos_Isurv_P2C24+aPos_IFT_P2C24+aPos_Ipv_P2C24+aPos_Iss_P2C24))
  Itotal_C2P24= c(Itotal_C2P24,(RPre_Igerm_C2P24+aPre_Isurv_C2P24+aPre_IFT_C2P24+aPre_Ipv_C2P24+aPre_Iss_C2P24+aPos_Igerm_C2P24+aPos_Isurv_C2P24+aPos_IFT_C2P24+aPos_Ipv_C2P24+aPos_Iss_C2P24))
}

#####TOTAL#####

##Bootstrap CIs##
set.seed(1)

aItotal_P2C9  = vector()
aItotal_C2P9  = vector()
aItotal_P2C14 = vector()
aItotal_C2P14 = vector()
aItotal_P2C19 = vector()
aItotal_C2P19 = vector()
aItotal_P2C24 = vector()
aItotal_C2P24 = vector()


for (i in 1:10000) {
  
  RPre_Igerm_P2C9 = (sample(Pre_Igerm_P2C9 , 1, replace = T))
  RPre_Igerm_C2P9 = (sample(Pre_Igerm_C2P9 , 1, replace = T))
  RPre_Igerm_P2C14= (sample(Pre_Igerm_P2C14, 1, replace = T))
  RPre_Igerm_C2P14= (sample(Pre_Igerm_C2P14, 1, replace = T))
  RPre_Igerm_P2C19= (sample(Pre_Igerm_P2C19, 1, replace = T))
  RPre_Igerm_C2P19= (sample(Pre_Igerm_C2P19, 1, replace = T))
  RPre_Igerm_P2C24= (sample(Pre_Igerm_P2C24, 1, replace = T))
  RPre_Igerm_C2P24= (sample(Pre_Igerm_C2P24, 1, replace = T))
  
  RPre_Isurv_P2C9  = (sample(Pre_Isurv_P2C9  , 1, replace = T))
  RPre_Isurv_C2P9  = (sample(Pre_Isurv_C2P9  , 1, replace = T))
  RPre_Isurv_P2C14 = (sample(Pre_Isurv_P2C14 , 1, replace = T))
  RPre_Isurv_C2P14 = (sample(Pre_Isurv_C2P14 , 1, replace = T))
  RPre_Isurv_P2C19 = (sample(Pre_Isurv_P2C19 , 1, replace = T))
  RPre_Isurv_C2P19 = (sample(Pre_Isurv_C2P19 , 1, replace = T))
  RPre_Isurv_P2C24 = (sample(Pre_Isurv_P2C24 , 1, replace = T))
  RPre_Isurv_C2P24 = (sample(Pre_Isurv_C2P24 , 1, replace = T))
  
  APre_Isurv_P2C9  = (1-RPre_Igerm_P2C9 )*RPre_Isurv_P2C9 
  APre_Isurv_C2P9  = (1-RPre_Igerm_C2P9 )*RPre_Isurv_C2P9 
  APre_Isurv_P2C14 = (1-RPre_Igerm_P2C14)*RPre_Isurv_P2C14
  APre_Isurv_C2P14 = (1-RPre_Igerm_C2P14)*RPre_Isurv_C2P14
  APre_Isurv_P2C19 = (1-RPre_Igerm_P2C19)*RPre_Isurv_P2C19
  APre_Isurv_C2P19 = (1-RPre_Igerm_C2P19)*RPre_Isurv_C2P19
  APre_Isurv_P2C24 = (1-RPre_Igerm_P2C24)*RPre_Isurv_P2C24
  APre_Isurv_C2P24 = (1-RPre_Igerm_C2P24)*RPre_Isurv_C2P24
  
  RPre_IFT_9   = (sample(Pre_IFT_9  , 1, replace=T))
  RPre_IFT_14  = (sample(Pre_IFT_14 , 1, replace=T))
  RPre_IFT_19  = (sample(Pre_IFT_19 , 1, replace=T))
  RPre_IFT_24  = (sample(Pre_IFT_24 , 1, replace=T))
  
  APre_IFT_P2C9 = (1-(APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPre_IFT_9
  APre_IFT_C2P9 = (1-(APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPre_IFT_9
  APre_IFT_P2C14= (1-(APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPre_IFT_14
  APre_IFT_C2P14= (1-(APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPre_IFT_14
  APre_IFT_P2C19= (1-(APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPre_IFT_19
  APre_IFT_C2P19= (1-(APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPre_IFT_19
  APre_IFT_P2C24= (1-(APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPre_IFT_24
  APre_IFT_C2P24= (1-(APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPre_IFT_24
  
  RPre_Ipv_P2C9 = (sample(Pre_Ipv_P2C9 ,1, replace= T))
  RPre_Ipv_C2P9 = (sample(Pre_Ipv_C2P9 ,1, replace= T))
  RPre_Ipv_P2C14= (sample(Pre_Ipv_P2C14,1, replace= T))
  RPre_Ipv_C2P14= (sample(Pre_Ipv_C2P14,1, replace= T))
  RPre_Ipv_P2C19= (sample(Pre_Ipv_P2C19,1, replace= T))
  RPre_Ipv_C2P19= (sample(Pre_Ipv_C2P19,1, replace= T))
  RPre_Ipv_P2C24= (sample(Pre_Ipv_P2C24,1, replace= T))
  RPre_Ipv_C2P24= (sample(Pre_Ipv_C2P24,1, replace= T))
  
  APre_Ipv_P2C9 =  (1-(APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPre_Ipv_P2C9 
  APre_Ipv_C2P9 =  (1-(APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPre_Ipv_C2P9 
  APre_Ipv_P2C14=  (1-(APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPre_Ipv_P2C14
  APre_Ipv_C2P14=  (1-(APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPre_Ipv_C2P14
  APre_Ipv_P2C19=  (1-(APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPre_Ipv_P2C19
  APre_Ipv_C2P19=  (1-(APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPre_Ipv_C2P19
  APre_Ipv_P2C24=  (1-(APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPre_Ipv_P2C24
  APre_Ipv_C2P24=  (1-(APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPre_Ipv_C2P24
  
  RPre_Iss_P2C9 = (sample(Pre_Iss_P2C9 ,1, replace= T)) 
  RPre_Iss_C2P9 = (sample(Pre_Iss_C2P9 ,1, replace= T)) 
  RPre_Iss_P2C14= (sample(Pre_Iss_P2C14,1, replace= T)) 
  RPre_Iss_C2P14= (sample(Pre_Iss_C2P14,1, replace= T)) 
  RPre_Iss_P2C19= (sample(Pre_Iss_P2C19,1, replace= T)) 
  RPre_Iss_C2P19= (sample(Pre_Iss_C2P19,1, replace= T)) 
  RPre_Iss_P2C24= (sample(Pre_Iss_P2C24,1, replace= T))
  RPre_Iss_C2P24= (sample(Pre_Iss_C2P24,1, replace= T))
  
  APre_Iss_P2C9 = (1-(APre_Ipv_P2C9 +APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPre_Iss_P2C9 
  APre_Iss_C2P9 = (1-(APre_Ipv_C2P9 +APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPre_Iss_C2P9 
  APre_Iss_P2C14= (1-(APre_Ipv_P2C14+APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPre_Iss_P2C14
  APre_Iss_C2P14= (1-(APre_Ipv_C2P14+APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPre_Iss_C2P14
  APre_Iss_P2C19= (1-(APre_Ipv_P2C19+APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPre_Iss_P2C19
  APre_Iss_C2P19= (1-(APre_Ipv_C2P19+APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPre_Iss_C2P19
  APre_Iss_P2C24= (1-(APre_Ipv_P2C24+APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPre_Iss_P2C24
  APre_Iss_C2P24= (1-(APre_Ipv_C2P24+APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPre_Iss_C2P24
  
  RPos_Igerm_P2C9 = (sample(Pos_Igerm_P2C9 , 1, replace = T))
  RPos_Igerm_C2P9 = (sample(Pos_Igerm_C2P9 , 1, replace = T))
  RPos_Igerm_P2C14= (sample(Pos_Igerm_P2C14, 1, replace = T))
  RPos_Igerm_C2P14= (sample(Pos_Igerm_C2P14, 1, replace = T))
  RPos_Igerm_P2C19= (sample(Pos_Igerm_P2C19, 1, replace = T))
  RPos_Igerm_C2P19= (sample(Pos_Igerm_C2P19, 1, replace = T))
  RPos_Igerm_P2C24= (sample(Pos_Igerm_P2C24, 1, replace = T))
  RPos_Igerm_C2P24= (sample(Pos_Igerm_C2P24, 1, replace = T))
  
  APos_Igerm_P2C9 = (1-(APre_Iss_P2C9 +APre_Ipv_P2C9 +APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Igerm_P2C9 
  APos_Igerm_C2P9 = (1-(APre_Iss_C2P9 +APre_Ipv_C2P9 +APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Igerm_C2P9 
  APos_Igerm_P2C14= (1-(APre_Iss_P2C14+APre_Ipv_P2C14+APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Igerm_P2C14
  APos_Igerm_C2P14= (1-(APre_Iss_C2P14+APre_Ipv_C2P14+APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Igerm_C2P14
  APos_Igerm_P2C19= (1-(APre_Iss_P2C19+APre_Ipv_P2C19+APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Igerm_P2C19
  APos_Igerm_C2P19= (1-(APre_Iss_C2P19+APre_Ipv_C2P19+APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Igerm_C2P19
  APos_Igerm_P2C24= (1-(APre_Iss_P2C24+APre_Ipv_P2C24+APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Igerm_P2C24
  APos_Igerm_C2P24= (1-(APre_Iss_C2P24+APre_Ipv_C2P24+APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Igerm_C2P24
  
  RPos_Isurv_P2C9  = (sample(Pos_Isurv_P2C9  , 1, replace = T))
  RPos_Isurv_C2P9  = (sample(Pos_Isurv_C2P9  , 1, replace = T))
  RPos_Isurv_P2C14 = (sample(Pos_Isurv_P2C14 , 1, replace = T))
  RPos_Isurv_C2P14 = (sample(Pos_Isurv_C2P14 , 1, replace = T))
  RPos_Isurv_P2C19 = (sample(Pos_Isurv_P2C19 , 1, replace = T))
  RPos_Isurv_C2P19 = (sample(Pos_Isurv_C2P19 , 1, replace = T))
  RPos_Isurv_P2C24 = (sample(Pos_Isurv_P2C24 , 1, replace = T))
  RPos_Isurv_C2P24 = (sample(Pos_Isurv_C2P24 , 1, replace = T))
  
  APos_Isurv_P2C9  =  (1-(APos_Igerm_P2C9 +APre_Iss_P2C9 +APre_Ipv_P2C9 +APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Isurv_P2C9 
  APos_Isurv_C2P9  =  (1-(APos_Igerm_C2P9 +APre_Iss_C2P9 +APre_Ipv_C2P9 +APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Isurv_C2P9 
  APos_Isurv_P2C14 =  (1-(APos_Igerm_P2C14+APre_Iss_P2C14+APre_Ipv_P2C14+APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Isurv_P2C14
  APos_Isurv_C2P14 =  (1-(APos_Igerm_C2P14+APre_Iss_C2P14+APre_Ipv_C2P14+APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Isurv_C2P14
  APos_Isurv_P2C19 =  (1-(APos_Igerm_P2C19+APre_Iss_P2C19+APre_Ipv_P2C19+APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Isurv_P2C19
  APos_Isurv_C2P19 =  (1-(APos_Igerm_C2P19+APre_Iss_C2P19+APre_Ipv_C2P19+APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Isurv_C2P19
  APos_Isurv_P2C24 =  (1-(APos_Igerm_P2C24+APre_Iss_P2C24+APre_Ipv_P2C24+APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Isurv_P2C24
  APos_Isurv_C2P24 =  (1-(APos_Igerm_C2P24+APre_Iss_C2P24+APre_Ipv_C2P24+APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Isurv_C2P24
  
  RPos_IFT_P2C9 = (sample(Pos_IFT_P2C9 , 1, replace=T))
  RPos_IFT_C2P9 = (sample(Pos_IFT_C2P9 , 1, replace=T))
  RPos_IFT_P2C14= (sample(Pos_IFT_P2C14, 1, replace=T))
  RPos_IFT_C2P14= (sample(Pos_IFT_C2P14, 1, replace=T))
  RPos_IFT_P2C19= (sample(Pos_IFT_P2C19, 1, replace=T))
  RPos_IFT_C2P19= (sample(Pos_IFT_C2P19, 1, replace=T))
  RPos_IFT_P2C24= (sample(Pos_IFT_P2C24, 1, replace=T))
  RPos_IFT_C2P24= (sample(Pos_IFT_C2P24, 1, replace=T))
  
  APos_IFT_P2C9 = (1-(APos_Isurv_P2C9 +APos_Igerm_P2C9 +APre_Iss_P2C9 +APre_Ipv_P2C9 +APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_IFT_P2C9 
  APos_IFT_C2P9 = (1-(APos_Isurv_C2P9 +APos_Igerm_C2P9 +APre_Iss_C2P9 +APre_Ipv_C2P9 +APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_IFT_C2P9 
  APos_IFT_P2C14= (1-(APos_Isurv_P2C14+APos_Igerm_P2C14+APre_Iss_P2C14+APre_Ipv_P2C14+APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_IFT_P2C14
  APos_IFT_C2P14= (1-(APos_Isurv_C2P14+APos_Igerm_C2P14+APre_Iss_C2P14+APre_Ipv_C2P14+APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_IFT_C2P14
  APos_IFT_P2C19= (1-(APos_Isurv_P2C19+APos_Igerm_P2C19+APre_Iss_P2C19+APre_Ipv_P2C19+APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_IFT_P2C19
  APos_IFT_C2P19= (1-(APos_Isurv_C2P19+APos_Igerm_C2P19+APre_Iss_C2P19+APre_Ipv_C2P19+APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_IFT_C2P19
  APos_IFT_P2C24= (1-(APos_Isurv_P2C24+APos_Igerm_P2C24+APre_Iss_P2C24+APre_Ipv_P2C24+APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_IFT_P2C24
  APos_IFT_C2P24= (1-(APos_Isurv_C2P24+APos_Igerm_C2P24+APre_Iss_C2P24+APre_Ipv_C2P24+APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_IFT_C2P24
  
  RPos_Ipv_P2C9 = (sample(Pos_Ipv_P2C9 ,1, replace= T))
  RPos_Ipv_C2P9 = (sample(Pos_Ipv_C2P9 ,1, replace= T))
  RPos_Ipv_P2C14= (sample(Pos_Ipv_P2C14,1, replace= T))
  RPos_Ipv_C2P14= (sample(Pos_Ipv_C2P14,1, replace= T))
  RPos_Ipv_P2C19= (sample(Pos_Ipv_P2C19,1, replace= T))
  RPos_Ipv_C2P19= (sample(Pos_Ipv_C2P19,1, replace= T))
  RPos_Ipv_P2C24= (sample(Pos_Ipv_P2C24,1, replace= T))
  RPos_Ipv_C2P24= (sample(Pos_Ipv_C2P24,1, replace= T))
  
  APos_Ipv_P2C9 =  (1-(APos_IFT_P2C9 +APos_Isurv_P2C9 +APos_Igerm_P2C9 +APre_Iss_P2C9 +APre_Ipv_P2C9 +APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Ipv_P2C9 
  APos_Ipv_C2P9 =  (1-(APos_IFT_C2P9 +APos_Isurv_C2P9 +APos_Igerm_C2P9 +APre_Iss_C2P9 +APre_Ipv_C2P9 +APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Ipv_C2P9 
  APos_Ipv_P2C14=  (1-(APos_IFT_P2C14+APos_Isurv_P2C14+APos_Igerm_P2C14+APre_Iss_P2C14+APre_Ipv_P2C14+APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Ipv_P2C14
  APos_Ipv_C2P14=  (1-(APos_IFT_C2P14+APos_Isurv_C2P14+APos_Igerm_C2P14+APre_Iss_C2P14+APre_Ipv_C2P14+APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Ipv_C2P14
  APos_Ipv_P2C19=  (1-(APos_IFT_P2C19+APos_Isurv_P2C19+APos_Igerm_P2C19+APre_Iss_P2C19+APre_Ipv_P2C19+APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Ipv_P2C19
  APos_Ipv_C2P19=  (1-(APos_IFT_C2P19+APos_Isurv_C2P19+APos_Igerm_C2P19+APre_Iss_C2P19+APre_Ipv_C2P19+APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Ipv_C2P19
  APos_Ipv_P2C24=  (1-(APos_IFT_P2C24+APos_Isurv_P2C24+APos_Igerm_P2C24+APre_Iss_P2C24+APre_Ipv_P2C24+APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Ipv_P2C24
  APos_Ipv_C2P24=  (1-(APos_IFT_C2P24+APos_Isurv_C2P24+APos_Igerm_C2P24+APre_Iss_C2P24+APre_Ipv_C2P24+APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Ipv_C2P24
  
  RPos_Iss_P2C9 = (sample(Pos_Iss_P2C9 ,1, replace= T)) 
  RPos_Iss_C2P9 = (sample(Pos_Iss_C2P9 ,1, replace= T)) 
  RPos_Iss_P2C14= (sample(Pos_Iss_P2C14,1, replace= T)) 
  RPos_Iss_C2P14= (sample(Pos_Iss_C2P14,1, replace= T)) 
  RPos_Iss_P2C19= (sample(Pos_Iss_P2C19,1, replace= T)) 
  RPos_Iss_C2P19= (sample(Pos_Iss_C2P19,1, replace= T)) 
  RPos_Iss_P2C24= (sample(Pos_Iss_P2C24,1, replace= T))
  RPos_Iss_C2P24= (sample(Pos_Iss_C2P24,1, replace= T))
  
  APos_Iss_P2C9 = (1-(APos_Ipv_P2C9 +APos_IFT_P2C9 +APos_Isurv_P2C9 +APos_Igerm_P2C9 +APre_Iss_P2C9 +APre_Ipv_P2C9 +APre_IFT_P2C9 +APre_Isurv_P2C9 +RPre_Igerm_P2C9 ))*RPos_Iss_P2C9 
  APos_Iss_C2P9 = (1-(APos_Ipv_C2P9 +APos_IFT_C2P9 +APos_Isurv_C2P9 +APos_Igerm_C2P9 +APre_Iss_C2P9 +APre_Ipv_C2P9 +APre_IFT_C2P9 +APre_Isurv_C2P9 +RPre_Igerm_C2P9 ))*RPos_Iss_C2P9 
  APos_Iss_P2C14= (1-(APos_Ipv_P2C14+APos_IFT_P2C14+APos_Isurv_P2C14+APos_Igerm_P2C14+APre_Iss_P2C14+APre_Ipv_P2C14+APre_IFT_P2C14+APre_Isurv_P2C14+RPre_Igerm_P2C14))*RPos_Iss_P2C14
  APos_Iss_C2P14= (1-(APos_Ipv_C2P14+APos_IFT_C2P14+APos_Isurv_C2P14+APos_Igerm_C2P14+APre_Iss_C2P14+APre_Ipv_C2P14+APre_IFT_C2P14+APre_Isurv_C2P14+RPre_Igerm_C2P14))*RPos_Iss_C2P14
  APos_Iss_P2C19= (1-(APos_Ipv_P2C19+APos_IFT_P2C19+APos_Isurv_P2C19+APos_Igerm_P2C19+APre_Iss_P2C19+APre_Ipv_P2C19+APre_IFT_P2C19+APre_Isurv_P2C19+RPre_Igerm_P2C19))*RPos_Iss_P2C19
  APos_Iss_C2P19= (1-(APos_Ipv_C2P19+APos_IFT_C2P19+APos_Isurv_C2P19+APos_Igerm_C2P19+APre_Iss_C2P19+APre_Ipv_C2P19+APre_IFT_C2P19+APre_Isurv_C2P19+RPre_Igerm_C2P19))*RPos_Iss_C2P19
  APos_Iss_P2C24= (1-(APos_Ipv_P2C24+APos_IFT_P2C24+APos_Isurv_P2C24+APos_Igerm_P2C24+APre_Iss_P2C24+APre_Ipv_P2C24+APre_IFT_P2C24+APre_Isurv_P2C24+RPre_Igerm_P2C24))*RPos_Iss_P2C24
  APos_Iss_C2P24= (1-(APos_Ipv_C2P24+APos_IFT_C2P24+APos_Isurv_C2P24+APos_Igerm_C2P24+APre_Iss_C2P24+APre_Ipv_C2P24+APre_IFT_C2P24+APre_Isurv_C2P24+RPre_Igerm_C2P24))*RPos_Iss_C2P24
  
  aItotal_P2C9 = c(aItotal_P2C9 ,(RPre_Igerm_P2C9 +APre_Isurv_P2C9 +APre_IFT_P2C9 +APre_Ipv_P2C9 +APre_Iss_P2C9 +APos_Igerm_P2C9 +APos_Isurv_P2C9 +APos_IFT_P2C9 +APos_Ipv_P2C9 +APos_Iss_P2C9 ))
  aItotal_C2P9 = c(aItotal_C2P9 ,(RPre_Igerm_C2P9 +APre_Isurv_C2P9 +APre_IFT_C2P9 +APre_Ipv_C2P9 +APre_Iss_C2P9 +APos_Igerm_C2P9 +APos_Isurv_C2P9 +APos_IFT_C2P9 +APos_Ipv_C2P9 +APos_Iss_C2P9 ))
  aItotal_P2C14= c(aItotal_P2C14,(RPre_Igerm_P2C14+APre_Isurv_P2C14+APre_IFT_P2C14+APre_Ipv_P2C14+APre_Iss_P2C14+APos_Igerm_P2C14+APos_Isurv_P2C14+APos_IFT_P2C14+APos_Ipv_P2C14+APos_Iss_P2C14))
  aItotal_C2P14= c(aItotal_C2P14,(RPre_Igerm_C2P14+APre_Isurv_C2P14+APre_IFT_C2P14+APre_Ipv_C2P14+APre_Iss_C2P14+APos_Igerm_C2P14+APos_Isurv_C2P14+APos_IFT_C2P14+APos_Ipv_C2P14+APos_Iss_C2P14))
  aItotal_P2C19= c(aItotal_P2C19,(RPre_Igerm_P2C19+APre_Isurv_P2C19+APre_IFT_P2C19+APre_Ipv_P2C19+APre_Iss_P2C19+APos_Igerm_P2C19+APos_Isurv_P2C19+APos_IFT_P2C19+APos_Ipv_P2C19+APos_Iss_P2C19))
  aItotal_C2P19= c(aItotal_C2P19,(RPre_Igerm_C2P19+APre_Isurv_C2P19+APre_IFT_C2P19+APre_Ipv_C2P19+APre_Iss_C2P19+APos_Igerm_C2P19+APos_Isurv_C2P19+APos_IFT_C2P19+APos_Ipv_C2P19+APos_Iss_C2P19))
  aItotal_P2C24= c(aItotal_P2C24,(RPre_Igerm_P2C24+APre_Isurv_P2C24+APre_IFT_P2C24+APre_Ipv_P2C24+APre_Iss_P2C24+APos_Igerm_P2C24+APos_Isurv_P2C24+APos_IFT_P2C24+APos_Ipv_P2C24+APos_Iss_P2C24))
  aItotal_C2P24= c(aItotal_C2P24,(RPre_Igerm_C2P24+APre_Isurv_C2P24+APre_IFT_C2P24+APre_Ipv_C2P24+APre_Iss_C2P24+APos_Igerm_C2P24+APos_Isurv_C2P24+APos_IFT_C2P24+APos_Ipv_C2P24+APos_Iss_C2P24))
  
}

summary(sort(aItotal_P2C9 )[c(250, 9750)])
summary(sort(aItotal_C2P9 )[c(250, 9750)])
summary(sort(aItotal_P2C14)[c(250, 9750)])
summary(sort(aItotal_C2P14)[c(250, 9750)])
summary(sort(aItotal_P2C19)[c(250, 9750)])
summary(sort(aItotal_C2P19)[c(250, 9750)])
summary(sort(aItotal_P2C24)[c(250, 9750)])
summary(sort(aItotal_C2P24)[c(250, 9750)])


sd(aItotal_P2C9 )
sd(aItotal_C2P9 )
sd(aItotal_P2C14)
sd(aItotal_C2P14)
sd(aItotal_P2C19)
sd(aItotal_C2P19)
sd(aItotal_P2C24)
sd(aItotal_C2P24)
