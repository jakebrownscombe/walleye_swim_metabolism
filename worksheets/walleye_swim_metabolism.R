#Walleye swim tunnel respirometry data analysis
#Jake Brownscombe
#5-Feb-2024

#required packages ----
library(dplyr)
library(descr)
library(ggplot2)
library(tidyverse)
library(nlme)
library(lmerTest)
library(lme4)
library(Matrix)
library(AICcmodavg)
library(patchwork) 
source("~/github/walleye_swim_metabolism/worksheets/functions.R") #R2 calc for lme


#load data ----

swim.data <- readRDS("~/github/walleye_swim_metabolism/data/swim.data.rds")
swim.meta <- readRDS("~/github/walleye_swim_metabolism/data/swim.meta.rds")
swim.data.unfiltered <- readRDS("~/github/walleye_swim_metabolism/data/swim.data.unfiltered.rds")





#data removed ----

#This initial set of trials at 21C, fish showed clear signs of stress/exhaustion. This is illustrated in the paper appendix:
ggplot(swim.data.unfiltered %>% filter(!is.na(group)), aes(speedcm.s, massMO2, col=temp.cat, linetype=group))+geom_point()+geom_smooth(method="lm")+
  scale_color_viridis_d(option="plasma", begin=0, end=0.9, name="ºC")+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))+theme_bw()
#     




#data summary ----
#numbers, sizes at each temp
head(swim.data)
swim.data %>% group_by(temp.cat2) %>% summarise(n=length(unique(FishID)), mean(weight.kg), sd(weight.kg), min(weight.kg),max(weight.kg),
                                                mean(length.mm), sd(length.mm), min(length.mm),max(length.mm) ) %>% as.data.frame()
#








#acceleration - metabolic rate ----
swim.data.acc <- swim.data %>% filter(!is.na(gs))
head(swim.data.acc)

ggplot(swim.data.acc, aes(gs, log_massMO2, col=temp.cat))+geom_point()+geom_smooth(method="lm")+theme_bw()+
  theme_bw()+xlab("Acceleration (g)")+ ylab(bquote(~MO[2]* ' (mg'* ~O[2] ~hour^-1 ~kg^-1*')'))+
  ggplot(swim.data.acc, aes(log_weight.kg, log_massMO2))+geom_point()+geom_smooth(method="lm")+theme_bw()


#model 
M1<-lme(log_massMO2 ~ gs + temp + log_weight.kg + sex +
          gs*temp + gs*log_weight.kg + log_weight.kg*sex + gs*sex,
        random = ~1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc) #random intercept
M2<-lme(log_massMO2 ~ gs + temp + log_weight.kg + sex +
          gs*temp + gs*log_weight.kg + log_weight.kg*sex + gs*sex,
        random = ~1 + gs | FishID, method = "REML", na.action = na.omit, data = swim.data.acc) #random slope+int
M3<-lme(log_massMO2 ~ gs + temp + log_weight.kg + sex +
          gs*temp + gs*log_weight.kg + log_weight.kg*sex + gs*sex,
        random = ~gs - 1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc) #random slope     

anova(M1,M2,M3)

rsquared.lme(list(M1)) 
rsquared.lme(list(M2)) 
rsquared.lme(list(M3)) 

#intercepts the best. 
MO2.acc.model <-lme(log_massMO2 ~ gs + temp + log_weight.kg + sex +
                      gs*temp + gs*log_weight.kg + log_weight.kg*sex + gs*sex,
                    random = ~1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc)
drop1(MO2.acc.model, test="Chisq")
#drop interactions

MO2.acc.model2 <-lme(log_massMO2 ~ gs + temp + log_weight.kg + sex,
                     random = ~1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc)
drop1(MO2.acc.model2, test="Chisq")
#drop sex

MO2.acc.model3 <-lme(log_massMO2 ~ gs + temp + log_weight.kg,
                     random = ~1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc)
drop1(MO2.acc.model3, test="Chisq")

#FINAL model 
MO2.acc.model.final <-lme(log_massMO2 ~ gs + temp + log_weight.kg,
                          random = ~1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc)
summary(MO2.acc.model.final)
rsquared.lme(list(MO2.acc.model.final))


#model val
E <- resid(MO2.acc.model.final)
F <- fitted(MO2.acc.model.final)
plot(E~F)
plot(E ~ swim.data.acc$gs)
plot(E ~ swim.data.acc$temp)
plot(E ~ swim.data.acc$log_weight.kg)
plot(F ~ swim.data.acc$log_massMO2)
#decent model. 


#predictions
mass <- expand.grid(weight.kg= seq(0.3, 3.5, 0.1))
temps <- expand.grid(temp= seq(1, 30, 1))
acc <- expand.grid(gs= seq(0.0, 2, 0.1))

preds.acc <- merge(mass, temps)
preds.acc <- merge(preds.acc, acc)
preds.acc$log_weight.kg <- log(preds.acc$weight.kg)

preds.acc$MO2.pred <- exp(predict(MO2.acc.model.final, preds.acc, level = 0))
head(preds.acc)

sum <- swim.data.acc %>% group_by(temp.cat) %>% summarise(temp.cat2=round(mean(temp),0))
swim.data.acc$temp.cat2 <- sum$temp.cat2[match(swim.data.acc$temp.cat, sum$temp.cat)]
head(swim.data.acc)


#plots. 
preds.acc2 <- preds.acc %>% filter(temp==5 & weight.kg==0.3 | temp==5 & weight.kg==3.5 |
                                     temp==11 & weight.kg==0.3 | temp==11 & weight.kg==3.5 |
                                     temp==16 & weight.kg==0.3 | temp==16 & weight.kg==3.5 |
                                     temp==21 & weight.kg==0.3 | temp==21 & weight.kg==3.5 )

ggplot(swim.data.acc, aes(gs, massMO2, col=as.factor(temp.cat2), pch=as.factor(temp.cat2)))+geom_point()+
  geom_smooth(data=preds.acc2, aes(gs, MO2.pred, col=as.factor(temp), linetype=as.factor(weight.kg)), fill="NA", inherit.aes = F)+
  scale_y_continuous(limits=c(0,500))+scale_x_continuous(limits=c(0,2.0))+
  scale_linetype(name="Mass (kg)")+
  scale_color_viridis_d(option="plasma", begin=0, end=0.9, name="Temperature (ºC)")+
  scale_shape(name="Temperature (ºC)")+
  theme_bw()+
  labs(x="Acceleration (g)", y=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2),
                                                ' (mg O'['2'], ' kg'^'-1', ' hour'^'-1', ')')), 
       title="Acceleration - Metabolism", subtitle="A) Acceleration + Temperature")+
  
  ggplot(preds.acc %>% filter(temp==15), aes(gs, weight.kg,  fill=MO2.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", 
                       name=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2))))+
  theme_classic()+
  labs(x="Acceleration (g)", y="Mass (kg)", subtitle="B) Acceleration + Mass (15ºC)")+
  
  ggplot(preds.acc %>% filter(weight.kg==1), aes(gs, temp, fill=MO2.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", 
                       name=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2))))+
  theme_classic()+
  labs(x="Acceleration (g)", y="Temperature (ºC)", subtitle="C) Acceleration + Temperature (1 kg fish)")+
  
  ggplot(preds.acc %>% filter(gs==1), aes(weight.kg, temp, fill=MO2.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", 
                       name=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2))))+
  theme_classic()+
  labs(x="Mass (kg)", y="Temperature (ºC)", subtitle="D) Mass + Temperature (1 g acceleration)")

#











#speed - metabolic rate  ----
head(swim.data)
swim.data$temp.cat2 <- sum$temp.cat2[match(swim.data$temp.cat, sum$temp.cat)]
ggplot(swim.data, aes(speed.corr, log_massMO2, col=as.factor(temp.cat2)))+geom_point()+geom_smooth(method="lm")


#MO2 modeling 
M1<-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + sex +
          speed.corr*temp + speed.corr*log_weight.kg + log_weight.kg*sex + speed.corr*sex,
        random = ~1 | FishID, method = "REML", na.action = na.omit, data = swim.data) #random intercept
M2<-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + sex +
          speed.corr*temp + speed.corr*log_weight.kg + log_weight.kg*sex + speed.corr*sex,
        random = ~1 + speed.corr | FishID, method = "REML", na.action = na.omit, data = swim.data) #random slope+int
M3<-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + sex +
          speed.corr*temp + speed.corr*log_weight.kg + log_weight.kg*sex + speed.corr*sex,
        random = ~speed.corr - 1 | FishID, method = "REML", na.action = na.omit, data = swim.data) #random slope     

anova(M1,M2,M3)

rsquared.lme(list(M1)) 
rsquared.lme(list(M2)) 
rsquared.lme(list(M3)) 
#chose both int + slope

M3.1 <-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + sex +
             speed.corr*temp + speed.corr*log_weight.kg + log_weight.kg*sex + speed.corr*sex,
           random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data) #random slope     
drop1(M3.1, test="Chisq")
#retain speed x weight 

M3.2 <-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + sex +
             speed.corr*log_weight.kg,
           random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data) #random slope     
drop1(M3.2, test="Chisq")
#drop sex

M3.3 <-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + speed.corr*log_weight.kg,
           random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data) #random slope     
drop1(M3.3, test="Chisq")


#final model
swim.final <-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg,
                 random = ~speed.corr - 1 | FishID, method = "REML", na.action = na.omit, data = swim.data) #random slope   
summary(swim.final)
rsquared.lme(list(swim.final)) 

#model val
E <- resid(swim.final)
F <- fitted(swim.final)
plot(E~F)
plot(E ~ swim.data$speed.corr)
plot(E ~ swim.data$temp)
plot(E ~ swim.data$log_weight.kg)
plot(F ~ swim.data$log_massMO2)
#decent model. 


#predictions
mass <- expand.grid(weight.kg= seq(0.3, 3.5, 0.1))
temps <- expand.grid(temp= seq(1, 30, 1))
speed <- expand.grid(speed.corr= seq(0.0, 100, 5))

preds.swim <- merge(mass, temps)
preds.swim <- merge(preds.swim, speed)
preds.swim$log_weight.kg <- log(preds.swim$weight.kg)
head(preds.swim)

preds.swim$MO2.pred <- exp(predict(swim.final, preds.swim, level = 0))
head(preds.swim)

preds.swim2 <- preds.swim %>% filter(temp==5 & weight.kg==0.3 | temp==5 & weight.kg==3.5 |
                                       temp==11 & weight.kg==0.3 | temp==11 & weight.kg==3.5 |
                                       temp==16 & weight.kg==0.3 | temp==16 & weight.kg==3.5 |
                                       temp==21 & weight.kg==0.3 | temp==21 & weight.kg==3.5 )


ggplot(swim.data, aes(speed.corr, massMO2, col=as.factor(temp.cat2), pch=as.factor(temp.cat2)))+geom_point()+
  geom_smooth(data=preds.swim2, aes(speed.corr, MO2.pred, col=as.factor(temp), linetype=as.factor(weight.kg)), fill=NA, inherit.aes = F)+
  scale_y_continuous(limits=c(0,650))+scale_x_continuous(limits=c(0,100))+scale_linetype(name="Mass (kg)")+
  scale_color_viridis_d(option="plasma", begin=0, end=0.9, name="Temperature (ºC)")+theme_bw()+
  scale_shape(name="Temperature (ºC)")+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2),
                                                                  ' (mg O'['2'], ' kg'^'-1', ' hour'^'-1', ')')), 
       title="Swimming Speed - Metabolism", subtitle="A) Speed + Temperature")+
  
  ggplot(preds.swim %>% filter(temp==15), aes(speed.corr, weight.kg,  fill=MO2.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", 
                       name=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2))))+
  theme_classic()+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y="Mass (kg)", subtitle="B) Speed + Mass (15ºC)")+
  
  ggplot(preds.swim %>% filter(weight.kg==1), aes(speed.corr, temp, fill=MO2.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", 
                       name=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2))))+
  theme_classic()+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y="Temperature (ºC)", subtitle="C) Speed + Temperature (1 kg fish)")+
  
  ggplot(preds.swim %>% filter(speed.corr==50), aes(weight.kg, temp, fill=MO2.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", 
                       name=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2))))+
  theme_classic()+
  labs(x="Mass (kg)", y="Temperature (ºC)", 
       subtitle=expression('D) Mass + Temperature (50'~ cm~s^-1~')'))
#













#Tagging effects ----
#focus on hatchery fish, which have untagged individuals to compare
head(swim.data)
tag.data <- swim.data %>% filter(group=="hatchery control" | group=="hatchery tag")

#model 
Mt1 <-lme(log_massMO2 ~ speed.corr + temp + log_weight.kg + tagged,
          random = ~speed.corr - 1 | FishID, method = "REML", na.action = na.omit, data = tag.data) #random slope   
summary(Mt1) 
#no clear tagging effect

#plot
ggplot(tag.data, aes(speed.corr, massMO2, col=tagged))+geom_point()+geom_smooth(method="lm")+
  theme_bw()+labs(x=expression(Swimming~Speed~(cm~s^-1)), y=bquote(~MO[2]* ' (mg'* ~O[2] ~hour^-1 ~kg^-1*')'))+
  scale_color_discrete(name="Tagged")+scale_y_continuous(limits=c(0,500))+scale_x_continuous(limits=c(35,90))
#









#Hatchery effects ----
swim.data$group2 <- ifelse(swim.data$group=="hatchery control" | swim.data$group=="hatchery tag", "hatchery", "wild")
ggplot(swim.data, aes(group2, weight.kg))+geom_boxplot()
ggplot(swim.data %>% filter(weight.kg<1.8), aes(group2, weight.kg))+geom_boxplot()
ggplot(swim.data %>% filter(weight.kg<1.8), aes(temp, massMO2, col=group2))+geom_point()
#don't have data for this question due to lack of size/temp overlap 













#acceleration - swimming speed ----
swim.data.acc$loggs <- log(swim.data.acc$gs)
ggplot(swim.data.acc, aes(speed.corr, loggs, col=temp.cat))+geom_point()+geom_smooth(method="lm")+theme_bw()


#figure out a general intercept to know when fish are at 'rest'
SM.base <-lme(loggs ~ speed.corr,
              random = ~speed.corr -1 | FishID, method = "REML", 
              data = swim.data.acc)# random slope
summary(SM.base)
exp(-1.1815350)
#0.31 is the predicted at rest value


#full model
SM1<-lme(loggs ~ speed.corr + log_weight.kg + temp + sex +
           speed.corr*log_weight.kg + speed.corr*temp + speed.corr*sex,
         random = ~1 | FishID, method = "REML", 
         data = swim.data.acc)# random intercept
SM2<-lme(loggs ~ speed.corr + log_weight.kg + temp + sex +
           speed.corr*log_weight.kg + speed.corr*temp + speed.corr*sex,
         random = ~speed.corr -1 | FishID, method = "REML", 
         data = swim.data.acc)# random slope
SM3<-lme(loggs ~ speed.corr + log_weight.kg + temp + sex +
           speed.corr*log_weight.kg + speed.corr*temp + speed.corr*sex,
         random = ~1 + speed.corr | FishID, method = "REML", 
         data = swim.data.acc)# both
anova(SM1,SM2, SM3)
rsquared.lme(list(SM1))
rsquared.lme(list(SM2))
rsquared.lme(list(SM3))
#go with slope 

swim.acc.model <-lme(loggs ~ speed.corr + log_weight.kg + temp + sex +
                       speed.corr*log_weight.kg + speed.corr*temp + speed.corr*sex,
                     random = ~speed.corr -1 | FishID, method = "ML", 
                     data = swim.data.acc)
drop1(swim.acc.model, test="Chisq")
#drop interactions

swim.acc.model2 <-lme(loggs ~ speed.corr + log_weight.kg + temp + sex,
                      random = ~speed.corr -1 | FishID, method = "ML", 
                      data = swim.data.acc)
drop1(swim.acc.model2, test="Chisq")
#drop temp and sex

#final model
swim.acc.model.final <-lme(loggs ~ speed.corr + log_weight.kg,
                           random = ~speed.corr -1 | FishID, method = "REML", 
                           data = swim.data.acc)# both

summary(swim.acc.model.final)
rsquared.lme(list(swim.acc.model.final))
exp(-1.1500124)
#0.32 is the predicted at rest value. not too different from base model


#model validation
E <- resid(swim.acc.model.final)
F <- fitted(swim.acc.model.final)
plot(E~F)
plot(E ~ swim.data.acc$speed.corr)
plot(E ~ swim.data.acc$temp)
plot(E ~ swim.data.acc$log_weight.kg)
plot(F ~ swim.data.acc$loggs)
#not great but OK

#predictions
mass <- expand.grid(weight.kg= seq(0.3, 3.5, 0.1))
speed <- expand.grid(speed.corr=seq(0,100,1))
preds.acc.swim <- merge(mass, speed)
preds.acc.swim$log_weight.kg <- log(preds.acc.swim$weight.kg)
preds.acc.swim$gs.pred <- exp(predict(swim.acc.model.final, preds.acc.swim, level = 0))
head(preds.acc.swim)

#plot
ggplot(swim.data.acc, aes(speed.corr, gs, col=weight.kg))+geom_point()+
  geom_smooth(data=preds.acc.swim %>% filter(weight.kg==0.3 | weight.kg==1.0 | weight.kg==3.5),
              aes(speed.corr, gs.pred, group=as.factor(weight.kg)))+theme_bw()+
  scale_y_continuous(limits=c(0,3))+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y="Acceleration (g)")+scale_color_viridis_c(name="Mass (kg)")
#








#Ucrit ---- 
head(swim.data)
Ucrit <- swim.data %>% group_by(Trial, FishID, temp.cat, group) %>% summarise(Ucrit=mean(ucrit), temp=mean(temp), weight.kg=mean(weight.kg)) %>% 
  filter(!is.na(Ucrit))
Ucrit$log_weight <- log(Ucrit$weight.kg)

ggplot(Ucrit, aes(temp, Ucrit))+geom_point()+
  geom_smooth(method="lm", formula=y~poly(x,2))+
  theme_bw()+scale_y_continuous(limits=c(0,100))+
  labs(x="Temperature (ºC)", y="Critical Swimming Speed (cm/sec)")


#model
Ucrit1<-lme(Ucrit ~ poly(temp,2) + log_weight + poly(temp,2)*log_weight,
            random = ~1 | FishID, method = "REML", na.action = na.omit, data = Ucrit) #random intercept
Ucrit2<-lme(Ucrit ~ poly(temp,2) + log_weight + poly(temp,2)*log_weight,
            random = ~1 + temp | FishID, method = "REML", na.action = na.omit, data = Ucrit) #random slope+int
Ucrit3<-lme(Ucrit ~ poly(temp,2) + log_weight + poly(temp,2)*log_weight,
            random = ~temp - 1 | FishID, method = "REML", na.action = na.omit, data = Ucrit) #random slope     

anova(Ucrit1,Ucrit2,Ucrit3)

rsquared.lme(list(Ucrit1)) 
rsquared.lme(list(Ucrit2)) 
rsquared.lme(list(Ucrit3)) 
#hmm tough call. intercept  

Ucrit1.1<-lme(Ucrit ~ poly(temp,2) + log_weight + poly(temp,2)*log_weight,
              random = ~1 | FishID, method = "ML", na.action = na.omit, data = Ucrit) #random intercept
drop1(Ucrit1.1, test="Chisq")

Ucrit1.1<-lme(Ucrit ~ poly(temp,2) + log_weight,
              random = ~1 | FishID, method = "ML", na.action = na.omit, data = Ucrit) #random intercept
drop1(Ucrit1.1, test="Chisq")

Ucrit1.1<-lme(Ucrit ~ poly(temp,2),
              random = ~1 | FishID, method = "ML", na.action = na.omit, data = Ucrit) #random intercept
drop1(Ucrit1.1, test="Chisq")

Ucrit.final <-lme(Ucrit ~ poly(temp,2),
                  random = ~1 | FishID, method = "REML", na.action = na.omit, data = Ucrit) #random intercept
summary(Ucrit.final)
rsquared.lme(list(Ucrit.final)) 



#Ucrit tagging effects
head(Ucrit)
Ucrit$group <- swim.data$group[match(Ucrit$Trial, swim.data$Trial)]
Ucrit.tag <- Ucrit %>% filter(group=="hatchery control" | group=="hatchery tag")

ggplot(Ucrit.tag, aes(temp, Ucrit, col=group))+geom_point()
#seems alright. remove a couple of 13 temps. 

Ucrit.tag <- Ucrit.tag %>% filter(temp<12 | temp>15)
Ucrit.tag$temp.cat <- ifelse(Ucrit.tag$temp<12, 11, 16)
Ucrit.tag$tagged <- ifelse(Ucrit.tag$group=="hatchery tag", "yes", "no")
#plot

ggplot(Ucrit.tag, aes(as.factor(temp.cat), Ucrit, col=tagged))+geom_boxplot()+
  theme_bw()+scale_y_continuous(limits=c(0,80))+scale_color_discrete(name="Tagged")+
  labs(x="Temperature (ºC)", y=expression(Swimming~Speed~(cm~s^-1)))+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(fill = "white", colour = "black"))


#model
tag.ucrit.lme = lme(Ucrit ~ as.factor(temp.cat) + tagged + as.factor(temp.cat)*tagged, 
                    random = ~1 | FishID, method = "ML", na.action = na.omit, data=Ucrit.tag)
summary(tag.ucrit.lme)
drop1(tag.ucrit.lme)

tag.ucrit.lme = lme(Ucrit ~ as.factor(temp.cat) + tagged, 
                    random = ~1 | FishID, method = "REML", na.action = na.omit, data=Ucrit.tag)
summary(tag.ucrit.lme)



#tagging effects plots together
ggplot(tag.data, aes(speed.corr, massMO2, col=tagged))+geom_point()+geom_smooth(method="lm")+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), 
       y=expression(paste(italic(dot(M)), scriptstyle(O), scriptscriptstyle(2), 
                          ' (mg O'['2'], ' kg'^'-1', ' hour'^'-1', ')')), 
       title="Tagging Effects", subtitle="A) Metabolic Rate")+
  scale_color_discrete(name="Tagged")+scale_y_continuous(limits=c(0,500))+scale_x_continuous(limits=c(35,90))+theme_bw()+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(fill = "white", colour = "black"))+
  
  ggplot(Ucrit.tag, aes(as.factor(temp.cat), Ucrit, col=tagged))+geom_boxplot()+
  theme_bw()+scale_y_continuous(limits=c(0,80))+scale_color_discrete(name="Tagged")+
  labs(x="Temperature (ºC)", y=expression(Swimming~Speed~(cm~s^-1)), subtitle="B) Critical Swimming Speed")+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(fill = "white", colour = "black"))+
  
  plot_layout(ncol=1)
#










#acceleration - tailbeat frequency ----
head(swim.data.acc)
swim.data.acc$log_TB <- log(swim.data.acc$TB.s)
ggplot(swim.data.acc, aes(gs, log_TB, col=temp))+geom_point()+geom_smooth(method="lm")+
  theme_bw()+scale_y_continuous(limits=c(0,1.5))+scale_x_continuous(limits=c(0,2))


#model
TBF1<-lme(log_TB ~ gs + log_weight.kg + temp,
          random = ~1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc %>% filter(!is.na(TB.s))) #random intercept
TBF2<-lme(log_TB ~ gs + log_weight.kg + temp,
          random = ~1 + gs | FishID, method = "REML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope+int
TBF3<-lme(log_TB ~ gs + log_weight.kg + temp,
          random = ~gs - 1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     

anova(TBF1,TBF2,TBF3)
rsquared.lme(list(TBF1)) 
rsquared.lme(list(TBF2)) 
rsquared.lme(list(TBF3)) 

#intercept 
TBF3.1<-lme(log_TB ~ gs + log_weight.kg + temp,
            random = ~1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     
summary(TBF3.1)
rsquared.lme(list(TBF3.1)) 
drop1(TBF3.1, test="Chisq")

#final model:
summary(TBF1)
rsquared.lme(list(TBF1)) 


#preds
acc <- data.frame(gs=seq(0.3, 2, 0.1))
preds.TBF <- merge(acc, temps)
preds.TBF <- merge(preds.TBF, mass)
preds.TBF$log_weight.kg <- log(preds.TBF$weight.kg)
preds.TBF$TBF.pred <- exp(predict(TBF1, preds.TBF, level = 0))
head(preds.TBF)

preds.TBF2 <- preds.TBF %>% filter(temp==5 & weight.kg==0.3 | temp==5 & weight.kg==3.5 |
                                     temp==11 & weight.kg==0.3 | temp==11 & weight.kg==3.5 |
                                     temp==16 & weight.kg==0.3 | temp==16 & weight.kg==3.5 |
                                     temp==21 & weight.kg==0.3 | temp==21 & weight.kg==3.5)
head(preds.TBF2)

ggplot(swim.data.acc, aes(gs, TB.s, col=as.factor(temp.cat2)))+geom_point()+
  geom_line(data=preds.TBF2, aes(gs, TBF.pred, col=as.factor(temp), linetype=as.factor(weight.kg)))+
  theme_bw()+scale_y_continuous(limits=c(0,4.2))+scale_x_continuous(limits=c(0,2.2))+
  scale_color_viridis_d(option="plasma", end=0.9, name="ºC")+scale_linetype(name="Weight (kg)")+
  labs(x="Acceleration (g)", y=expression(Tailbeat~Frequency~(beats~sec^-1)), title="Tailbeat Frequency", subtitle="A) Acceleration")
#






#swim speed - tailbeat frequency ----
ggplot(swim.data.acc, aes(speed.corr, log_TB, col=temp))+geom_point()+geom_smooth(method="lm")
ggplot(swim.data.acc, aes(speed.corr, log_TB, col=weight.kg))+geom_point()+geom_smooth(method="lm")
ggplot(swim.data.acc, aes(weight.kg, temp))+geom_point()

TBFS1<-lme(log_TB ~ speed.corr + log_weight.kg + temp + speed.corr*temp + speed.corr*log_weight.kg,
           random = ~1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc %>% filter(!is.na(TB.s))) #random intercept
TBFS2<-lme(log_TB ~ speed.corr + log_weight.kg + temp + speed.corr*temp + speed.corr*log_weight.kg,
           random = ~1 + speed.corr| FishID, method = "REML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope+int
TBFS3<-lme(log_TB ~ speed.corr + log_weight.kg + temp + speed.corr*temp + speed.corr*log_weight.kg,
           random = ~speed.corr- 1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     

anova(TBFS1,TBFS2,TBFS3)
rsquared.lme(list(TBFS1)) 
rsquared.lme(list(TBFS2)) 
rsquared.lme(list(TBFS3)) 
#slope

TBFS3.1<-lme(log_TB ~ speed.corr + log_weight.kg + temp + speed.corr*temp + speed.corr*log_weight.kg,
             random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     
drop1(TBFS3.1, test="Chisq")

TBFS3.2<-lme(log_TB ~ speed.corr + log_weight.kg + temp,
             random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     
drop1(TBFS3.2, test="Chisq")

#issue is temp and mass are correlated in a strange way that seems to be only a major issue with Swimming speed. 
#choose weight
TBFS3.3<-lme(log_TB ~ speed.corr + log_weight.kg,
             random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     
summary(TBFS3.3)
rsquared.lme(list(TBFS3.3)) 

TBFS3.temp<-lme(log_TB ~ speed.corr + temp,
                random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     
summary(TBFS3.temp)
rsquared.lme(list(TBFS3.temp)) 


#final:
TBFS3.final<-lme(log_TB ~ speed.corr + log_weight.kg,
                 random = ~speed.corr - 1 | FishID, method = "REML", na.action = na.omit, data = swim.data.acc%>% filter(!is.na(TB.s))) #random slope     
summary(TBFS3.final)
rsquared.lme(list(TBFS3.final)) 

#preds
acc <- data.frame(speed.corr=seq(40, 95, 1))
mass <- data.frame(weight.kg=seq(0.3,3.5,0.2))
preds.TBFS <- merge(acc, mass)
preds.TBFS$log_weight.kg <- log(preds.TBFS$weight.kg)
preds.TBFS$TBFS.pred <- exp(predict(TBFS3.3, preds.TBFS, level = 0))


ggplot(swim.data.acc, aes(speed.corr, TB.s, col=weight.kg))+geom_point()+
  geom_line(data=preds.TBFS, aes(speed.corr, TBFS.pred, col=weight.kg, group=weight.kg))+
  theme_bw()+scale_y_continuous(limits=c(0,4.2))+scale_x_continuous(limits=c(35,100))+
  scale_color_viridis_c(end=0.9, name="Mass (kg)")+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y=expression(Tailbeat~Frequency~(beats~sec^-1)), subtitle="B) Swimming Speed")


#both tailbeat plots
ggplot(swim.data.acc, aes(gs, TB.s, col=as.factor(temp.cat2)))+geom_point()+
  geom_line(data=preds.TBF2, aes(gs, TBF.pred, col=as.factor(temp), linetype=as.factor(weight.kg)))+
  theme_bw()+scale_y_continuous(limits=c(0,4.4))+scale_x_continuous(limits=c(0,2.2))+
  scale_color_viridis_d(option="plasma", end=0.9, name="ºC")+scale_linetype(name="Mass (kg)")+
  labs(x="Acceleration (g)", y=expression(Tailbeat~Frequency~(beats~sec^-1)), title="Tailbeat Frequency", subtitle="A) Acceleration")+
  
  ggplot(swim.data.acc, aes(speed.corr, TB.s, col=weight.kg))+geom_point()+
  geom_line(data=preds.TBFS, aes(speed.corr, TBFS.pred, col=weight.kg, group=weight.kg))+
  theme_bw()+scale_y_continuous(limits=c(0,4.4))+scale_x_continuous(limits=c(35,100))+
  scale_color_viridis_c(end=0.9, name="Mass (kg)")+
  labs(x=expression(Swimming~Speed~(cm~s^-1)), y=expression(Tailbeat~Frequency~(beats~sec^-1)), subtitle="B) Swimming Speed")+
  
  plot_layout(ncol=1)
#








#Cost of Transport ----
head(swim.data)
swim.data$COT <- swim.data$massMO2/3600/(swim.data$speed.corr/100)
ggplot(swim.data, aes(speed.corr, COT, col=weight.kg))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2))
ggplot(swim.data, aes(speed.corr, COT, col=temp))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2))
ggplot(swim.data, aes(as.factor(speedcm.s), COT))+geom_boxplot()

#look at hatchery fish out of curiosity
ggplot(swim.data %>% filter(group=="hatchery control"|group=="hatchery tag"),
       aes(speed.corr, COT, col=weight.kg))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2))

#just big fish 
ggplot(swim.data %>% filter(weight.kg>2),
       aes(speed.corr, COT, col=weight.kg))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2))



#might need to look at this by temp and or size.. problem is they're correlated. 
COTm1<-lme(COT~ poly(speed.corr, 2) + log_weight.kg + temp + poly(speed.corr, 2)*temp + poly(speed.corr, 2)*log_weight.kg,
           random = ~1 | FishID, method = "REML", na.action = na.omit, data = swim.data) #random intercept
COTm2<-lme(COT ~ poly(speed.corr, 2) + log_weight.kg + temp + poly(speed.corr, 2)*temp + poly(speed.corr, 2)*log_weight.kg,
           random = ~1 + speed.corr | FishID, method = "REML", na.action = na.omit, data = swim.data) #random slope+int
COTm3<-lme(COT ~ poly(speed.corr, 2) + log_weight.kg + temp + poly(speed.corr, 2)*temp + poly(speed.corr, 2)*log_weight.kg,
           random = ~speed.corr - 1 | FishID, method = "REML", na.action = na.omit, data = swim.data) #random slope     

anova(COTm1,COTm2,COTm3)
rsquared.lme(list(COTm1)) 
rsquared.lme(list(COTm2)) 
rsquared.lme(list(COTm3)) 
#both 


COTm3.1<-lme(COT ~ poly(speed.corr, 2) + log_weight.kg + temp + poly(speed.corr, 2)*temp + poly(speed.corr, 2)*log_weight.kg,
             random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data) #random slope     
drop1(COTm3.1, test="Chisq")
#drop interactions

COTm3.1<-lme(COT ~ poly(speed.corr, 2) + log_weight.kg + temp,
             random = ~speed.corr - 1 | FishID, method = "ML", na.action = na.omit, data = swim.data) #random slope
drop1(COTm3.1, test="Chisq")

#final model
COTm3.final <-lme(COT ~ poly(speed.corr, 2) + log_weight.kg + temp,
                  random = ~speed.corr- 1 | FishID, method = "REML", na.action = na.omit, data = swim.data) #random slope
summary(COTm3.final)
rsquared.lme(list(COTm3.final)) 
exp(-0.052)

#predictions 
speeds <- data.frame(speed.corr=seq(40, 100, 1))
preds.COT <- merge(speeds, temps)
preds.COT <- merge(preds.COT, mass)
preds.COT$log_weight.kg <- log(preds.COT$weight.kg)
preds.COT$COT.pred <- predict(COTm3.final, preds.COT, level = 0)
head(preds.COT)

preds.COT2 <- preds.COT %>% filter(temp==5 & weight.kg==0.3 | temp==5 & weight.kg==3.5 |
                                     temp==11 & weight.kg==0.3 | temp==11 & weight.kg==3.5 |
                                     temp==16 & weight.kg==0.3 | temp==16 & weight.kg==3.5 |
                                     temp==21 & weight.kg==0.3 | temp==21 & weight.kg==3.5)
head(preds.COT2)

ggplot(swim.data, aes(speed.corr,COT,  col=as.factor(temp.cat2)))+geom_point()+
  geom_line(data=preds.COT2, aes(speed.corr, COT.pred, col=as.factor(temp), linetype=as.factor(weight.kg)))+
  theme_bw()+scale_y_continuous(limits=c(0,0.22))+
  scale_color_viridis_d(option="plasma", end=0.9, name="ºC")+scale_linetype(name="Weight (kg)")+
  labs(y="Cost of Transport", x=expression('Swimming Speed'~(cm~sec^-1)), title="A) Cost of Transport")

#







#plot Ucrit and COT together
ggplot(Ucrit, aes(temp, Ucrit))+geom_point()+
  geom_smooth(method="lm", formula=y~poly(x,2))+
  theme_bw()+scale_y_continuous(limits=c(0,100))+
  labs(x="Temperature (ºC)", y=expression(paste(italic(U),scriptstyle(crit),' ',( cm~s^-1))), title="A) Critical Swimming Speed")+
  
  ggplot(swim.data, aes(speed.corr,COT,  col=as.factor(temp.cat2)))+geom_point()+
  geom_line(data=preds.COT2, aes(speed.corr, COT.pred, col=as.factor(temp), linetype=as.factor(weight.kg)))+
  theme_bw()+scale_y_continuous(limits=c(0,0.22))+
  scale_color_viridis_d(option="plasma", end=0.9, name="ºC")+scale_linetype(name="Mass (kg)")+
  labs(y=expression(paste('Costs of Transport', ' (mg O'['2'], ' kg'^'-1', ' m'^'-1', ')')),
       x=expression('Swimming Speed'~(cm~sec^-1)), title="B) Costs of Transport ")+
  
  plot_layout(ncol=1)
#







#Q10 ---- 
#could try within fish that swam amongst a few temps..
swim.f <- swim.data %>% group_by(FishID, temp.cat2) %>% summarise(count=length(Trial)) %>% 
  pivot_wider(values_from=count, names_from=temp.cat2) %>% filter(!is.na(`11`) & !is.na(`16`))
head(swim.f)

swim.f2 <- swim.data %>% filter(FishID %in% swim.f$FishID)
head(swim.f2)

swim.sum <- swim.f2 %>% filter(speedcm.s==40) %>% group_by(FishID, temp.cat2) %>% summarise(massMO2=mean(massMO2)) %>% 
  pivot_wider(values_from=massMO2, names_from=temp.cat2) %>% 
  mutate(Q10=(`16`/`11`)^(10/(16-11)))
head(swim.sum)

swim.sum %>% group_by() %>% filter(!is.na(Q10)) %>%  summarise(mean(Q10))
1.8 


#just look at hatchery fish 
head(swim.data)
swim.hatch <- swim.data %>% filter(group2=="hatchery")
swim.hatch %>% filter(group2=="hatchery") %>% group_by(temp.cat2) %>% summarise(mean(massMO2))

(289 / 201)^(10/(16-11))
#2.1 for hatchery fish. 

#slow speeds
swim.hatch %>% filter(group2=="hatchery" & speedcm.s<60) %>% group_by(temp.cat2) %>% summarise(mean(massMO2))
(230 / 183)^(10/(16-11))
#1.579623

#fast
swim.hatch %>% filter(group2=="hatchery" & speedcm.s>50) %>% group_by(temp.cat2) %>% summarise(mean(massMO2))
(332 / 225)^(10/(16-11))
#2.177264















