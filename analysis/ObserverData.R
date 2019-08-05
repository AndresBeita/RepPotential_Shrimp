######################################################
######## Estimation of size of first transition ######
########            observer data               ######
########    Andres Beita-Jiménez, July 2019     ######
######################################################

#load libraries
library(rgdal)
library(tidyr)
library(dplyr)
library(ggplot2)
library(marmap)
library(raster)
library(mapproj)
library(reshape2)
library(mgcv)
library(splitstackshape)
library(broom)
library( MuMIn)
library(ggpubr)
#load data
explor<-read.csv("data/ObserverData_shrimp.csv",header = TRUE,check.names=FALSE)
str(explor)
years<-1990:2018

#create copy of the data to work with
dat<-explor


#########################################################
################       data management      #############
#########################################################

#subset fall data in sfas 5 and 6
#dat2<-subset(dat,Id %in% c("5","6"))
#dat2<-subset(dat2,season=="Fall")

#create code to later join
#dat2$code <- paste("SFA",dat2$Id,"-",dat2$year)
#unique(dat2$code)

#correct total by newratio
#dat2$total.corr <- dat2$total*dat2$newratio

#reshape data
dim(dat)
dat2<-gather(dat,"size_mm","count",15:94)
str(dat2)
dat2$size_mm<-as.numeric(dat2$size_mm)
dat2$trip_id<-paste(dat2$year,dat2$month,dat2$day,dat2$starthr,dat2$startmin,dat2$sfa,sep="-")

#calculate total by set
explor4<-dat2 %>%
  dplyr::select(year,sfa,trip_id,long,lat,sex,size_mm,count) %>%
  group_by(year,trip_id,sex,size_mm) %>%
  mutate(tot=sum(count))

#select variables of interest
explor5<-explor4 %>%
  group_by(year,sfa,trip_id,long,lat,sex,size_mm) %>%
  summarize(total=mean(tot))

#estimate density by year
#dens <- explor5 %>%
#  group_by(year,Id,trip_id,code,newratio) %>%
#  summarize(total.trip=sum(total))

#dens$total.corr=dens$total.trip*dens$newratio

#dens <- dens %>%
#  group_by(year,Id,code) %>%
#  summarise(dens.year=mean(total.corr))

#expand to proper long format
explor6<-subset(explor5,sfa %in% c("4","5","6") & total > 0)
expand.explor<-expandRows(explor6, "total")
expand.explor$sex.binom<-ifelse(expand.explor$sex=="female",1,0)


#plot
#ggplot(expand.explor,aes(x=size_mm,y=sex.binom))+
#  geom_point()+facet_wrap("year")

####### environmental variables ########
#load data
#sal<-read.csv("data/salinity_anomalies.csv",header = TRUE)
#temp <- read.csv("data/temperature_anomalies.csv",header = TRUE)

#join environmental variables
#env <- left_join(sal,temp,by="code")
#str(env)


#data2<-left_join(expand.explor,env,by="code")

#join density
#den<-dens %>%
#  ungroup() %>%
#  dplyr::select(code,dens.year)
#den
#data2 <- left_join(data2,den,by="code")
#data2$year<-as.factor(data2$year)
#str(data2)

#plot density
#png("plots/gam/density2.png",width = 800, height = 400)
#dd<-ggplot(subset(dens,year %in% years),aes(x=year, y=dens.year,colour=Id ))+
#  geom_line(size=1)+theme_bw()+
#  theme(text = element_text(size=16),legend.position="none")+
#  ylab("Individuals per tow")+xlab("Year")
#dd
#dev.off()

#########################################
######## GLM transition #################
#########################################
data2<-expand.explor
str(data2)

data3 <- data2 %>%
  #subset(year %in% 1990:2016) %>%
  subset(sfa == 6)

data3$year<-as.factor(data3$year)
options(na.action = "na.fail")
l_mat_linear3 = glm(sex.binom ~ size_mm +year,#*sfa,#+Sal.all+Temp.all+dens.year,
                    family=binomial(link = "logit"),
                    data = data3)
summary(l_mat_linear3)
