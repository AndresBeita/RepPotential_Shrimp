######################################################
######## Estimation of size of first transition ######
########    Andres Beita-Jiménez, May 2019      ######
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
explor<-read.csv("data/pandalus_borealis_length_freq_2017-10-27.csv",header = TRUE)
years<-1996:2016

#create copy of the data to work with
spatial<-explor


#read shapefile
sfas <- readOGR("./sfa2", "SFA4to7_polygons")
#plot to check
#plot to check
plot(sfas)
text(sfas,sfas$SFA) #SFA 4 and 5 are inverted
#change it
sfas$SFA<-ifelse(sfas$SFA==5,4,ifelse(sfas$SFA==4,5,sfas$SFA))
#plot to check
plot(sfas)
text(sfas,sfas$SFA)#now is good

# specify that long and lat are coordinates
coordinates(spatial) <- ~ long_dec + lat_dec

# Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
proj4string(spatial) <- proj4string(sfas)

#create vector of clasification
start_time <- Sys.time()

sfa<-over(spatial, sfas)

end_time <- Sys.time()
end_time - start_time

#create data binding both
dat<-cbind(explor,sfa)

str(dat)
unique(dat$SFA)
dat$Id<-as.factor(dat$SFA)

#########################################################
################       data management      #############
#########################################################

#subset fall data in sfas 5 and 6
dat2<-subset(dat,Id %in% c("5","6"))
dat2<-subset(dat2,season=="Fall")

#create code to later join
dat2$code <- paste("SFA",dat2$Id,"-",dat2$year)
unique(dat2$code)

#correct total by newratio
#dat2$total.corr <- dat2$total*dat2$newratio

#calculate total by set
explor4<-dat2 %>%
  dplyr::select(year,Id,trip_id,long_dec,lat_dec,sex,size_mm,count,code,newratio) %>%
  group_by(year,trip_id,sex,size_mm) %>%
  mutate(tot=sum(count))
  
#select variables of interest
explor5<-explor4 %>%
  group_by(year,Id,trip_id,long_dec,lat_dec,sex,size_mm,code,newratio) %>%
  summarize(total=mean(tot))

#estimate density by year
dens <- explor5 %>%
  group_by(year,Id,trip_id,code,newratio) %>%
  summarize(total.trip=sum(total))

dens$total.corr=dens$total.trip*dens$newratio

dens <- dens %>%
  group_by(year,Id,code) %>%
  summarise(dens.year=mean(total.corr))

#expand to proper long format
expand.explor<-expandRows(explor5, "total")
expand.explor$sex.binom<-ifelse(expand.explor$sex=="female",1,0)


#plot
#ggplot(expand.explor,aes(x=size_mm,y=sex.binom))+
#  geom_point()+facet_wrap("year")

####### environmental variables ########
#load data
sal<-read.csv("data/salinity_anomalies.csv",header = TRUE)
temp <- read.csv("data/temperature_anomalies.csv",header = TRUE)

#join environmental variables
env <- left_join(sal,temp,by="code")
str(env)


data2<-left_join(expand.explor,env,by="code")

#join density
den<-dens %>%
  ungroup() %>%
  dplyr::select(code,dens.year)
den
data2 <- left_join(data2,den,by="code")
data2$year<-as.factor(data2$year)
str(data2)

#plot density
#png("plots/gam/density2.png",width = 800, height = 400)
dd<-ggplot(subset(dens,year %in% years),aes(x=year, y=dens.year,colour=Id ))+
  geom_line(size=1)+theme_bw()+
  theme(text = element_text(size=16),legend.position="none")+
  ylab("Individuals per tow")+xlab("Year")
dd
#dev.off()

#########################################
######## GLM transition #################
#########################################
str(data2)

data3 <- data2 %>%
  subset(year %in% 1996:2016) %>%
  subset(Id == "6")

options(na.action = "na.fail")
l_mat_linear2 = glm(sex.binom ~ size_mm +year,#*Id,#+Sal.all+Temp.all+dens.year,
                    family=binomial(link = "logit"),
                    data = data3)
summary(l_mat_linear2)

install.packages("glmmTMB")
library(glmmTMB)

prueba<-glmmTMB(sex.binom ~ size_mm +year +ar1(1/year),
             family=binomial,
             data = data3)
