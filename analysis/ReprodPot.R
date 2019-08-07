#####################################################
####### Estimation of reproductive potential ########
#######    Andrés Beita-Jiménez, May 2019    ########
#####################################################

#read data
calib.fec<-read.csv("data/calib_Fec.csv",header = TRUE)
fec<-read.csv("data/fec2.csv",header=TRUE)
abundance<-read.csv("data/abundance.csv",header=TRUE)

#run "ObserverData" firs and then run "L50" 
modl50<-l_mat_linear2
modl50.2<-l_mat_linear3

#calibration
str(calib.fec)
calib<-lm(eggs~DryW,data=calib.fec)
calib
summary(calib)
confint(calib)

lm_eqn <- function(calib.fec){
  m <- lm(eggs~DryW,data=calib.fec);
  eq <- substitute(italic(E) == a + b %.% italic(DW)*","~~italic(r)^2~"="~r2, 
                   list(a = as.numeric(format(coef(m)[1], digits = 2)),
                        b = as.numeric(format(coef(m)[2], digits = 2)),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ggplot(calib.fec,aes(x=DryW,y=eggs))+
  geom_smooth(method = "lm",colour="blue")+geom_point(size=2,alpha=0.6)+
  theme_bw()+
  geom_text(x = 0.25, y = 1700, label = lm_eqn(calib.fec), parse = TRUE)+
  ylab("Number of Eggs")+xlab("Dry Weight (g)")+
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/Calibration.png", width = 20, height = 12, units = "cm")

#fecundity
eggs<-predict(calib,fec)
fec2<-cbind(fec,eggs)

str(fec2)
#subset(fec2,is.nan(eggs))
mod.fec<-lm(log(eggs)~log(LC),data=fec2)

summary(mod.fec)
confint(mod.fec)

lm_eqn <- function(fec2){
  m <- lm(eggs~LC,data=fec2);
  eq <- substitute(italic(E) == a + b %.% italic(LC)*","~~italic(r)^2~"="~r2, 
                   list(a = as.numeric(format(coef(m)[1], digits = 2)),
                        b = as.numeric(format(coef(m)[2], digits = 2)),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ggplot(fec2,aes(x=LC,y=eggs))+
  geom_smooth(method = "lm",colour="blue")+geom_point(size=2,alpha=0.6)+
  theme_bw()+
  geom_text(x = 22, y = 1750, label = lm_eqn(fec2), parse = TRUE)+
  ylab("Number of Eggs")+xlab("LC (mm)")+
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/fec.png", width = 20, height = 12, units = "cm")

#fecundity excluding first site
fec3<-subset(fec2,Shrimp>68)

str(fec2)
mod.fec2<-lm(log10(eggs)~log10(LC),data=fec3)


summary(mod.fec2)
confint(mod.fec2)

a<-10^mod.fec2$coefficients[1]
b<-mod.fec2$coefficients[2]

lm_eqn <- function(fec3){
  m <- lm(log10(eggs)~log10(LC),data=fec3);
  eq <- substitute(italic(Log[10](E)) == a + b %.% italic(Log[10](CL))*","~~italic(r)^2~"="~r2, 
                   list(a = as.numeric(format(coef(m)[1], digits = 3)),
                        b = as.numeric(format(coef(m)[2], digits = 3)),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ggplot(fec3,aes(x=LC,y=eggs))+
  geom_smooth(method = "lm",formula= y~I(a*(x^b)), colour="blue")+geom_point(size=2,alpha=0.6)+
  theme_bw()+
  geom_text(x = 22.5, y = 1820, label = lm_eqn(fec3), parse = TRUE)+
  ylab("Number of Eggs")+xlab("LC (mm)")+
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/fecundityDef.png", width = 20, height = 12, units = "cm")


#length weight relation
str(fec3)
mod.lw<-lm(log10(W)~log10(LC),data=fec3,na.action=na.omit)

summary(mod.lw)
confint(mod.lw)

a<-10^mod.lw$coefficients[1]
b<-mod.lw$coefficients[2]

lm_eqn <- function(fec3){
  m <- lm(log10(W)~log10(LC),data=fec3,na.action=na.omit);
  eq <- substitute(italic(log[10](W)) == a + b %.% italic(log[10](CL))*","~~italic(r)^2~"="~r2, 
                   list(a = as.numeric(format(coef(m)[1], digits = 2)),
                        b = as.numeric(format(coef(m)[2], digits = 2)),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ggplot(fec3,aes(x=LC,y=W))+
  geom_smooth(method = "lm",formula= y~I(a*(x^b)),colour="blue")+geom_point(size=2,alpha=0.6)+
  theme_bw()+
  geom_text(x = 22.5, y = 14.3, label = lm_eqn(fec3), parse = TRUE)+
  ylab("Weight (g)")+xlab("LC (mm)")+
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/LW_def.png", width = 20, height = 12, units = "cm")


#estimate proportion at size
str(expand.explor)

prop.size<-expand.explor %>%
  subset(Id == "6") %>%
  group_by(year,Id,size_mm) %>%
  summarise(conteo=n()) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(tot=sum(conteo))

prop.size$prop<-prop.size$conteo/prop.size$tot
prop.size$LC<-prop.size$size_mm

#join abundance data
str(abundance)
abundance$Id<-as.factor(abundance$Id)
prop.size<-left_join(prop.size,abundance,by=c("Id","year"))

#function to estimate fecundity in 1974 
fec1974<-function(x){
  log10_res<-3.4614*log10(x)-1.6670
  res<-10^log10_res
  res
}

fec1974nov<-function(x){
  log10_res<-3.0106*log10(x)-1.0147
  res<-10^log10_res
  res
}

summary(mod.fec2)

fec2018<-function(x){
  log10_res<-mod.fec2$coefficients[2]*log10(x)+mod.fec2$coefficients[1]
  res<-10^log10_res
  res
}
#SSB
prop.size<-subset(prop.size,year %in% 1996:2016)
prop.size$year<-as.factor(prop.size$year)
transition<-predict(modl50,prop.size, type="response")
prop.size$P<-transition
WatL<-predict(mod.lw,prop.size, type="response")
prop.size$WatL<-10^WatL
prop.size$rel.ssb<-prop.size$prop*prop.size$P*prop.size$WatL
prop.size$ssb<-prop.size$rel.ssb*prop.size$abundance
Fec<-predict(mod.fec2,prop.size, type="response")
fec1<-fec2018(prop.size$size_mm)
prop.size$Fec<-10^Fec
fec2<-fec1974(prop.size$size_mm)
prop.size$Fec2<-fec2
fec3<-fec1974nov(prop.size$size_mm)
prop.size$Fec3<-fec3
prop.size$EP<-prop.size$prop*prop.size$P*prop.size$Fec*prop.size$abundance
prop.size$EP2<-prop.size$prop*prop.size$P*prop.size$Fec2*prop.size$abundance
prop.size$EP3<-prop.size$prop*prop.size$P*prop.size$Fec3*prop.size$abundance

rep.pot<-prop.size %>%
  group_by(year,Id) %>%
  summarise(SSB=sum(ssb),TEP=sum(EP),TEP2=sum(EP2),TEP3=sum(EP3))

ggplot(rep.pot,aes(x=as.integer(as.character(year)),SSB, group=1))+
  geom_line()+
  #facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("SSB")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/SSB.png", width = 20, height = 18, units = "cm")

ggplot(rep.pot,aes(x=as.integer(as.character(year)),TEP, group=1))+
  geom_line()+
  #facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("Total Egg Prod")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/TotalEggProd.png", width = 20, height = 18, units = "cm")

str(prop.size)

#jpeg("plots/ReprodPotential.png", width = 600, height = 350)
tiff("plots/ReprodPotential.tiff", width = 9, height = 5, units = 'in', res = 300)
par(mar = c(5, 5, 3, 5))
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$SSB/1000000000, type="l",
     xlab="Year",ylab=expression(SSB (x10^9)), col = "black",lwd=2,cex.lab=1.5,cex.axis=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP/1000000000, type="l",
     xaxt="n",yaxt="n", col = "red",lty=2,
     ylab = "", xlab = "",lwd=2,cex.axis=1.5)
axis(side = 4,cex.axis=1.5)
mtext(expression(TEP (x10^11)), side = 4, line = 3, cex=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP2/1000000000, type="l",
     xaxt="n",yaxt="n", col = "blue",lty=2,
     ylab = "", xlab = "",lwd=2)
legend("topright", c("SSB","TEP 1974", "TEP 2018"),
       col = c("black","blue", "red"), lty = c(1, 2, 2), 
       cex = 1.5,lwd=2,bty="n")

dev.off()
#reference points

#install.packages("EnvStats")
library(EnvStats)
refpoin<-subset(rep.pot,year %in% 1996:2003)
#SSB
gmssb<- geoMean(refpoin$SSB)
USRssb<- gmssb*0.8
LRPssb<- gmssb*0.3


#TEP1978
gm<- geoMean(refpoin$TEP2)
USRtep1978<- gm*0.8
LRPtep1978<- gm*0.3

#TEP2018
gm<- geoMean(refpoin$TEP)
USRtep2018<- gm*0.8
LRPtep2018<- gm*0.3

rp.18 <- ggplot(rep.pot,aes(TEP/100000000000,SSB/1000000000))+
  geom_point(size=2,alpha=0.5)+
  theme_bw()+
  geom_hline(yintercept = USRssb/1000000000, colour="green")+
  geom_hline(yintercept = LRPssb/1000000000, colour="red")+
  geom_vline(xintercept = USRtep2018/100000000000, colour="green")+
  geom_vline(xintercept = LRPtep2018/100000000000, colour="red")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab(expression("SSB "(x10^9)))+xlab(expression("TEP "(x10^11)))

rp.18

ggsave("plots/refpoint2018.png", width = 20, height = 12, units = "cm")

rp.78 <- ggplot(rep.pot,aes(TEP2/100000000000,SSB/1000000000))+
  geom_point(size=2,alpha=0.5)+
  theme_bw()+
  geom_hline(yintercept = USRssb/1000000000, colour="green")+
  geom_hline(yintercept = LRPssb/1000000000, colour="red")+
  geom_vline(xintercept = USRtep1978/100000000000, colour="green")+
  geom_vline(xintercept = LRPtep1978/100000000000, colour="red")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab(expression("SSB "(x10^9)))+xlab(expression("TEP "(x10^11)))

rp.78

ggsave("plots/refpoint1978.png", width = 20, height = 12, units = "cm")



####### Include L50 from observer data #####
#SSB

#prop.size<-subset(prop.size,year %in% 1996:2016)
#prop.size$year<-as.factor(prop.size$year)
transition<-predict(modl50.2,prop.size, type="response")
prop.size$P<-transition
#WatL<-predict(mod.lw,prop.size, type="response")
#prop.size$WatL<-WatL
prop.size$rel.ssb<-prop.size$prop*prop.size$P*prop.size$WatL
prop.size$ssb.2<-prop.size$rel.ssb*prop.size$abundance
#Fec<-predict(mod.fec2,prop.size, type="response")
#fec1<-fec2018(prop.size$size_mm)
#prop.size$Fec<-Fec
#fec2<-fec1974(prop.size$size_mm)
#prop.size$Fec2<-fec2
#fec3<-fec1974nov(prop.size$size_mm)
#prop.size$Fec3<-fec3
prop.size$EP.2<-prop.size$prop*prop.size$P*prop.size$Fec*prop.size$abundance
prop.size$EP2.2<-prop.size$prop*prop.size$P*prop.size$Fec2*prop.size$abundance
prop.size$EP3.2<-prop.size$prop*prop.size$P*prop.size$Fec3*prop.size$abundance

rep.pot<-prop.size %>%
  group_by(year,Id) %>%
  summarise(SSB=sum(ssb),TEP=sum(EP),TEP2=sum(EP2),TEP3=sum(EP3),
            SSB.2=sum(ssb.2),TEP.2=sum(EP.2),TEP2.2=sum(EP2.2),TEP3.2=sum(EP3.2))

ggplot(rep.pot,aes(x=as.integer(as.character(year)),SSB.2, group=1))+
  geom_line()+
  #facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("SSB")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/SSB_obser.png", width = 20, height = 18, units = "cm")

ggplot(rep.pot,aes(x=as.integer(as.character(year)),TEP.2, group=1))+
  geom_line()+
  #facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("Total Egg Prod")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/TotalEggProd_observ.png", width = 20, height = 18, units = "cm")

str(prop.size)

#jpeg("plots/ReprodPotential.png", width = 600, height = 350)
tiff("plots/ReprodPotential_observ.tiff", width = 9, height = 5, units = 'in', res = 300)
par(mar = c(5, 5, 3, 5))
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$SSB.2/1000000000, type="l",
     xlab="Year",ylab=expression(SSB (x10^9)), col = "black",lwd=2,cex.lab=1.5,cex.axis=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP.2/1000000000, type="l",
     xaxt="n",yaxt="n", col = "red",lty=2,
     ylab = "", xlab = "",lwd=2,cex.axis=1.5)
axis(side = 4,cex.axis=1.5)
mtext(expression(TEP (x10^11)), side = 4, line = 3, cex=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP2.2/1000000000, type="l",
     xaxt="n",yaxt="n", col = "blue",lty=2,
     ylab = "", xlab = "",lwd=2)
legend("topright", c("SSB","TEP 1974", "TEP 2018"),
       col = c("black","blue", "red"), lty = c(1, 2, 2), 
       cex = 1.5,lwd=2,bty="n")

dev.off()

tiff("plots/SSB_observ-survey.tiff", width = 9, height = 5, units = 'in', res = 300)
par(mar = c(5, 5, 3, 5))
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$SSB.2/1000000000, type="l",
     xlab="Year",ylab=expression(SSB (x10^9)), col = "black",lwd=2,cex.lab=1.5,cex.axis=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$SSB/1000000000, type="l",
     xaxt="n",yaxt="n", col = "red",lty=1,
     ylab = "", xlab = "",lwd=2,cex.axis=1.5)
legend("topright", c("Survey","Observer"),
       col = c("black", "red"), lty = c(1, 1), 
       cex = 1.5,lwd=2,bty="n",title = "Data")

dev.off()


#reference points acording to datatype
refpoin2<-subset(rep.pot,year %in% 1996:2003)

#SSB observer data
gmssb.2<- geoMean(refpoin2$SSB.2)
USRssb.2<- gmssb.2*0.8
LRPssb.2<- gmssb.2*0.3

rp.ssb <- ggplot(rep.pot,aes(SSB/1000000000,SSB.2/1000000000))+
  geom_point(size=2,alpha=0.5)+
  theme_bw()+
  geom_hline(yintercept = USRssb.2/1000000000, colour="green")+
  geom_hline(yintercept = LRPssb.2/1000000000, colour="red")+
  geom_vline(xintercept = USRssb/1000000000, colour="green")+
  geom_vline(xintercept = LRPssb/1000000000, colour="red")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab(expression("SSB"[o] (x10^9)))+xlab(expression("SSB"[s] (x10^9)))

rp.ssb

ggsave("plots/refpointSSB_datatype.png", width = 20, height = 12, units = "cm")

ggarrange(rp.78,rp.18,rp.ssb,labels=c("A","B","C"),ncol = 1, nrow=3)

ggsave("plots/refpointall1col.png", width = 13, height = 24, units = "cm")

ggarrange(rp.78,rp.18,rp.ssb,labels=c("A","B","C"))
ggsave("plots/refpointall2col.png", width = 30, height = 21, units = "cm")

ggarrange(rp.18,rp.ssb,labels=c("A","B"),ncol = 1, nrow=2)
ggsave("plots/refpoint2018_ssb.png", width = 17, height = 22, units = "cm")

ggarrange(rp.78,rp.18,labels=c("A","B"),ncol = 1, nrow=2)
ggsave("plots/refpoint1978_2018.png", width = 17, height = 22, units = "cm")

##################################################

par(mar = c(5, 5, 3, 5))
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP2, type="l",
     xlab="Year",ylab="TEP", col = "black",lwd=2)
par(new=TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP, type="l",
     xaxt="n",yaxt="n", col = "red",lty=2,
     ylab = "", xlab = "",lwd=2)
legend("topleft", c("1974", "2018"),
       col = c("black", "red"), lty = c(1, 2))

summary(mod.fec2)

fec2018<-function(x){
  log10_res<-mod.fec2$coefficients[2]*log10(x)+mod.fec2$coefficients[1]
  res<-10^log10_res
  res
}

summary(mod.fec2)
a<-10^mod.fec2$coefficients[1]
b<-mod.fec2$coefficients[2]

db<-crossing(LC=seq(19,30,by=0.1))
eggs<-predict(mod.fec2,db)
preds <- predict(mod.fec2, newdata = db, type = "response", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- 10^(preds$fit + (critval * preds$se.fit))
lwr <- 10^(preds$fit - (critval * preds$se.fit))
fit <- 10^(preds$fit)
dat.fec<-data.frame(LC=db$LC,eggs=fit,upr=upr,lwr=lwr)

ggplot(dat.fec,aes(x=LC,y=eggs))+
  geom_line(linetype = "dashed",colour="blue", size=1)+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=LC), fill = "blue", alpha = 0.3)+
  stat_function(fun=fec1974, geom="line",linetype = "dashed",colour="red", size=1)+
  xlab("Size (mm)") + ylab("Number of eggs")+
  theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_text(x = 28, y = 1300, label = "2018",colour="blue", parse = TRUE)+
  geom_text(x = 25, y = 1750, label = "1978",colour="red", parse = TRUE)
  
ggsave("plots/comp_fec.png", width = 20, height = 12, units = "cm") 
  

ggplot(data.frame(x=c(19, 30)), aes(x=x)) + 
  stat_function(fun=fec2018, geom="line",linetype = "dashed",colour="blue", size=1) +
  #geom_smooth(,method = "lm",formula= y~I(a*(x^b)),colour="blue",linetype = "dashed")+
  stat_function(fun=fec1974, geom="line",linetype = "dashed",colour="red", size=1) +
  xlab("Size (mm)") + ylab("Number of eggs")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




#relative
#SSB
prop.size<-subset(prop.size,year %in% 1996:2016)
prop.size$year<-as.factor(prop.size$year)
transition<-predict(modl50,prop.size, type="response")
prop.size$P<-transition
WatL<-predict(mod.lw,prop.size, type="response")
prop.size$WatL<-10^WatL
prop.size$rel.ssb<-prop.size$prop*prop.size$P*prop.size$WatL
#prop.size$ssb<-prop.size$rel.ssb*prop.size$abundance
Fec<-predict(mod.fec2,prop.size, type="response")
#fec1<-fec2018(prop.size$size_mm)
prop.size$Fec<-10^Fec
fec2<-fec1974(prop.size$size_mm)
prop.size$Fec2<-fec2
fec3<-fec1974nov(prop.size$size_mm)
prop.size$Fec3<-fec3
prop.size$EP<-prop.size$prop*prop.size$P*prop.size$Fec#*prop.size$abundance
prop.size$EP2<-prop.size$prop*prop.size$P*prop.size$Fec2#*prop.size$abundance
prop.size$EP3<-prop.size$prop*prop.size$P*prop.size$Fec3#*prop.size$abundance

rep.pot<-prop.size %>%
  group_by(year,Id) %>%
  summarise(SSB=sum(rel.ssb),TEP=sum(EP),TEP2=sum(EP2),TEP3=sum(EP3))

ggplot(rep.pot,aes(x=as.integer(as.character(year)),SSB, group=1))+
  geom_line()+
  #facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("SSB")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/relSSB.png", width = 20, height = 18, units = "cm")

ggplot(rep.pot,aes(x=as.integer(as.character(year)),TEP, group=1))+
  geom_line()+
  #facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("Total Egg Prod")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/relTotalEggProd.png", width = 20, height = 18, units = "cm")

str(prop.size)

#jpeg("plots/ReprodPotential.png", width = 600, height = 350)
tiff("plots/RelativeReprodPotential.tiff", width = 9, height = 5, units = 'in', res = 300)
par(mar = c(5, 5, 3, 5))
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$SSB, type="l",
     xlab="Year",ylab="Relative SSB", col = "black",lwd=2,cex.lab=1.5,cex.axis=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP, type="l",
     xaxt="n",yaxt="n", col = "red",lty=2,
     ylab = "", xlab = "",lwd=2,cex.axis=1.5)
axis(side = 4,cex.axis=1.5)
mtext("Relative EP", side = 4, line = 3, cex=1.5)
par(new = TRUE)
plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$TEP2, type="l",
     xaxt="n",yaxt="n", col = "blue",lty=2,
     ylab = "", xlab = "",lwd=2)
legend("bottomright", c("SSB","EP 1974", "EP 2018"),
       col = c("black","blue", "red"), lty = c(1, 2, 2), 
       cex = 1.5,lwd=2,bty="n")

dev.off()

rep.pot$dif<-rep.pot$TEP2-rep.pot$TEP

plot(x=as.integer(as.character(rep.pot$year)),y=rep.pot$dif, type="l")

