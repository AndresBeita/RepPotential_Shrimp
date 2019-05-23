#####################################################
####### Estimation of reproductive potential ########
#######    Andrés Beita-Jiménez, May 2019    ########
#####################################################

#read data
calib.fec<-read.csv("data/calib_Fec.csv",header = TRUE)
fec<-read.csv("data/fec2.csv",header=TRUE)

#run "L50" first
modl50<-l_mat_linear2

#calibration
str(calib.fec)
calib<-lm(eggs~DryW,data=calib.fec)
calib
summary(calib)

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
mod.fec<-lm(eggs~LC,data=fec2)

summary(mod.fec)

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

#length weight relation
str(fec2)
mod.lw<-lm(W~LC,data=fec2,na.action=na.omit)

summary(mod.lw)

lm_eqn <- function(fec2){
  m <- lm(W~LC,data=fec2,na.action=na.omit);
  eq <- substitute(italic(W) == a + b %.% italic(LC)*","~~italic(r)^2~"="~r2, 
                   list(a = as.numeric(format(coef(m)[1], digits = 2)),
                        b = as.numeric(format(coef(m)[2], digits = 2)),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ggplot(fec2,aes(x=LC,y=W))+
  geom_smooth(method = "lm",colour="blue")+geom_point(size=2,alpha=0.6)+
  theme_bw()+
  geom_text(x = 22, y = 14, label = lm_eqn(fec2), parse = TRUE)+
  ylab("Weight (g)")+xlab("LC (mm)")+
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/LW2.png", width = 20, height = 12, units = "cm")


#estimate proportion at size
str(expand.explor)
prop.size<-expand.explor %>%
  group_by(year,Id,size_mm) %>%
  summarise(conteo=n()) %>%
  ungroup() %>%
  group_by(year,Id) %>%
  mutate(tot=sum(conteo))

prop.size$prop<-prop.size$conteo/prop.size$tot*100
prop.size$LC<-prop.size$size_mm


#SSB
prop.size<-subset(prop.size,year %in% 1996:2016)
prop.size$year<-as.factor(prop.size$year)
transition<-predict(modl50,prop.size, type="response")
prop.size$P<-transition
WatL<-predict(mod.lw,prop.size, type="response")
prop.size$WatL<-WatL
prop.size$rel.ssb<-prop.size$prop*prop.size$P*prop.size$WatL
Fec<-predict(mod.fec,prop.size, type="response")
prop.size$Fec<-Fec
prop.size$EP<-prop.size$prop*prop.size$P*prop.size$Fec

rep.pot<-prop.size %>%
  group_by(year,Id) %>%
  summarise(r.ssb=sum(rel.ssb),r.EP=sum(EP))

ggplot(rep.pot,aes(x=as.integer(as.character(year)),r.ssb, group=1))+
  geom_line()+
  facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("Relative SSB")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/rel.SSB.png", width = 20, height = 18, units = "cm")

ggplot(rep.pot,aes(x=as.integer(as.character(year)),r.EP, group=1))+
  geom_line()+
  facet_wrap("Id",ncol=1)+
  theme_bw()+
  ylab("Relative Egg Prod")+xlab("Year")+
  theme(text = element_text(size=16),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("plots/rel.EggProd.png", width = 20, height = 18, units = "cm")
