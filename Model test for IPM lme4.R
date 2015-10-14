##########################################################################################
##########################################################################################
###MODEL TEST GROWTH

library(lme4)
library(MuMIn)

#Load/select/format data
#dropbox <- "C:/Users/DMAR/Dropbox/Rainfor data edited"
dropbox <- "C:/Users/rozendad/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/vande202/Dropbox/Rainfor data edited"

alldata<-read.table(paste(dropbox,"All data for analysis.txt",sep="/"),h=T)

data<-alldata[,c("region","PlotCode","subplotID","subplot.area",
                 "subplot.no.trees","TreeID","family",
                 "binomial","region","IntervalLength",
                 "CensusDate","dbh0","dbh1","dbhgrowth","BA0","dead",
                 "recruit","WD","subplotBA.m2ha")]

gdata<-data[!is.na(data$dbhgrowth) & !is.na(data$dbh0) & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead) & data$dbh0>0,]

gdata2<-gdata

#gdata$subplotID<-factor(gdata$subplotID)
#gdata$subplotID<-as.numeric(gdata$subplotID)
#gdata$PlotCode<-factor(gdata$PlotCode)
#gdata$PlotCode<-as.numeric(gdata$PlotCode)

#Standardize predictors
gdata$ln.dbh0<-log(gdata$dbh0)
gdata$s.ln.dbh0<-(gdata$ln.dbh0-mean(gdata$ln.dbh0))/sd(gdata$ln.dbh0)
gdata$ln.dbh0_2<-gdata$ln.dbh0^2
gdata$s.ln.dbh0_2<-(gdata$ln.dbh0_2-mean(gdata$ln.dbh0_2))/sd(gdata$ln.dbh0_2)   
gdata$s.dbh0<-(gdata$dbh0-mean(gdata$dbh0))/sd(gdata$dbh0)                           
gdata$s.WD<-(gdata$WD-mean(gdata$WD))/sd(gdata$WD)

##########################################################################################
##########################################################################################
##DATA EXPLORATION
#Look at calendar dates and plot locations.
#Include all plots.

library(maps)
library(mapdata)

cal.year<-aggregate(gdata$CensusDate,list(gdata$PlotCode),mean,na.rm=T)
names(cal.year)<-c("PlotCode","CensusDate")

##Reformat data formatting script still to include all 209 plots
#now 181 plots, min=1976, max=2009, in total 4 plots before 1980, 42 before 1990,
#109 before 2000, 133 before 2005) 

##No. plots and spp per continent
#Africa: 62 plots, 700 spp
length(unique(gdata[gdata$region=="Africa",]$PlotCode))
length(unique(gdata[gdata$region=="Africa",]$binomial))

#South America: 127 plots, 3178 spp
length(unique(gdata[gdata$region=="South.America",]$PlotCode))
length(unique(gdata[gdata$region=="South.America",]$binomial))

#Count number of stems per species
agg<-aggregate(gdata$dbh0,list(gdata$region,gdata$PlotCode,gdata$binomial),function(x) length(na.omit(x)))
names(agg)<-c("region","PlotCode","binomial","no.stems")
agg2<-agg[order(-agg$no.stems),]

##Load plot data for maps
plotdata<-read.table("Plotdata.txt",h=T)
plotdata2<-plotdata[plotdata$PlotCode %in% gdata2$PlotCode,]
plotdata2$LongitudeDecimal<--plotdata2$LongitudeDecimal

x11(10,6)
par(mfrow=c(1,2),mar=c(0,0,0,0))
map("worldHires", xlim=c(-85,-30), ylim=c(-20,15), 
    col="gray90", fill=TRUE, mar=c(0,1,0,1))
points(plotdata2$LongitudeDecimal, plotdata2$LatitudeDecimal, pch=19, cex=1.1)

map("worldHires", xlim=c(-20,50), ylim=c(-20,15), 
    col="gray90", fill=TRUE, mar=c(0,1,0,1))
points(plotdata2$LongitudeDecimal, plotdata2$LatitudeDecimal, pch=19, cex=1.1)

##########################################################################################
##########################################################################################
##MODELS TO COMPARE
#Model 1: works
model1<-lmer(dbhgrowth~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
               (ln.dbh0 + WD|PlotCode/binomial), data=gdata)
AICc(model1)          
summary(model1)
r.squaredGLMM(model1)

#Model 2: does not converge
model2<-lmer(dbhgrowth~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
               (ln.dbh0 + ln.dbh0_2 + WD|PlotCode/binomial), data=gdata)
AICc(model2)          
summary(model2)
r.squaredGLMM(model2)

#Model 3: does not converge
model3<-lmer(dbhgrowth~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
               (ln.dbh0 + WD + WD * ln.dbh0|PlotCode/binomial), data=gdata)
AICc(model3)          
summary(model3)
r.squaredGLMM(model3)

#Model 4: does not converge
model4<-lmer(dbhgrowth~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
               (1|PlotCode) + (ln.dbh0 + WD + WD * ln.dbh0|binomial), data=gdata)
AICc(model4)          
summary(model4)
r.squaredGLMM(model4)

#Model 5: works
model5<-lmer(dbhgrowth~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
               (1|PlotCode) + (ln.dbh0 + WD|binomial), data=gdata)
AICc(model5)          
summary(model5)
r.squaredGLMM(model5)

#Model 6: works
model6<-lmer(dbhgrowth~s.ln.dbh0 + s.ln.dbh0_2 + (s.ln.dbh0 + s.ln.dbh0_2|PlotCode/binomial), data=gdata)
AICc(model6)          
summary(model6)
r.squaredGLMM(model6)

newdata<-data.frame(WD=c(rep(0.3,100),rep(0.6,100),rep(0.9,100)),
                    dbh0=c(rep(seq(10,200,length.out=100),3)))
newdata$ln.dbh0<-log(newdata$dbh0)
newdata$ln.dbh0_2<-newdata$ln.dbh0^2
newdata$pred<-predict(model6,newdat=newdata,re.form=NA)
plot(newdata$dbh0,newdata$pred)

#Extract species-specific coefficients for in IPM
coefficients <- coef(model6)$binomial
coefficients$sp.plotnames<-row.names(coefficients)
vector<-strsplit(coefficients$sp.plotnames,":") #fix still

write.table(coefficients,"Coefficients for IPM test growth.txt",row.names=F,quote=F,sep="\t")
##########################################################################################
##########################################################################################
###MODEL TEST MORTALITY

library(lme4)
library(MuMIn)

#Load/select/format data
#dropbox <- "C:/Users/DMAR/Dropbox/Rainfor data edited"
dropbox <- "C:/Users/rozendad/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/vande202/Dropbox/Rainfor data edited"

alldata<-read.table(paste(dropbox,"All data for analysis.txt",sep="/"),h=T)

data<-alldata[,c("region","PlotCode","subplotID","subplot.area",
                 "subplot.no.trees","TreeID",
                 "binomial","region","IntervalLength",
                 "CensusDate","dbh0","dbh1","dbhgrowth","BA0","dead",
                 "recruit","WD","subplotBA.m2ha")]

##data$dead==1, dead.
#4 trees with dbh0=0, now excluded.
mdata<-data[!is.na(data$dbh0) & data$dbh0>0 & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead),]  #check still!

#mdata$subplotID<-factor(mdata$subplotID)
#mdata$subplotID<-as.numeric(mdata$subplotID)
#mdata$PlotCode<-factor(mdata$PlotCode)
#mdata$PlotCode<-as.numeric(mdata$PlotCode)

#Standardize predictors
mdata$ln.dbh0<-log(mdata$dbh0)
mdata$s.ln.dbh0<-(mdata$ln.dbh0-mean(mdata$ln.dbh0))/sd(mdata$ln.dbh0)
mdata$ln.dbh0_2<-mdata$ln.dbh0^2
mdata$s.ln.dbh0_2<-(mdata$ln.dbh0_2-mean(mdata$ln.dbh0_2))/sd(mdata$ln.dbh0_2)   
mdata$s.dbh0<-(mdata$dbh0-mean(mdata$dbh0))/sd(mdata$dbh0)                           
mdata$s.WD<-(mdata$WD-mean(mdata$WD))/sd(mdata$WD)

#Add in survival
mdata$surv<-ifelse(mdata$dead==0,1,0)

#To correct for IntervalLength still.

##########################################################################################
##########################################################################################
##MODELS TO COMPARE
#Model 1: does not work with predictors on original scale
model1<-glmer(dead~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
               (1|PlotCode) + (ln.dbh0 + WD|binomial), family=binomial, data=mdata)
AICc(model1)          
summary(model1)
r.squaredGLMM(model1)

#Model 2: does not work with predictors on original scale
model2<-glmer(dead~ln.dbh0 + ln.dbh0_2 + WD + WD * ln.dbh0 + WD * ln.dbh0_2 + 
                (1|PlotCode/binomial), family=binomial, data=mdata)
AICc(model2)          
summary(model2)
r.squaredGLMM(model2)

#Model 3: does not work
model3<-glmer(dead~s.ln.dbh0 + s.ln.dbh0_2 + s.WD + s.WD * s.ln.dbh0 + s.WD * s.ln.dbh0_2 + 
                (1|PlotCode), family=binomial, data=mdata)
AICc(model3)          
summary(model3)
r.squaredGLMM(model3)

#Model 4: works
model4<-glmer(dead~s.ln.dbh0 + s.ln.dbh0_2 + s.WD + s.WD * s.ln.dbh0 +  
                (1|PlotCode), family=binomial, data=mdata)
AICc(model4)          
summary(model4)
r.squaredGLMM(model4)

#Model 5: works sometimes
model5<-glmer(dead~s.ln.dbh0 + s.ln.dbh0_2 + s.WD + s.WD * s.ln.dbh0 +  
                (1|PlotCode/binomial), family=binomial, data=mdata)
AICc(model5)          
summary(model5)
r.squaredGLMM(model5)

#Model 6: does not work
model6<-glmer(dead~s.ln.dbh0 + s.ln.dbh0_2 + (s.ln.dbh0 + s.ln.dbh0_2|PlotCode/binomial), 
              family=binomial, data=mdata)
AICc(model6)          
summary(model6)
r.squaredGLMM(model6)

#Model 7: does not work
model7<-glmer(surv~ln.dbh0 + (ln.dbh0|PlotCode/binomial), 
              family=binomial, data=mdata)
AICc(model7)          
summary(model7)
r.squaredGLMM(model7)

#Model 7: works
mdata.Africa<-mdata[mdata$region=="Africa",]
model7<-glmer(surv~ln.dbh0 + (ln.dbh0|PlotCode/binomial), 
              family=binomial, data=mdata.Africa)
AICc(model7)          
summary(model7)
r.squaredGLMM(model7)

#Model 8
mdata.South.America<-mdata[mdata$region=="South.America",]
model8<-glmer(surv~ln.dbh0 + (ln.dbh0|PlotCode/binomial), 
              family=binomial, data=mdata.South.America)
AICc(model8)          
summary(model8)
r.squaredGLMM(model8)

#
newdata2<-data.frame(WD=c(rep(0.3,100),rep(0.6,100),rep(0.9,100)),
                    dbh0=c(rep(seq(10,200,length.out=100),3)))
newdata2$ln.dbh0<-log(newdata2$dbh0)
newdata2$ln.dbh0_2<-newdata2$ln.dbh0^2
newdata2$s.ln.dbh0<-(newdata2$ln.dbh0-mean(newdata2$ln.dbh0))/sd(newdata2$ln.dbh0)
newdata2$s.ln.dbh0_2<-(newdata2$ln.dbh0_2-mean(newdata2$ln.dbh0_2))/sd(newdata2$ln.dbh0_2)   
newdata2$s.dbh0<-(newdata2$dbh0-mean(newdata2$dbh0))/sd(newdata2$dbh0)                           
newdata2$s.WD<-(newdata2$WD-mean(newdata2$WD))/sd(newdata2$WD)

newdata2$pred<-predict(model5,newdat=newdata2,re.form=NA,type=c("response"))
plot(newdata2$dbh0,newdata2$pred)

#Extract species-specific coefficients for in IPM
coefficients <- coef(model8)$binomial
coefficients$sp.plotnames<-row.names(coefficients)
names(coefficients)<-c("intercept","ln.dbh0","sp.plotnames")
vector<-strsplit(coefficients$sp.plotnames,":") #fix still

#write.table(coefficients,"Coefficients for IPM test mort.txt",row.names=F,quote=F,sep="\t")
write.table(coefficients,"Coefficients for IPM test surv SA.txt",row.names=F,quote=F,sep="\t")

##########################################################################################
##########################################################################################
#Select recruitment data

alldata<-read.table("All data for analysis.txt",h=T)

data<-alldata[,c("PlotCode","subplotID","subplot.area",
                 "subplot.no.trees","TreeID",
                 "binomial","region","IntervalLength",
                 "CensusDate","dbh0","dbh1","dbhgrowth","BA0","dead",
                 "recruit","WD","subplotBA.m2ha")]

#4 trees with dbh0=0, now excluded.
rdata<-data[!is.na(data$dbh1) & data$dbh1>0 & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees),]  #check still!
