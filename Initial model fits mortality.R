##MORTALITY MODELS

alldata<-read.table("All data for analysis.txt",h=T)

data<-alldata[alldata$subplot.type!="small",c("PlotCode","subplotID","subplot.area",
                                              "subplot.no.trees","TreeID",
                                              "binomial","region","IntervalLength",
                                              "CensusDate","dbh0","dbhgrowth","BA0","dead",
                                              "recruit","WD","subplotBA.m2ha")]

##data$dead==1, dead.
#4 trees with dbh0=0, now excluded.
mdata<-data[!is.na(data$dbh0) & data$dbh0>0 & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead),]  #check still!

mdata$subplotID<-factor(mdata$subplotID)
mdata$subplotID<-as.numeric(mdata$subplotID)
mdata$PlotCode<-factor(mdata$PlotCode)
mdata$PlotCode<-as.numeric(mdata$PlotCode)

min.fWD<-min(mdata$WD)
max.fWD<-max(mdata$WD)

##########################################################################################
##########################################################################################
##INITIAL MORTALITY MODEL

library(filzbach)

#function for predicted mortality rates
pred.mort<-function(m0,m1,s1,s3,s4,c1,c2){
  
  WD_slope <- (m1 - m0)/(max.fWD - min.fWD)
  WD_int <- m0 - (WD_slope * min.fWD)
  
  pot.longevity <- WD_int + (WD_slope * mdata$WD)
  
  #s2 <- log(99)/(s4 * (1-s3))
  #m.size <- ((mdata$dbh0/10)^s1)/(1 + exp(s2 * (mdata$dbh0 - s3 * s4)))
  m.size <- 1
  
  #m.comp <- c1 + (1 - c1) * exp(-c2 * mdata$subplotBA.m2ha)
  m.comp <- 1
  
  #k <- pot.mort * m.size * m.comp
  #p.ann <- 1 / (1+exp(-k))
  p.ann <- 1 / (1 + (pot.longevity * m.size * m.comp)) 
  p.int <- 1 - ((1-p.ann)^mdata$IntervalLength)
  
  p.int[p.int<0.0001] <- 0.0001
  p.int[p.int>0.9999] <- 0.9999
  
  return(p.int)
}

mort.ll<-function(m0,m1,s1,s3,s4,c1,c2){
  
  ll=numeric()
  
  #calculate likelihood for each tree
  for (i in nrow(mdata)){
    
    #probability of mortality
    p.int<-pred.mort(m0,m1,s1,s3,s4,c1,c2)
    
    #assign status (dead/alive) and calculate likelihood
    if (mdata[i,]$dead==1) ll[i]<-dbinom(data[i,]$dead,size=1,prob=p.int[i],log=T)
    else ll[i]<-dbinom(data[i,]$dead,size=1,prob=1-p.int[i],log=T)
    
    #sum likelihood per tree
    ll_tot<-sum(ll)
    
    if(is.na(ll_tot)) print(range(p.int))
    #if(is.na(p.int)) print("NA")
    
    return(ll_tot)
    
  }
}

fb.pars.m<-list(
  m0 = c(1e-3,1000,1,1,0,1),
  m1 = c(1e-3,1000,1,1,0,1),
  s1 = c(1e-3,10,1,1,1,1),
  s3 = c(1e-3,100,1,1,1,1),
  s4 = c(1e-3,1,1,0,1,1),
  c1 = c(1e-3,1,1,1,1,1),      #from 0 to 1 only
  c2 = c(1e-3,100,1,1,1,1)
)  

fb.out.m<-filzbach(40000,20000,mort.ll,nrow(mdata),fb.pars.m)

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##RECRUITMENT

#Number of stems/ha/yr per subplot
#Potential rate * landscape (plot) effect. A measure of WD? Spp not possible

library(filzbach)

#Calculate recruitment rate
rdata<-aggregate(data$recruit,list(data$PlotCode,data$subplotID,data$subplot.area,
                                   data$subplot.no.trees,data$region,data$IntervalLength,
                                   data$CensusDate,data$subplotBA.m2ha),sum,na.rm=T)
names(rdata)<-c("PlotCode","subplotID","subplot.area","subplot.no.trees",
                "region","IntervalLength","CensusDate","subplotBA.m2ha","sum.recruits")
rdata$recr.rate<-rdata$sum.recruits/rdata$IntervalLength

meanWD<-aggregate(data$WD,list(data$PlotCode),mean,na.rm=T)
names(meanWD)<-c("PlotCode","meanWD")

rdata2<-merge(rdata,meanWD,all.x=T)

hist(rdata2$recr.rate,breaks=100)

#function for predicted recruitment rates
pred.recr<-function(r1,w1){
  
  #potential recruitment could be a function of WD: no (aggregated over stems...)
  pot.recr = r1
  
  #landscape WD, ranging from 0 to 1, or absolute cwm?
  #Add still
  r.LWD = min(1,(w1 + (1-w1)))
  
  #subplot BA?
  #Add? Same structure as for WD, but probably not a large effect for 10 cm trees
  #Although it may work through in growth.
  
  pred = pot.recr * r.LWD
  
  return(pred)
}

#likelihood: GAMMA OR ZERO-INFLATED GAMMA?