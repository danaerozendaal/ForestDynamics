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

min.WD<-min(mdata$WD)
max.WD<-max(mdata$WD)

mdata$log.dbh0<-log(mdata$dbh0/20)

#Put values for all subplots into a list
split.mdata <- split(mdata,mdata$subplotID)

fDBH<-split(mdata$dbh0,mdata$subplotID)
fWD<-split(mdata$WD,mdata$subplotID)
PlotCode<-split(mdata$PlotCode,mdata$subplotID)
dbhgrowth<-split(mdata$dbhgrowth,mdata$subplotID)
log.fDBH<-split(mdata$log.dbh0,mdata$subplotID)

##########################################################################################
##########################################################################################
##INITIAL MORTALITY MODEL

library(filzbach)

comp.fun <- function(target,c2) {
  
  if (!is.finite(min(exp(c2*target)))) print(range(target))
  
  #return((target/20)^c2)
  return(exp(c2 * target))
  
}

subplot.comp <- function(dbhs,c0,c1,c2) {
  
  #competition effect matrix
  comp <- outer(dbhs,rep(c2,length(dbhs)),FUN=comp.fun)
  
  #growth effect on focal trees
  fcomp <- exp(-c0 * ((exp(dbhs))^c1) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  
}

#function for predicted mortality rates
pred.mort<-function(m0,m1,s1,s2,c0,c1,c2,E){
#pred.mort<-function(m0,m1,s1,s2,c0,c1,c2){
  
  #Potential longevity
  WD_slope <- (m1 - m0)/(max.WD - min.WD)
  WD_int <- m0 - (WD_slope * min.WD)
  
  pot.longevity <- WD_int + (WD_slope * mdata$WD)
  
  #Size effect
  #s2 <- log(99)/(s4 * (1-s3))
  #m.size <- ((mdata$dbh0/10)^s1)/(1 + exp(s2 * (mdata$dbh0 - s3 * s4)))
  m.size <- ((mdata$dbh0/202)^s1) * exp(-s2 * mdata$dbh0)
  #m.size <- 1
  
  #Competition effect
  #m.comp <- c1 + (1 - c1) * exp(-c2 * mdata$subplotBA.m2ha)
  #m.comp <- 1
  m.comp = unlist(lapply(log.fDBH,FUN=subplot.comp,c0,c1,c2))
  
  #k <- pot.mort * m.size * m.comp
  #p.ann <- 1 / (1+exp(-k))
  p.ann <- 1 / (1 + (pot.longevity * m.size * m.comp * E)) 
  #p.ann <- 1 / (1 + (pot.longevity * m.size * m.comp)) 
  p.int <- 1 - ((1-p.ann)^mdata$IntervalLength)
  
  p.int[p.int<0.0001] <- 0.0001
  p.int[p.int>0.9999] <- 0.9999
  
  return(p.int)
}

mort.ll<-function(m0,m1,s1,s2,c0,c1,c2,E_all,E_mean,E_sd){
#mort.ll<-function(m0,m1,s1,s2,c0,c1,c2){
  
  #probability of mortality
  p.int<-pred.mort(m0,m1,s1,s2,c0,c1,c2,E_all[mdata$PlotCode])
  #p.int<-pred.mort(m0,m1,s1,s2,c0,c1,c2)
  
  #assign status (dead/alive) and calculate likelihood
  ll <- log(ifelse(mdata$dead==1,p.int,1-p.int))
      
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  
  #sum likelihood per tree
  ll_tot<-sum(ll) + log_E_hier
  #ll_tot<-sum(ll)
  
  if(is.na(ll_tot)) print(range(p.int))
  #if(is.na(p.int)) print("NA")
  
  return(ll_tot)
  
}

fb.pars.m <- list(
  m0 = c(1e-3,1000,1,1,0,1),
  m1 = c(1e-3,1000,1,1,0,1),
  #s1_0 = c(1e-6,3,1,0,1,1),  
  #s1_1 = c(1e-6,3,1,0,1,1), 
  #s2_0 = c(1e-6,10,1,0,1,1),
  #s2_1 = c(1e-6,10,1,0,1,1),
  s1 = c(1e-6,3,1,0,0,1),
  s2 = c(1e-6,10,1,0,0,1),
  c0 = c(-50,50,1,0,0,1),      
  c1 = c(-5,5,0,0,0,1),
  c2 = c(-5,5,2,0,0,1),
  #c2_int = c(-5,5,2,0,0,1),
  #c2_slope = c(-5,5,2,0,0,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1)
)

fb.out.m<-filzbach(40000,20000,mort.ll,nrow(mdata),fb.pars.m)
df.fb.out<-as.data.frame(fb.out.m)
write.table(df.fb.out,"FB output mort model 1.txt",row.names=F,quote=F,sep="\t")

final_out <- paste(readLines("C:/Users/rozendad/Dropbox/Current projects/UofR/ForestDynamics model runs/workspace/Default_MCMC_final_out.txt"), collapse="\t")
write.table(final_out,"Final_out mort model 1.txt",row.names=F,quote=F,sep="\t")

#Convergence (converged)
pdf("Model mort 1.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(5,4,1,1))
mort.llvec<-function(x) mort.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7])
fb.out.ll<-apply(fb.out.m,1,mort.llvec)
plot(fb.out.ll,type="l",main="40000/20000")

#Calculate goodness of fit: how to calculate?? Per plot/subplot?
fb.pm<-colMeans(fb.out.m)
pred<-pred.mort(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6],fb.pm[7])
plot(pred,mdata$dead,main=paste("r2=",cor(pred,mdata$dead)^2))
abline(0,1)

dev.off()

#Calculate credible intervals
fb.ci<-apply(fb.out.m,2,FUN=quantile,probs=c(0.025,0.5,0.975))
df.fb.ci<-as.data.frame(fb.ci)
write.table(df.fb.ci,"Parameters mort model 1.txt",row.names=F,quote=F,sep="\t")

#  s1 = c(1e-3,10,1,1,1,1),
#  s3 = c(1e-3,100,1,1,1,1),
#  s4 = c(1e-3,1,1,0,1,1),

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