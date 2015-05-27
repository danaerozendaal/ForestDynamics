##############
#GROWTH MODELS

library(filzbach)

#Load/select/format data
alldata<-read.table("All data for analysis.txt",h=T)

data<-alldata[alldata$subplot.type!="small",c("PlotCode","subplotID","subplot.area",
                                              "subplot.no.trees","TreeID",
                                              "binomial","region","IntervalLength",
                                              "CensusDate","dbh0","dbhgrowth","BA0","dead",
                                              "recruit","WD","subplotBA.m2ha")]

gdata<-data[!is.na(data$dbhgrowth) & !is.na(data$dbh0) & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead) & data$dbh0>0,]  
gdata$subplotID<-factor(gdata$subplotID)
gdata$subplotID<-as.numeric(gdata$subplotID)
gdata$PlotCode<-factor(gdata$PlotCode)
gdata$PlotCode<-as.numeric(gdata$PlotCode)

#See below: divide by 20 already
gdata$log.dbh0<-log(gdata$dbh0/20)

#Put values for all subplots into a list
split.gdata <- split(gdata,gdata$subplotID)

fDBH<-split(gdata$dbh0,gdata$subplotID)
fWD<-split(gdata$WD,gdata$subplotID)
PlotCode<-split(gdata$PlotCode,gdata$subplotID)
dbhgrowth<-split(gdata$dbhgrowth,gdata$subplotID)
log.fDBH<-split(gdata$log.dbh0,gdata$subplotID)

#ul.fWD <- unlist(fWD)
#ul.PlotCode <- unlist(PlotCode)
#ul.dbhgrowth <- unlist(dbhgrowth)
#ul.fDBH <- unlist(fDBH)
#ul.log.fDBH <- unlist(log.fDBH)

min.WD <- min(gdata$WD)
max.WD <- max(gdata$WD)

#g0 = pot growth at minimum WD
#g1 = pot growth at maximum WD

##########################################################################################
##########################################################################################
##GROWTH MODEL (not working yet, competition dependent on WD)

#c0, c1 also as a function of WD
#s1 (just positive?), s2 (only positive) as a function of WD

comp.fun <- function(target,c2) {
  
  if (!is.finite(min(exp(c2*target)))) print(range(target))

  #return((target/20)^c2)
  return(exp(c2 * target))
  
}

subplot.comp <- function(dbhs,c0,c1,c2) {
#subplot.comp <- function(subplot,c0,c1,c2_int,c2_slope) {
#subplot.comp <- function(dbhs,c0,c1,c2_int,c2_slope) {
    
  #competition effect matrix
  comp <- outer(dbhs,rep(c2,length(dbhs)),FUN=comp.fun)
  print(length(dbhs))
  #comp <- outer(subplot$log.dbh0,c2_int+c2_slope*subplot$WD,FUN=comp.fun)
  #comp <- outer(subplot$log.dbh0,c2_int+c2_slope*(lapply(split.gdata,"[[","WD")),FUN=comp.fun)
  #comp <- outer(dbhs,c2_int+c2_slope*fWD,FUN=comp.fun)
  
  #growth effect on focal trees
  #fcomp <- exp(-c0 * ((dbhs/20)^c1) * rowSums(comp))
  #fcomp <- exp(-(c0_int + c0_slope*ul.fWD) * ((exp(dbhs))^(c1_int + c1_slope*ul.fWD)) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  #fcomp <- exp(-c0 * ((exp(subplot$log.dbh0))^c1) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  fcomp <- exp(-c0 * ((exp(dbhs))^c1) * rowSums(comp[,-1]))   #remove first column (is focal tree, I think?)
    
}

#function for predicted growth rates
pred.growth<-function(g0,g1,s1_0,s1_1,s2_0,s2_1,c0,c1,c2,E){
#pred.growth<-function(g0,g1,s1,s2,c0,c1,c2_int,c2_slope,E){
#pred.growth<-function(g0,g1,s1,s2,c0,c1,c2,E){
  
  WD_slope <- (g1 - g0)/(max.WD - min.WD)
  WD_int <- g0 - (WD_slope * min.WD)
  
  pot.growth = WD_int + (WD_slope * gdata$dbh)
  #pot.growth <- 1
  
  s1_slope <- (s1_1 - s1_0)/(max.WD - min.WD)
  s1_int <- s1_0 - (s1_slope * min.WD)
  s2_slope <- (s2_1 - s2_0)/(max.WD - min.WD)
  s2_int <- s2_0 - (s2_slope * min.WD)
  
  g.size = ((gdata$dbh0/202)^(s1_int + (s1_slope * gdata$WD))) * exp(-(s2_int + (s2_slope * gdata$WD)) * gdata$dbh0)
  #g.size = ((gdata$dbh0/202)^s1) * exp(-s2 * gdata$dbh0)
  #g.size <- 1
  
  #g.comp = unlist(lapply(split.gdata,FUN=subplot.comp,c0,c1,c2_int,c2_slope))
  g.comp = unlist(lapply(fDBH,FUN=subplot.comp,c0,c1,c2))
  #g.comp <- 1
  
  pred = pot.growth * E * g.size * g.comp
  
  #print(c(range(pot.growth),range(g.comp)))
  print(c(range(pot.growth),range(g.size),range(g.comp)))
  #print(c(range(pot.growth),g0,g1,WD_slope,WD_int))
  
  return(pred)
}

#growth.ll <- function(g0,g1,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma_int,sigma_slope) {
growth.ll <- function(g0,g1,s1_0,s1_1,s2_0,s2_1,c0,c1,c2,E_all,E_mean,E_sd,sigma_int,sigma_slope) {
#growth.ll <- function(g0,g1,s1_0,s1_1,s2_0,s2_1,c0,c1,c2_int,c2_slope,E_all,E_mean,E_sd,sigma_int,sigma_slope) {
  
  #g.pred <- pred.growth(g0,g1,s1_0,s1_1,s2_0,s2_1,c0,c2_int,c2_slope,E_all[gdata$PlotCode])
  g.pred <- pred.growth(g0,g1,s1_0,s1_1,s2_0,s2_1,c0,c1,c2,E_all[gdata$PlotCode])
  #g.pred <- pred.growth(g0,g1,s1,s2,c0,c1,c2,E_all[gdata$PlotCode])
  
  sigma <- sigma_int + sigma_slope * g.pred
  
  #likelihood
  g.ll <- sum(dnorm(gdata$dbhgrowth,g.pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  if(is.na(g.ll)) print(range(g.pred))
  
  print((cor(g.pred,gdata$dbhgrowth))^2)
  
  ll <- sum(g.ll) + sum(log_E_hier)
  
  #if(iter %% 100 == 0) print((cor(g.pred,gdata$dbhgrowth))^2) 
  
  return(ll)
  
}

fb.pars <- list(
  g0 = c(1e-3,100,1,1,0,1),
  g1 = c(1e-3,100,1,1,0,1),
  s1_0 = c(1e-6,3,1,0,1,1),  
  s1_1 = c(1e-6,3,1,0,1,1), 
  s2_0 = c(1e-6,10,1,0,1,1),
  s2_1 = c(1e-6,10,1,0,1,1),
  #s1 = c(1e-6,3,1,0,1,1), 
  #s2 = c(1e-6,10,1,0,1,1),
  c0 = c(-50,50,1,0,0,1),      
  c1 = c(-5,5,0,0,0,1),
  c2 = c(-5,5,2,0,0,1),
  #c2_int = c(-5,5,2,0,0,1),
  #c2_slope = c(-5,5,2,0,0,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1),
  sigma_int = c(1e-3,10,1,1,0,1),
  sigma_slope = c(1e-3,10,1,1,0,1)
)

fb.out<-filzbach(40000,20000,growth.ll,nrow(gdata),fb.pars)
df.fb.out<-as.data.frame(fb.out)
write.table(df.fb.out,"FB output model 1.txt",row.names=F,quote=F,sep="\t")

final_out <- paste(readLines("C:/Users/rozendad/Dropbox/Current projects/UofR/ForestDynamics model runs/workspace/Default_MCMC_final_out.txt"), collapse="\t")
write.table(final_out,"Final_out model 1.txt",row.names=F,quote=F,sep="\t")

#Convergence
pdf("Model 1.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(5,4,1,1))
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8:188],x[189],x[190],x[191],x[192])
fb.out.ll<-apply(fb.out,1,growth.llvec)
plot(fb.out.ll,type="l",main="40000/20000")

#Calculate goodness of fit
fb.pm<-colMeans(fb.out)
pred<-pred.growth(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6],fb.pm[7],
                  (fb.pm[8:188])[ul.PlotCode])
plot(pred,ul.dbhgrowth,main=paste("r2=",cor(pred,ul.dbhgrowth)^2))
abline(0,1)

dev.off()

#Calculate credible intervals
fb.ci<-apply(fb.out,2,FUN=quantile,probs=c(0.025,0.5,0.975))
df.fb.ci<-as.data.frame(fb.ci)
write.table(df.fb.ci,"Parameters model 1.txt",row.names=F,quote=F,sep="\t")

##############
#########TRIAL

#lapply(list2,`[`,-1)
#lapply(split.gdata, "[[", "WD")

#comp.fun <- function(target,c2) {
comp.fun <- function(target,c2) {
  
  print(target)
  return(exp(c2 * target))
  
}

subplot.comp <- function(dbhs,c0,c1,c2) {
  
  #competition effect matrix
  #comp <- outer(dbhs,rep(c2,length(dbhs)),FUN=comp.fun)
  comp <- outer(dbhs,rep(c2,length(dbhs)),FUN=comp.fun)
  #print(length(dbhs))
  print(comp)
  print(class(comp))
  #comp <- outer(subplot$log.dbh0,c2_int+c2_slope*subplot$WD,FUN=comp.fun)
  #comp <- outer(subplot$log.dbh0,c2_int+c2_slope*(lapply(split.gdata,"[[","WD")),FUN=comp.fun)
  
  #growth effect on focal trees
  #fcomp <- exp(-(c0_int + c0_slope*ul.fWD) * ((exp(dbhs))^(c1_int + c1_slope*ul.fWD)) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  #fcomp <- exp(-c0 * ((exp(subplot$log.dbh0))^c1) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  #fcomp <- exp(-c0 * ((exp(dbhs))^c1) * rowSums(comp[,-1]))   #remove first column (is focal tree?)
  fcomp <- exp(-c0 * ((exp(dbhs))^c1) * rowSums(comp))   #remove first column (is focal tree?)
  
  print(rowSums(comp))
  print(fcomp)
  return(fcomp)
  
}

dbhs<-list(p1=c(-0.5,-0.3,-0.01,1),p2=c(-0.2,0.2,-0.5),p3=c(1,-2,-0.5))
dbhs2<-cbind(p1=c(-0.5,-0.3,-0.01),p2=c(-0.2,0.2,-0.5),p3=c(1,-2,-0.5))
c0<-0.5
c1<-0
c2<-2

g.comp <- unlist(lapply(dbhs,FUN=subplot.comp,c0,c1,c2))

exp(c2 * -0.5)


