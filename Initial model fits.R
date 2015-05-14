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
gdata$PlotCode<-factor(gdata$PlotCode)
gdata$PlotCode<-as.numeric(gdata$PlotCode)

#See below: divide by 20 already
gdata$log.dbh0<-log(gdata$dbh0/20)

#Put values for all subplots into a list
fDBH<-split(gdata$dbh0,gdata$subplotID)
WD<-split(gdata$WD,gdata$subplotID)
PlotCode<-split(gdata$PlotCode,gdata$subplotID)
dbhgrowth<-split(gdata$dbhgrowth,gdata$subplotID)
log.fDBH<-split(gdata$log.dbh0,gdata$subplotID)

ul.WD <- unlist(WD)
ul.PlotCode <- unlist(PlotCode)
ul.dbhgrowth <- unlist(dbhgrowth)
ul.fDBH <- unlist(fDBH)
ul.log.fDBH <- unlist(log.fDBH)

##########################################################################################
##########################################################################################
##MODEL 1: potential growth (WD)

comp.fun <- function(focal,target,c2) {
  
  #return((target/20)^c2)
  return(exp(c2 * target))
  
}

subplot.comp <- function(fdbhs,c0,c1,c2) {
  
  #competition effect matrix
  comp <- outer(fdbhs,fdbhs,FUN=comp.fun,c2)
  
  #growth effect on focal trees
  #fcomp <- exp(-c0 * ((fdbhs/20)^c1) * rowSums(comp))
  fcomp <- exp(-c0 * ((exp(fdbhs))^c1) * rowSums(comp))
  
}

#function for predicted growth rates
pred.growth<-function(g1,g2,s1,s2,c0,c1,c2,E){
  
  pot.growth = g1 + g2*ul.WD
  #pot.growth <- 1
  
  g.size = ((ul.fDBH/202)^s1) * exp(-s2 * ul.fDBH)
  #g.size<- 1
  
  #g.comp = unlist(lapply(fDBH,FUN=subplot.comp,c0,c1,c2))
  g.comp = unlist(lapply(log.fDBH,FUN=subplot.comp,c0,c1,c2))
  #g.comp <- 1
  
  pred = pot.growth * E * g.size * g.comp
  
  return(pred)
}

#growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma_exp) {
growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma_int,sigma_slope) {
  #growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma_slope) {
  #growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma) {
  
  g.pred <- pred.growth(g1,g2,s1,s2,c0,c1,c2,E_all[ul.PlotCode])
  
  sigma <- sigma_int + sigma_slope * g.pred
  #sigma <- g.pred^sigma_exp
  #sigma <- sigma_slope * g.pred
  
  #likelihood
  g.ll <- sum(dnorm(ul.dbhgrowth,g.pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  if(is.na(log_E_hier)) print(range(E_all))
  if(is.na(g.ll)) print(range(g.pred))
  
  ll <- sum(g.ll) + sum(log_E_hier)
  #print(ll)
  
  return(ll)
  
}

fb.pars <- list(
  g1 = c(-100,100,1,0,0,1),
  g2 = c(-100,100,1,0,0,1),
  s1 = c(1e-3,100,0,0,1,1),       #s1=1 is a linear relation, s1<1=saturating, s1>1 (kind of exponential?)  
  s2 = c(1e-3,100,1,0,1,1),
  c0 = c(-100,100,1,0,0,1),      
  c1 = c(-100,100,0,0,0,1),
  c2 = c(-100,100,2,0,0,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1),
  #sigma = c(1e-3,10,1,1,0,1)
  sigma_int = c(1e-3,10,1,1,0,1),
  sigma_slope = c(1e-3,10,1,1,0,1)
  #sigma_exp = c(1,10,1,1,0,1)
)

fb.out<-filzbach(20000,5000,growth.ll,nrow(gdata),fb.pars)
df.fb.out<-as.data.frame(fb.out)
write.table(df.fb.out,"FB output model 1.txt",row.names=F,quote=F,sep="\t")

final_out<-read.table(paste("C:/Users/rozendad/Dropbox/Current projects/UofR/ForestDynamics model runs/workspace","Default_MCMC_final_out.txt",sep="/"),h=T)
write.table(final_out,"Final_out model 1.txt",row.names=F,quote=F,sep="\t")

#Convergence
pdf("Model 1.pdf")
par(mfrow=c(1,2),mar=c(5,4,1,1))
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8:188],x[189],x[190],x[191],x[192])
fb.out.ll<-apply(fb.out,1,growth.llvec)
plot(fb.out.ll,type="l",main="40000/20000")

#Calculate goodness of fit
fb.pm<-colMeans(fb.out)
pred<-pred.growth(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6],fb.pm[7],
                  (fb.pm[8:188])[ul.PlotCode])
plot(pred,ul.dbhgrowth,main=paste("r2=",cor(pred,ul.dbhgrowth)^2)
     abline(0,1)
     
     dev.off()
     
     #Calculate credible intervals
     fb.ci<-apply(fb.out,2,FUN=quantile,probs=c(0.025,0.5,0.975))
     df.fb.ci<-as.data.frame(fb.ci)
     write.table(df.fb.ci,"Parameters model 1.txt",row.names=F,quote=F,sep="\t")
     
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

#Initial models for growth, mortality, and recruitment

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
gdata$PlotCode<-factor(gdata$PlotCode)
gdata$PlotCode<-as.numeric(gdata$PlotCode)

meandbh<-mean(gdata$dbh0)  #now 21.6
maxdbh<-max(gdata$dbh0)  #now 202

##########################################################################################
##########################################################################################
##GROWTH MODEL TRIAL WITH SIZE-DEPENDENT COMPETITION EFFECTS

#fix c1 to 0, and c2 to 2. Then the competition effect reduces to total subplot BA.
#RUN ALSO with non-fixed parameters.

library(filzbach)

#Put values for all subplots into a list
fDBH<-split(gdata$dbh0,gdata$subplotID)

WD<-split(gdata$WD,gdata$subplotID)

PlotCode<-split(gdata$PlotCode,gdata$subplotID)

dbhgrowth<-split(gdata$dbhgrowth,gdata$subplotID)

comp.fun <- function(focal,target,c2) {
  
  #   #matrix algebra version
  #   B <- c(c0,c1,c2)
  #   X <- cbind(1,focal,target)
  #   return(plogis(X%*%B)*target)    
  
  return((target/20)^c2)
    
}

#comp.fun <- function(c2,target) {
   
#  return((target/20)^c2)
  
#}

subplot.comp <- function(dbhs,c0,c1,c2) {
  
  #competition effect matrix
  comp <- outer(dbhs,dbhs,FUN=comp.fun,c2)
  
  #growth effect on focal trees
  fcomp <- exp(-c0 * ((dbhs/20)^c1) * rowSums(comp))
  
}

ul.WD <- unlist(WD)
ul.PlotCode <- unlist(PlotCode)
ul.dbhgrowth <- unlist(dbhgrowth)
ul.fDBH <- unlist(fDBH)

#function for predicted growth rates
pred.growth<-function(g1,g2,s1,s2,c0,c1,c2,E){
  
  pot.growth = g1 + g2*ul.WD
  #pot.growth <- 1
  
  g.size = ((ul.fDBH/202)^s1) * exp(-s2 * ul.fDBH)
  #g.size<- 1
  
  g.comp = unlist(lapply(fDBH,FUN=subplot.comp,c0,c1,c2))
  #g.comp <- 1
  
  pred = pot.growth * E * g.size * g.comp
  
  return(pred)
}

growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma_int,sigma_slope) {
#growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,E_all,E_mean,E_sd,sigma) {
  
  g.pred <- pred.growth(g1,g2,s1,s2,c0,c1,c2,E_all[ul.PlotCode])
    
  sigma <- sigma_int + sigma_slope * g.pred
  
  #likelihood
  g.ll <- sum(dnorm(ul.dbhgrowth,g.pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  if(is.na(log_E_hier)) print(range(E_all))
  if(is.na(g.ll)) print(range(g.pred))
  
  ll <- sum(g.ll) + sum(log_E_hier)
  #print(ll)
  
  return(ll)
  
}

fb.pars2 <- list(
  g1 = c(-100,100,1,0,1,1),
  g2 = c(-100,100,1,0,1,1),
  s1 = c(1e-6,100,0,0,0,1),       #s1=1 is a linear relation, s1<1=saturating, s1>1 (kind of exponential?)  
  s2 = c(-100,100,0,0,0,1),
  c0 = c(-100,100,1,0,1,1),      
  c1 = c(-100,100,0,0,1,1),
  c2 = c(-100,100,2,0,1,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1),
  #sigma = c(1e-3,10,1.0,1,0,1)
  sigma_int = c(1e-3,10,1,1,0,1),
  sigma_slope = c(1e-3,10,1,1,0,1)
)

fb.out.g2<-filzbach(40000,20000,growth.ll,nrow(gdata),fb.pars2)
df.fb.out.g2<-as.data.frame(fb.out.g2)
write.table(df.fb.out.g2,"fb.out.g.size-dep-comp fixed.txt",row.names=F,quote=F,sep="\t")

#Convergence
x11(5,5)
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7:187],x[188],x[189],x[190])
fb.out.ll2<-apply(fb.out,1,growth.llvec)
plot(fb.out.ll2,type="l",main="size-dep-comp fixed")

#Calculate goodness of fit
x11(5,5)
fb.pm<-colMeans(fb.out)
pred<-pred.growth(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6],
                  (fb.pm[7:187])[ul.PlotCode])
plot(pred,ul.dbhgrowth,main="size-dep-comp fixed")
abline(0,1)
cor(pred,ul.dbhgrowth)^2

#Calculate credible intervals
fb.ci2<-apply(fb.out,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci2
fb.ci3<-as.data.frame(fb.ci2)
write.table(fb.ci3,"Parameter estimates growth size-dep-comp fixed.txt",row.names=F,quote=F,sep="\t")

##########################################################################################
##########################################################################################
##GROWTH MODEL TRIAL WITH SIZE-DEPENDENT COMPETITION EFFECTS

#ADJUST: ADD WD EFFECTS, all parameters as linear functions of WD?

library(filzbach)

#Put values for all subplots into a list
fDBH<-split(gdata$dbh0,gdata$subplotID)

WD<-split(gdata$WD,gdata$subplotID)

PlotCode<-split(gdata$PlotCode,gdata$subplotID)

dbhgrowth<-split(gdata$dbhgrowth,gdata$subplotID)

ul.WD <- unlist(WD)

comp.fun <- function(focal,target,c2a,c2b) {
  
  #   #matrix algebra version
  #   B <- c(c0,c1,c2)
  #   X <- cbind(1,focal,target)
  #   return(plogis(X%*%B)*target)    
  
  #c2<-c2a + c2b * WD
  return((target/21.6)^c2)
  
}

subplot.comp <- function(dbhs,c0a,c0b,c1a,c1b,c2a,c2b) {
  
  #competition effect matrix
  comp <- outer(dbhs,dbhs,FUN=comp.fun,c2a,c2b)
  
  #growth effect on focal trees
  #c0<-c0a + c0b * WD
  #c1<-c1a + c1b * WD
  fcomp <- exp(-c0 * ((dbhs/21.6)^c1) * rowSums(comp))
  
}

#ul.WD <- unlist(WD)
ul.PlotCode <- unlist(PlotCode)
ul.dbhgrowth <- unlist(dbhgrowth)
ul.fDBH <- unlist(fDBH)

#function for predicted growth rates
pred.growth<-function(g1,g2,s1,c0a,c0b,c1a,c1b,c2a,c2b,E){
  
  pot.growth = g1 + g2*ul.WD
  #pot.growth <- 1
  
  #g.size = (ul.DBH/202)^s1
  g.size<- 1
  
  g.comp = unlist(lapply(fDBH,FUN=subplot.comp,c0a,c0b,c1a,c1b,c2a,c2b))
  #g.comp <- 1
  
  pred = pot.growth * E * g.size * g.comp
  
  return(pred)
}

#growth.ll <- function(g1,g2,s1,s2,c0,c1,c2,c4,E_all,E_mean,E_sd,sigma_int,sigma_slope) {
growth.ll <- function(g1,g2,s1,c0a,c0b,c1a,c1b,c2a,c2b,E_all,E_mean,E_sd,sigma) {
  
  g.pred <- pred.growth(g1,g2,s1,c0a,c0b,c1a,c1b,c2a,c2b,E_all[ul.PlotCode])
  
  #sigma <- sigma_int + sigma_slope * pred
  
  #likelihood
  g.ll <- sum(dnorm(ul.dbhgrowth,g.pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  if(is.na(log_E_hier)) print(range(E_all))
  if(is.na(g.ll)) print(range(g.pred))
  
  ll <- sum(g.ll) + sum(log_E_hier)
  #print(ll)
  
  return(ll)
  
}

#c0 should be positive?
fb.pars2 <- list(
  g1 = c(-100,100,1,0,0,1),
  g2 = c(-100,100,1,0,0,1),
  s1 = c(1e-6,1,0,0,0,1),       #s1=1 is a linear relation, s1<1=saturating, s1>1 (kind of exponential?)  
  c0a = c(-100,100,1,0,0,1),
  c0b = c(-100,100,1,0,0,1),
  c1a = c(-100,100,0,0,0,1),
  c1b = c(-100,100,0,0,0,1),
  c2a = c(-100,100,2,0,0,1),
  c2b = c(-100,100,2,0,0,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1),
  sigma = c(1e-3,10,1.0,1,0,1)
  #sigma_int = c(1e-3,10,1.0,1,0,1),
  #sigma_slope = c(1e-3,10,1.0,1,0,1)
)

fb.out.g2<-filzbach(40000,20000,growth.ll,nrow(gdata),fb.pars2)
df.fb.out.g2<-as.data.frame(fb.out.g2)
write.table(df.fb.out.g2,"fb.out.g.size-dep-comp WD.txt",row.names=F,quote=F,sep="\t")

#Convergence
x11(5,5)
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7:187],x[188],x[189],x[190])
fb.out.ll2<-apply(fb.out,1,growth.llvec)
plot(fb.out.ll2,type="l",main="size-dep-comp WD")

#Calculate goodness of fit
x11(5,5)
fb.pm<-colMeans(fb.out)
pred<-pred.growth(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6],
                  (fb.pm[7:187])[ul.PlotCode])
plot(pred,ul.dbhgrowth,main="size-dep-comp WD")
abline(0,1)
cor(pred,ul.dbhgrowth)^2

#Calculate credible intervals
fb.ci2<-apply(fb.out,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci2
fb.ci3<-as.data.frame(fb.ci2)
write.table(fb.ci3,"Parameter estimates growth size-dep-comp WD.txt",row.names=F,quote=F,sep="\t")

##########################################################################################
##########################################################################################
#Growth model with random effect and subplot BA as competition effect 

library(filzbach)

#function for predicted growth rates
pred.growth<-function(g1,g2,E,s1,s2,c1,c2){
  
  pot.growth = g1 + g2*gdata$WD
  g.size = exp(-0.5*(log(gdata$dbh0/s1)/s2)^2)
  g.comp = c1 + (1-c1) * exp(-c2*gdata$subplotBA.m2ha)
  
  pred = pot.growth * E * g.size * g.comp
  
  return(pred)
}

#log-likelihood function
#growth.ll <- function (g1,g2,E_all,E_mean,E_sd,s1,s2,c1,c2,sigma_int,sigma_slope){
growth.ll <- function (g1,g2,E_all,E_mean,E_sd,s1,s2,c1,c2,sigma){
  
  pred = pred.growth(g1,g2,E_all[gdata$PlotCode],s1,s2,c1,c2)
  
  #sigma <- sigma_int + sigma_slope * pred
  
  #likelihood
  loglike<-sum(dnorm(gdata$dbhgrowth,pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  if(is.na(log_E_hier)) print(range(E_all))
  if(is.na(loglike)) print(range(pred))
  
  return(loglike + log_E_hier)
  
}

fb.pars <- list(
  g1 = c(-100,100,1,0,0,1),
  g2 = c(-100,100,1,0,0,1),
  E_all = c(1e-3,1,0.5,1,0,1,181),
  E_mean = c(1e-3,1,0.5,1,0,1),
  E_sd = c(1e-3,1,0.5,1,0,1),
  s1 = c(1e-3,1000,1.0,1,0,1),
  s2 = c(1e-3,1000,1.0,1,0,1),
  c1 = c(1e-3,1.0,0.5,1,0,1),      #from 0 to 1 only
  c2 = c(1e-3,100.0,1.0,1,0,1),
  sigma = c(1e-3,10,1.0,1,0,1)
  #sigma_int = c(1e-3,10,1.0,1,0,1),
  #sigma_slope = c(1e-3,10,1.0,1,0,1)
)

fb.out<-filzbach(200000,200000,growth.ll,nrow(gdata),fb.pars)
#write.table(df.fb.out.g,"fb.out.g.txt",row.names=F,quote=F,sep="\t")

#Convergence
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3:183],x[184],x[185],x[186],x[187],
                                    x[188],x[189],x[190])
fb.out.ll2<-apply(fb.out,1,growth.llvec)
plot(fb.out.ll2,type="l")

#Calculate goodness of fit
fb.pm<-colMeans(fb.out)
pred<-pred.growth(fb.pm[1],fb.pm[2],(fb.pm[3:183])[gdata$PlotCode],fb.pm[186],
                  fb.pm[187],fb.pm[188],fb.pm[189])
plot(pred,gdata$dbhgrowth)
cor(pred,gdata$dbhgrowth)^2
abline(0,1)

#Calculate credible intervals
fb.ci2<-apply(fb.out,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci2
fb.ci3<-as.data.frame(fb.ci2)
write.table(fb.ci3,"Parameter estimates growth.txt",row.names=F,quote=F,sep="\t")

##########################################################################################
##########################################################################################
##Simple growth model: just a size, competition (subplot BA), and WD effect

library(filzbach)

#function for predicted growth rates
pred.growth<-function(g1,g2,s1,s2,c1,c2){
  
  pot.growth = g1 + g2*gdata$WD
  #pot.growth = 1
  
  #g.size = exp(-0.5*(log(gdata$dbh0/s1)/s2)^2)
  g.size = 1
  
  g.comp = c1 + (1-c1) * exp(-c2*gdata$subplotBA.m2ha)
  
  pred = pot.growth * g.size * g.comp
  
  return(pred)
}

#log-likelihood function
growth.ll <- function (g1,g2,s1,s2,c1,c2,sigma){
  
  pred = pred.growth(g1,g2,s1,s2,c1,c2)  
  
  #likelihood
  loglike<-sum(dnorm(gdata$dbhgrowth,pred,sigma,log=T))
  
  return(loglike)
  
}

fb.pars.g <- list(
  g1 = c(-10,10,1,0,0,1),
  g2 = c(-10,10,1,0,0,1),         
  #s1 = c(1e-3,100.0,1.0,1,0,1),
  s1 = c(1e-3,1000,1.0,1,1,1),
  s2 = c(1e-3,1000,1.0,1,1,1),
  c1 = c(1e-3,1.0,0.5,1,0,1),      #from 0 to 1 only
  c2 = c(1e-3,100.0,1.0,1,0,1),
  sigma = c(1e-3,10,1.0,1,0,1)
)

fb.out.g<-filzbach(60000,20000,growth.ll,nrow(gdata),fb.pars.g)
#Save output?

#assess convergence (no trend)
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7])
fb.out.ll<-apply(fb.out.g,1,growth.llvec)
plot(fb.out.ll,type="l")

#calculate goodness of fit
fb.pm<-colMeans(fb.out.g)
pred<-pred.growth(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6])
plot(pred,gdata$dbhgrowth)
#plot(pred,gdata$dbh0)
cor(pred,gdata$dbhgrowth)^2
abline(0,1)

#calculate credible intervals on parameters
fb.ci<-apply(fb.out.g,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci
fb.ci3<-as.data.frame(fb.ci)
write.table(fb.ci3,"Parameter estimates growth comp only.txt",row.names=F,quote=F,sep="\t")

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##OTHER MODELS IN PROGRESS

library(filzbach)

##MORTALITY
##simple model based on dbh only
##data$dead==1, dead.
#4 trees with dbh0=0, now excluded.

mdata<-data[!is.na(data$dbh0) & data$dbh0>0 & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead),]  #check still!

#function for predicted mortality rates
pred.mort<-function(pm,s1,s3,s4,c1,c2){
  
  pot.mort <- pm
  
  #s2 <- log(99)/(s4 * (1-s3))
  #m.size <- ((mdata$dbh0/10)^s1)/(1 + exp(s2 * (mdata$dbh0 - s3 * s4)))
  m.size <- 1
  
  m.comp <- c1 + (1 - c1) * exp(-c2 * mdata$subplotBA.m2ha)
  #m.comp <- 1
  
  k <- pot.mort * m.size * m.comp
  p.ann <- 1 / (1+exp(-k))
  p.int <- 1 - ((1-p.ann)^mdata$IntervalLength)
  
  p.int[p.int<0.0001] <- 0.0001
  p.int[p.int>0.9999] <- 0.9999
  
  return(p.int)
}

mort.ll<-function(pm,s1,s3,s4,c1,c2){
  
  ll=numeric()
  
  #calculate likelihood for each tree
  for (i in nrow(mdata)){
    
    #probability of mortality
    p.int<-pred.mort(pm,s1,s3,s4,c1,c2)
    
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
  pm = c(1e-3,1000,1,1,0,1),
  s1 = c(1e-3,10,1,1,0,1),
  s3 = c(1e-3,100,1,1,0,1),
  s4 = c(1e-3,1,1,0,0,1),
  c1 = c(1e-3,1.0,1,1,1,1),      #from 0 to 1 only
  c2 = c(1e-3,100,1,1,1,1)
)  

fb.out.m<-filzbach(20000,20000,mort.ll,nrow(mdata),fb.pars.m)

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



