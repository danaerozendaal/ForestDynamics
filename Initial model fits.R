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

#Calulate ln(dbh/20)
gdata$log.dbh0<-log(gdata$dbh0/20)

#Put values for all subplots into dataframes in a list
split.gdata <- split(gdata,gdata$subplotID)

min.WD <- min(gdata$WD)
max.WD <- max(gdata$WD)

##########################################################################################
##########################################################################################
##GROWTH MODEL (not working yet, competition dependent on WD)

comp.fun <- function(c2,target) {
  
  #if (!is.finite(min(exp(c2*target)))) print(range(target))
  return(exp(c2 * target))
  
}

subplot.comp <- function(subplot,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD) {
  
  c0_slope <- (c0_maxWD - c0_minWD)/(max.WD - min.WD)
  c0_int <- c0_minWD - (c0_slope * min.WD)
  
  c1_slope <- (c1_maxWD - c1_minWD)/(max.WD - min.WD)
  c1_int <- c1_minWD - (c1_slope * min.WD) 
  
  c2_slope <- (c2_maxWD - c2_minWD)/(max.WD - min.WD)
  c2_int <- c2_minWD - (c2_slope * min.WD)
  
  #competition effect matrix
  c2_WD <- c2_int+c2_slope*subplot$WD  
  comp <- outer(c2_WD,subplot$log.dbh0,FUN=comp.fun)
  
  diag(comp) <- 0
  #print(comp)
  
  #growth effect on focal trees
  c0_WD <- c0_int + c0_slope*subplot$WD
  c1_WD <- c1_int + c1_slope*subplot$WD
  fcomp <- exp(-(c0_WD) * (exp(subplot$log.dbh0*c1_WD)) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  
}

#function for predicted growth rates
pred.growth<-function(pg_minWD,pg_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E){
  
  pg_slope <- (pg_maxWD - pg_minWD)/(max.WD - min.WD)
  pg_int <- pg_minWD - (pg_slope * min.WD)
  
  pot.growth = pg_int + (pg_slope * gdata$dbh0)
  #pot.growth <- 1
  
  s1_slope <- (s1_maxWD - s1_minWD)/(max.WD - min.WD)
  s1_int <- s1_minWD - (s1_slope * min.WD)
  s2_slope <- (s2_maxWD - s2_minWD)/(max.WD - min.WD)
  s2_int <- s2_minWD - (s2_slope * min.WD)
  
  s1_WD <- s1_int + (s1_slope * gdata$WD)
  s2_WD <- s2_int + (s2_slope * gdata$WD)
  g.size = ((gdata$dbh0/202)^s1_WD) * exp(-(s2_WD) * gdata$dbh0)
  #g.size <- 1
  
  g.comp = unlist(lapply(split.gdata,FUN=subplot.comp,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD))
  #g.comp <- 1
  
  pred = pot.growth * E * g.size * g.comp
  
  #print(c(range(pot.growth),range(g.size),range(g.comp)))
  
  return(pred)
  
}

growth.ll <- function(pg_minWD,pg_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E_all,E_mean,E_sd,sigma_int,sigma_slope) {  
  
  g.pred <- pred.growth(pg_minWD,pg_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E_all[gdata$PlotCode])
  
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
  pg_minWD = c(1e-3,100,1,1,0,1),
  pg_maxWD = c(1e-3,100,1,1,0,1),
  s1_minWD = c(1e-6,3,1,0,1,1),  
  s1_maxWD = c(1e-6,3,1,0,1,1), 
  s2_minWD = c(1e-6,10,1,0,1,1),
  s2_maxWD = c(1e-6,10,1,0,1,1),
  c0_minWD = c(-50,50,1,0,0,1),
  c0_maxWD = c(-50,50,1,0,0,1),
  c1_minWD = c(-5,5,0,0,0,1),
  c1_maxWD = c(-5,5,0,0,0,1),
  c2_minWD = c(-5,5,2,0,0,1),
  c2_maxWD = c(-5,5,2,0,0,1),
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
