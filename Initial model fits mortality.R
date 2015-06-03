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

##########################################################################################
##########################################################################################
##INITIAL MORTALITY MODEL

library(filzbach)

comp.fun <- function(target,c2) {
  
  if (!is.finite(min(exp(c2*target)))) print(range(target))
  
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

#function for predicted mortality rates
pred.mort<-function(pl_minWD,pl_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E){
  
  #Potential longevity
  pl_slope <- (pl_maxWD - pl_minWD)/(max.WD - min.WD)
  pl_int <- pl_minWD - (pl_slope * min.WD)
  
  pot.longevity <- pl_int + (pl_slope * mdata$WD)
  #pot.longevity <- 1
  
  #Size effect
  s1_slope <- (s1_maxWD - s1_minWD)/(max.WD - min.WD)
  s1_int <- s1_minWD - (s1_slope * min.WD)
  s2_slope <- (s2_maxWD - s2_minWD)/(max.WD - min.WD)
  s2_int <- s2_minWD - (s2_slope * min.WD)
  
  s1_WD <- s1_int + (s1_slope * mdata$WD)
  s2_WD <- s2_int + (s2_slope * mdata$WD)
  m.size = ((mdata$dbh0/202)^s1_WD) * exp(-(s2_WD) * mdata$dbh0)
  #m.size <- 1
  
  #Competition effect
  m.comp = unlist(lapply(split.mdata,FUN=subplot.comp,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD))
  #m.comp <- 1
  
  #k <- pot.mort * m.size * m.comp
  #p.ann <- 1 / (1+exp(-k))
  p.ann <- 1 / (1 + (pot.longevity * m.size * m.comp * E)) 
  p.int <- 1 - ((1-p.ann)^mdata$IntervalLength)
  
  p.int[p.int<0.0001] <- 0.0001
  p.int[p.int>0.9999] <- 0.9999
  
  return(p.int)
}

mort.ll<-function(pl_minWD,pl_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E_all,E_mean,E_sd){
  
  #probability of mortality
  p.int<-pred.mort(pl_minWD,pl_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E_all[mdata$PlotCode])
  
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
  pl_minWD = c(1e-3,1000,1,1,0,1),
  pl_maxWD = c(1e-3,1000,1,1,0,1),
  s1_minWD = c(1e-6,10,1,0,0,1),  
  s1_maxWD = c(1e-6,10,1,0,0,1), 
  s2_minWD = c(1e-6,10,1,0,0,1),
  s2_maxWD = c(1e-6,10,1,0,0,1),
  c0_minWD = c(-10,10,1,0,0,1),
  c0_maxWD = c(-10,10,1,0,0,1),
  c1_minWD = c(-5,5,0,0,0,1),
  c1_maxWD = c(-5,5,0,0,0,1),
  c2_minWD = c(-5,5,2,0,0,1),
  c2_maxWD = c(-5,5,2,0,0,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1)
)

fb.out.m<-filzbach(40000,20000,mort.ll,nrow(mdata),fb.pars.m)
df.fb.out<-as.data.frame(fb.out.m)
write.table(df.fb.out,"FB output mort model 1.txt",row.names=F,quote=F,sep="\t")

final_out <- paste(readLines("C:/Users/DMAR/Dropbox/Current projects/UofR/ForestDynamics model runs/workspace/Default_MCMC_final_out.txt"), collapse="\t")
write.table(final_out,"Final_out mort model 1.txt",row.names=F,quote=F,sep="\t")

#Convergence
pdf("Model mort 1.pdf",width=8,height=4)
mort.llvec<-function(x) mort.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7])
fb.out.ll<-apply(fb.out.m,1,mort.llvec)
plot(fb.out.ll,type="l",main="40000/20000")

dev.off()

#Calculate credible intervals
fb.ci<-apply(fb.out.m,2,FUN=quantile,probs=c(0.025,0.5,0.975))
df.fb.ci<-as.data.frame(fb.ci)
write.table(df.fb.ci,"Parameters mort model 1.txt",row.names=F,quote=F,sep="\t")
