#################################
#GROWTH MODEL WITH ALL WD EFFECTS

library(data.table)
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

#Calculate ln(dbh/20)
gdata$log.dbh0<-log(gdata$dbh0/20)

#uses data.table library to do aggregation of competition effects
gdata.dt <- data.table(gdata)
setkey(gdata.dt,subplotID) #note: this re-sorts gdata by subplotID; doesn't seem to speed things up much

gdata.dt[,subplot.ba := pi*sum(dbh0^2)/40000/subplot.area ,by=subplotID]

#Put values for all subplots into dataframes in a list
split.gdata <- split(gdata,gdata$subplotID)

min.WD <- min(gdata$WD)
max.WD <- max(gdata$WD)
range.wd <- max.WD - min.WD

##########################################################################################
##########################################################################################
##GROWTH MODEL (competition dependent on WD)

subplot.comp <- function(ldbh,wd,c0_int,c0_slope,c1_int,c1_slope,c2_int,c2_slope) {
  
  #competition effect matrix
  c2_WD <- c2_int+c2_slope*wd
  comp <- exp(c2_WD %*% t(ldbh)) #same as outer with "*"; replaces comp.fun (now deleted)
  
  diag(comp) <- 0
  #print(comp)
  
  #growth effect on focal trees
  c0_WD <- c0_int + c0_slope*wd
  c1_WD <- c1_int + c1_slope*wd
  fcomp <- exp(-c0_WD * exp(ldbh*c1_WD) * rowSums(comp))   #log dbh to avoid exponent in comp.fun
  
  #print(c(c2_int,c2_slope))
  #print(c2_WD)
  #print(ldbh)
  #print(comp)
  #readline()
  #if (!is.finite(sum(fcomp))) cat(paste(sum(!is.finite(fcomp)), " bad g.comp\n"))
  
  return(fcomp) #save small amount of time by commenting out
}

#function for predicted growth rates
pred.growth<-function(pg_minWD,pg_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E){
  
  pg_slope <- (pg_maxWD - pg_minWD)/range.wd
  pg_int <- pg_minWD - (pg_slope * min.WD)
  
  pot.growth = pg_int + (pg_slope * gdata.dt$WD)
  #pot.growth <- 1
  
  s1_slope <- (s1_maxWD - s1_minWD)/range.wd
  s1_int <- s1_minWD - (s1_slope * min.WD)
  s2_slope <- (s2_maxWD - s2_minWD)/range.wd
  s2_int <- s2_minWD - (s2_slope * min.WD)
    
  s1_WD <- s1_int + (s1_slope * gdata.dt$WD)
  s2_WD <- s2_int + (s2_slope * gdata.dt$WD)
  g.size = ((gdata.dt$dbh0/202)^s1_WD) * exp(-s2_WD * gdata.dt$dbh0)
  #g.size <- 1
  
  c0_slope <- (c0_maxWD - c0_minWD)/range.wd
  c0_int <- c0_minWD - (c0_slope * min.WD)
  
  c1_slope <- (c1_maxWD - c1_minWD)/range.wd
  c1_int <- c1_minWD - (c1_slope * min.WD) 
  
  c2_slope <- (c2_maxWD - c2_minWD)/range.wd
  c2_int <- c2_minWD - (c2_slope * min.WD)
      
  #g.comp = unlist(lapply(split.gdata,FUN=subplot.comp,c0_int,c0_slope,c1_int,c1_slope,c2_int,c2_slope))
  #g.comp <- gdata.dt[,subplot.comp(log.dbh0,WD,c0_int,c0_slope,c1_int,c1_slope,c2_int,c2_slope),by=subplotID][,V1]

  g.comp <- exp(-c0_minWD * exp(gdata.dt$log.dbh0*c1_int) * gdata.dt$subplot.ba)

  #g.comp <- 1
  
  if (!is.finite(sum(g.comp))) print(c(c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD))
  
  pred = pot.growth * E * g.size * g.comp
  
  if (!is.finite(sum(pred))) print(c(range(pot.growth),range(g.size),range(g.comp)))
  
  return(pred)
  
}

growth.ll <- function(pg_minWD,pg_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E_all,E_mean,E_sd,sigma_int,sigma_slope) {  
  
  g.pred <- pred.growth(pg_minWD,pg_maxWD,s1_minWD,s1_maxWD,s2_minWD,s2_maxWD,c0_minWD,c0_maxWD,c1_minWD,c1_maxWD,c2_minWD,c2_maxWD,E_all[gdata$PlotCode])
  
  sigma <- sigma_int + sigma_slope * g.pred
  
  #likelihood
  g.ll <- sum(dnorm(gdata.dt$dbhgrowth,g.pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  if(is.na(g.ll)) print(range(g.pred))
  
  #print((cor(g.pred,gdata$dbhgrowth))^2)
  
  ll <- sum(g.ll) + sum(log_E_hier)
  
  #print(ll)
  
  return(ll)
  
}

fb.pars <- list(
  pg_minWD = c(1e-3,100,1,1,0,1),
  pg_maxWD = c(1e-3,100,1,1,0,1),
  s1_minWD = c(1e-6,15,1,1,0,1),  
  s1_maxWD = c(1e-6,15,1,1,0,1), 
  s2_minWD = c(1e-6,15,1,1,0,1),
  s2_maxWD = c(1e-6,15,1,1,0,1),
  c0_minWD = c(1e-6,5,1,1,0,1),
  c0_maxWD = c(1e-6,5,1,1,0,1),
  c1_minWD = c(-3,3,0,0,0,1),
  c1_maxWD = c(-3,3,0,0,0,1),
  c2_minWD = c(-2,3,2,0,0,1),
  c2_maxWD = c(-2,3,2,0,0,1),
  E_all = c(1e-3,1,1,0,1,0,181),
  E_mean = c(1e-3,1,1,1,1,1),
  E_sd = c(1e-3,1,1,1,1,1),
  sigma_int = c(1e-3,10,1,1,0,1),
  sigma_slope = c(1e-3,10,1,1,0,1)
)

system.time(for (i in 1:10) growth.ll(.69723,.71524,.18028,.97117,9.119e-5,3.4242e-3,
                                      1e-5,1e-5,1,1,1,1,
                                      rep(1,181),1,1,.06815,.99879))


t1<-Sys.time()

fb.out<-filzbach(50000,20000,growth.ll,nrow(gdata),fb.pars)
df.fb.out<-as.data.frame(fb.out)
write.table(df.fb.out,"FB output model 11.txt",row.names=F,quote=F,sep="\t")

Sys.time() - t1

final_out <- paste(readLines("C:/Users/DMAR/Dropbox/Current projects/UofR/ForestDynamics model runs 2/workspace/Default_MCMC_final_out.txt"), collapse="\t")
write.table(final_out,"Final_out model 11.txt",row.names=F,quote=F,sep="\t")

#Convergence
pdf("Model 11.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(5,4,1,1))
growth.llvec<-function(x) growth.ll(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],
                                    x[11],x[12],x[13:193],x[194],x[195],x[196],x[197])
fb.out.ll<-apply(fb.out,1,growth.llvec)
plot(fb.out.ll,type="l",main="80000/20000")

#Calculate goodness of fit
fb.pm<-colMeans(fb.out)
pred<-pred.growth(fb.pm[1],fb.pm[2],fb.pm[3],fb.pm[4],fb.pm[5],fb.pm[6],fb.pm[7],
                  fb.pm[8],fb.pm[9],fb.pm[10],fb.pm[11],fb.pm[12],(fb.pm[13:193])[gdata$PlotCode])
plot(pred,gdata$dbhgrowth,main=paste("r2=",cor(pred,gdata$dbhgrowth)^2))
abline(0,1)

dev.off()

#Calculate credible intervals
fb.ci<-apply(fb.out,2,FUN=quantile,probs=c(0.025,0.5,0.975))
df.fb.ci<-as.data.frame(fb.ci)
write.table(df.fb.ci,"Parameters model 11.txt",row.names=F,quote=F,sep="\t")
