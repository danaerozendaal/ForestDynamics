##########################################################################################
##RECRUITMENT MODEL

library(filzbach)

alldata<-read.table("All data for analysis.txt",h=T)

data<-alldata[alldata$subplot.type!="small",c("PlotCode","subplotID","subplot.area",
                                              "subplot.no.trees","TreeID",
                                              "binomial","region","IntervalLength",
                                              "CensusDate","dbh0","dbhgrowth","BA0","dead",
                                              "recruit","WD","subplotBA.m2ha")]

#4 trees with dbh0=0, now excluded.
rdata<-data[!is.na(data$dbh0) & data$dbh0>0 & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees),]  #check still!

#Calculate recruitment rate  #not needed
rdata2<-aggregate(data$recruit,list(data$PlotCode,data$subplotID,data$subplot.area,
                                    data$subplot.no.trees,data$region,data$IntervalLength,
                                    data$CensusDate,data$subplotBA.m2ha),sum,na.rm=T)
names(rdata2)<-c("PlotCode","subplotID","subplot.area","subplot.no.trees",
                 "region","IntervalLength","CensusDate","subplotBA.m2ha","sum.recruits")

rdata2$subplotID<-factor(rdata2$subplotID)
rdata2$subplotID<-as.numeric(rdata2$subplotID)
rdata2$PlotCode<-factor(rdata2$PlotCode)
rdata2$PlotCode<-as.numeric(rdata2$PlotCode)

rdata3<-rdata2

##INITIAL RECRUITMENT MODEL

#function for predicted recruitment rates
pred.recr<-function(pot.recr,c0,E){
  
  #Competition effect (just subplot BA, focal/neighbour WD can be incorporated)
  r.comp <- exp(-c0 * rdata3$subplotBA.m2ha)
  
  #Predicted annual recruitment
  pred.ann <- pot.recr * r.comp * E
  pred.int <- (pred.ann * rdata3$IntervalLength) * rdata3$subplot.area
  
  return(pred.int)
}

recr.ll<-function(pot.recr,c0,E_all,k_overdisp,E_mean,E_sd){
  
  #predicted number of recruits
  pred.int<-pred.recr(pot.recr,c0,E_all[rdata3$PlotCode])
  
  #calculate likelihood
  ll<-sum(dnbinom(rdata3$sum.recruits, size = k_overdisp, mu = pred.int, log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  
  #sum likelihood per tree
  ll_tot<-ll + log_E_hier
  
  if(is.na(ll_tot)) print(range(pred.int))
  #if(is.na(p.int)) print("NA")
  
  return(ll_tot)
  
}

fb.pars.r <- list(
  pot.recr = c(1e-3,100,1,1,0,1),
  c0 = c(1e-3,5,1,0,0,1),
  E_all = c(1e-3,5,1,0,0,1,181),
  k_overdisp = c(1e-3,10,0,0,0,1),
  E_mean = c(1e-3,5,1,1,0,1),
  E_sd = c(1e-3,5,1,1,0,1)
)

t1<-Sys.time()

fb.out.r<-filzbach(40000,20000,recr.ll,nrow(rdata3),fb.pars.r)
df.fb.out<-as.data.frame(fb.out.r)
write.table(df.fb.out,"FB output recr model 1.txt",row.names=F,quote=F,sep="\t")

Sys.time() - t1

final_out <- paste(readLines("C:/Users/DMAR/Dropbox/Current projects/UofR/ForestDynamics model runs recruitment/workspace/Default_MCMC_final_out.txt"), collapse="\t")
write.table(final_out,"Final_out recr model 1.txt",row.names=F,quote=F,sep="\t")

#Convergence
pdf("Model recr 1.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(5,4,1,1))
recr.llvec<-function(x) recr.ll(x[1],x[2],x[3:183],x[184],x[185],x[186])
fb.out.ll<-apply(fb.out.r,1,recr.llvec)
plot(fb.out.ll,type="l",main="40000/20000")

dev.off()

#Calculate credible intervals
fb.ci<-apply(fb.out.r,2,FUN=quantile,probs=c(0.025,0.5,0.975))
df.fb.ci<-as.data.frame(fb.ci)
write.table(df.fb.ci,"Parameters recr model 1.txt",row.names=F,quote=F,sep="\t")
