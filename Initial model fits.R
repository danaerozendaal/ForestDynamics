#Initial models for growth, mortality, and recruitment

#############
#############
#GROWTH MODEL

library(filzbach)

#function for predicted growth rates
pred.growth<-function(pg1,pg2,E,s1,s2,a,b,c1,c2){
  
  pot.growth = pg1 + pg2*gdata$WD
  g.size = exp(-0.5*(log(gdata$dbh0/s1)/s2)^2)
  c.tree = exp(a + b*gdata$dbh0)/((exp(a + b*gdata$dbh0))+1)
  sel.subplotID = gdata$subplotID
  BA.adj = (sum(c.tree * gdata[gdata$subplotID==sel.subplotID,]$BA0))/
                gdata$subplot.area
  g.comp = c1 + (1-c1) * exp(-c2*BA.adj)
    
  pred = pot.growth * E * g.size * g.comp
  
  return(pred)
}

#likelihood function
growth.ll <- function (pg1,pg2,E_all,E_mean,E_sd,s1,s2,a,b,c1,c2,sigma) {
  
  pred = pred.growth(pg1,pg2,E_all[gdata$PlotCode],s1,s2,a,b,c1,c2)
  
  #likelihood
  loglike<-sum(dnorm(gdata$dbhgrowth.yr,pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  
  return(loglike + log_E_hier)
  
}

fb.pars <- list(
  pg1 = c(1e-3,15,1.0,1,0,1),
  pg2 = c(-10,10,0,0,1,1),         #now fixed at 0
  E_all = c(1e-3,1,0.5,1,0,1,183),
  E_mean = c(1e-3,1,0.5,1,0,1),
  E_sd = c(1e-3,10,0.5,1,0,1),
  s1 = c(1e-3,100.0,1.0,1,0,1),
  s2 = c(1e-3,1000,1.0,1,0,1),
  a = c(-10,10,0,0,1,1),           #now fixed at 0
  b = c(-10,10,0,0,1,1),           #now fixed at 0
  c1 = c(1e-3,1.0,0.5,1,0,1),      #from 0 to 1 only
  c2 = c(1e-3,100.0,1.0,1,0,1),
  sigma = c(1e-3,10,1.0,1,0,1)
)

filzbach(100000,100000,growth.ll,nrow(gdata),fb.pars)

###########################
###########################
#GROWTH MODEL with just a size effect and random plot effect

library(filzbach)

#function for predicted growth rates
pred.growth<-function(pg,E,s1,s2){
  
  g.size = exp(-0.5*(log(gdata$dbh0/s1)/s2)^2)
    
  pred = pg * E * g.size
  
  return(pred)
}

#log-likelihood function
growth.ll <- function (pg,E_all,E_mean,E_sd,s1,s2,sigma){
  
  pred = pred.growth(pg,E_all[gdata$PlotCode],s1,s2)  
  
  #likelihood
  loglike<-sum(dnorm(gdata$dbhgrowth.yr,pred,sigma,log=T))
  
  #parameter hierarchy
  log_E_hier<-sum(dnorm(E_all,E_mean,E_sd,log=T))
  
  return(loglike + log_E_hier)
  
}

fb.pars <- list(
  pg = c(1e-3,15,1.0,1,0,1),
  E_all = c(1e-3,1,0.5,1,0,1,183),
  E_mean = c(1e-3,1,0.5,1,0,1),
  E_sd = c(1e-3,10,0.5,1,0,1),
  s1 = c(1e-3,100.0,1.0,1,0,1),
  s2 = c(1e-3,1000,1.0,1,0,1),
  sigma = c(1e-3,10,1.0,1,0,1)
)

filzbach(100000,100000,growth.ll,nrow(gdata),fb.pars)

