#################################
#MORT MODEL WITH ALL WD EFFECTS

#Load/select/format data
dropbox <- "C:/Users/DMAR/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/rozendad/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/vande202/Dropbox/Rainfor data edited"

alldata<-read.table(paste(dropbox,"All data for analysis.txt",sep="/"),h=T)

data<-alldata[,c("PlotCode","subplotID","subplot.area",
                 "subplot.no.trees","TreeID","family",
                 "binomial","region","IntervalLength",
                 "CensusDate","dbh0","dbhgrowth","BA0","dead",
                 "recruit","WD","subplotBA.m2ha")]

mdata<-data[!is.na(data$dbh0) & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead) & data$dbh0>0,]

#pars <- read.csv("D:\\code\\io\\RainforMortPars.csv")
#pars <- read.csv("D:\\code\\io\\RainforGrowthPars.csv")
pars <- read.csv(paste(dropbox,"RainforMortPars.csv",sep="/"),h=T)

pars <- as.matrix(pars[pars$burnin=="sampling",-(1:2)])

pm <- colMeans(pars)

p.mat <- matrix(nrow=189, ncol=10)
for (i in 0:9) {
  p.mat[,i+1] <- pm[grep(paste("p",i,"[0-9]",sep="_"), names(pm))]
}


wd.seq <- seq(0.3,0.7,by=0.01)

#p.wd is pars x wd x plots array
p.wd <- array(dim=c(5,length(wd.seq),189))
for(i in 1:5) {
  p.wd[i,,] <- cbind(1-wd.seq,wd.seq) %*% t(p.mat[,(2*i-1):(2*i)])
}

#plot each par across range of WD values for each plot
layout(matrix(1:6,nrow=2,byrow=T))
for (i in 1:5) {
  matplot(wd.seq,p.wd[i,,],type="l",lty=1,col=rgb(0,0,0,0.2))
  title(main=i)
}

dbh.seq <- 10:100

spline.pred <- function (dbh,y,endslope) {
  
  #breakpoints
  x <- c(10,20,40,70)    
  
  #intermediate variables
  h <- diff(x)
  d <- diff(y)/h
  u <- 6*diff(d)
  
  #clamped at right side
  u[length(u)] <- u[length(u)] - 3*(endslope-d[length(d)])
  
  #left side matrix constructed for 4 knots, need to change if more
  ls <- matrix(c(
    2*(h[1]+h[2]),h[2],
    h[2],2*h[2]+1.5*h[3]
  ),nrow=2,byrow=T)
  
  m <- c(0,solve(ls) %*% u)
  
  #add last m from Table 5.8 ("Numerical Methods Using Matlab")
  m[length(m)+1] <- 3*(endslope-d[length(d)])/h[length(h)] - 0.5*m[length(m)]
  
  len <- length(y)
  
  #create regression coefficients
  s <- cbind(y[-len],d-h*(2*m[-len]+m[-1])/6, m[-len]/2, (m[-1]-m[-len])/(6*h))
  
  #add left side curve (same int and slope as  next one)
  s <- rbind(c(s[1,1:2],0,0),s,c(y[len],endslope,0,0))
  
  i <- findInterval(dbh,x)+1
  x.obs2 <- dbh - c(x[1],x)[i] #leftmost segment has same reference x as next one
  
  #create preds
  k.pred <- s[i,1] + s[i,2]*x.obs2 + s[i,3]*x.obs2^2 + s[i,4]*x.obs2^3
  
  y.pred <- 1/(1+exp(-k.pred))
  
  return(y.pred)
  
} #spline.pred

#p.preds is wd x dbh x plot array
p.preds <- array(dim=c(3,length(dbh.seq),dim(p.wd)[3]))

for (iWD in 1:3) { #loop across WD values (min, median, max)
  
  wd.pos <- c(1,round(length(wd.seq)/2,0),length(wd.seq))[iWD]
  
  for (iPlot in 1:dim(p.wd)[3]) { #loop across plots
    p.preds[iWD,,iPlot] <- spline.pred(dbh.seq,p.wd[1:4,wd.pos,iPlot],p.wd[5,wd.pos,iPlot])
    
  } #for iPlot
} #for iWD

p.preds.logit <- log(p.preds/(1-p.preds))

layout(matrix(1:3,nrow=1))
matplot(dbh.seq,p.preds.logit[1,,],type="l",lty=1,col=rgb(0,0,0,0.2),xlim=c(10,100),ylim=c(-10,0),log="",xlab="DBH",ylab="Annual Mortality")
lines(dbh.seq,rowMeans(p.preds[1,,]),col="red",lwd=3)
title(main="WD 0.3")
matplot(dbh.seq,p.preds.logit[2,,],type="l",lty=1,col=rgb(0,0,0,0.2),xlim=c(10,100),ylim=c(-10,0),log="",xlab="DBH",ylab="Annual Mortality")
lines(dbh.seq,rowMeans(p.preds[2,,]),col="red",lwd=3)
title(main="WD 0.5")
matplot(dbh.seq,p.preds.logit[3,,],type="l",lty=1,col=rgb(0,0,0,0.2),xlim=c(10,100),ylim=c(-10,0),log="",xlab="DBH",ylab="Annual Mortality")
lines(dbh.seq,rowMeans(p.preds[3,,]),col="red",lwd=3)
title(main="WD 0.7")


#plot observed vs predicted accross dbh bins

mdata$sp.num <- match(mdata$binomial,unique(mdata$binomial))
mdata$plot.num <- match(mdata$PlotCode,unique(mdata$PlotCode))

mdata$pred <- spline.pred(mdata$dbh0,mdata$g.const.term,mdata$dbh.term,mdata$ldbh.term)

mdata$pred.census <- 1 - (1 - mdata$pred)^mdata$IntervalLength

layout(1)

dbin <- 10*round(mdata$dbh0/10,0)
obs.dbin <- tapply(mdata$dead,dbin,FUN=mean)
pred.dbin <- tapply(mdata$pred.census,dbin,FUN=mean)

d.ord <- order(unique(dbin))

plot(unique(dbin)[d.ord],obs.dbin,type="l")
lines(unique(dbin)[d.ord],pred.dbin,col=rgb(1,0,0,0.5))

##################################################################

