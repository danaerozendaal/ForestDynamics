
#################################
#GROWTH MODEL WITH ALL WD EFFECTS

#Load/select/format data
#dropbox <- "C:/Users/DMAR/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/rozendad/Dropbox/Rainfor data edited"
dropbox <- "C:/Users/vande202/Dropbox/Rainfor data edited"

alldata<-read.table(paste(dropbox,"All data for analysis.txt",sep="/"),h=T)

data<-alldata[,c("PlotCode","subplotID","subplot.area",
                 "subplot.no.trees","TreeID","family",
                 "binomial","region","IntervalLength",
                 "CensusDate","dbh0","dbhgrowth","BA0","dead",
                 "recruit","WD","subplotBA.m2ha")]

gdata<-data[!is.na(data$dbhgrowth) & !is.na(data$dbh0) & !is.na(data$WD) &
              !is.na(data$subplot.area) & !is.na(data$subplot.no.trees) & 
              !is.na(data$dead) & data$dbh0>0,]


#set knots
knots <- c(10,20,30,40,50,60,70,80,90,100)

#read in parameters
pars <- read.csv("D:\\code\\io\\RainforGrowthPars.csv")
pars <- as.matrix(pars[pars$burnin=="sampling",-(1:2)])

#posterior means
pm <- colMeans(pars)

#matrix of plot-level parameters
p.mat <- matrix(nrow=189, ncol=2+2*length(knots))
for (i in 0:(2*length(knots)+1)) {
  p.mat[,i+1] <- pm[grep(paste("p",i,"[0-9]",sep="_"), names(pm))]
}

#create matrix showing how plot parameters vary with WD
wd.seq <- seq(0.3,0.7,by=0.01)

#p.wd is pars x wd x plots array
p.wd <- array(dim=c(length(knots)+1,length(wd.seq),189))
for(i in 1:(length(knots)+1)) {
  p.wd[i,,] <- cbind(1-wd.seq,wd.seq) %*% t(p.mat[,(2*i-1):(2*i)])
}

#PLOT #1: plot each par across range of WD values for each plot
layout(matrix(1:10,nrow=2,byrow=T))
for (i in 1:(length(knots)+1)) {
  matplot(wd.seq,p.wd[i,,],type="l",lty=1,col=rgb(0,0,0,0.2))
  title(main=i)
}


#spline functino for generating predictions from knots (for a particular WD)
spline.pred <- function (dbh,y,endslope) {
  #(end slope not currently used for natural spline)
  
  #essentially 2 different versions here: 1 clamped at right end, 1 natural at rigth end
  #clamped version is currently commented out
  
  len <- length(y)
  
  #breakpoints
  x <- knots
  
  #intermediate variables
  h <- diff(x)
  d <- diff(y)/h
  u <- 6*diff(d)
  
  #clamped at right side
  #u[length(u)] <- u[length(u)] - 3*(endslope-d[length(d)])
  #natural at right side 
  #(leave u as above)
  
  
  #left side matrix
  ls <- matrix(0,nrow=len-2,ncol=len-2)
  #first row
  ls[1,1:2] <- c(2*(h[1]+h[2]),h[2]) 
  #last row (clamped)
  #ls[len-2,len-(3:2)] <- c(h[len-2],2*h[len-2]+1.5*h[len-1]) 
  #last row (natural)
  ls[len-2,len-(3:2)] <- c(h[len-2],2*(h[len-2]+h[len-1])) 
  
  #middle rows
  if (len>4) for (i in 2:(len-3)) {
    ls[i,i+(-1:1)] <- c(h[i-1],2*(h[i-1]+h[i]),h[i])  
  }
  
  m <- c(0,solve(ls) %*% u)
  
  #clamped right side: add last m from Table 5.8 ("Numerical Methods Using Matlab")
  #m[length(m)+1] <- 3*(endslope-d[length(d)])/h[length(h)] - 0.5*m[length(m)]
  #natural
  m[length(m)+1] <- 0
  
  #create regression coefficients
  s <- cbind(y[-len],d-h*(2*m[-len]+m[-1])/6, m[-len]/2, (m[-1]-m[-len])/(6*h))
  
  #clamped spline: add right side curve (same int and slope as  next one)
  #s <- rbind(c(s[1,1:2],0,0),s,c(y[len],endslope,0,0))
  #natural spline
  s <- rbind(c(s[1,1:2],0,0),s,c(y[len],d[len-1]+h[len-1]*(+m[len-1])/6,0,0))
  
  i <- findInterval(dbh,x)+1
  x.obs2 <- dbh - c(x[1],x)[i] #leftmost segment has same reference x as next one
  
  #create preds
  k.pred <- s[i,1] + s[i,2]*x.obs2 + s[i,3]*x.obs2^2 + s[i,4]*x.obs2^3
  
  y.pred <- exp(k.pred)
  
  return(y.pred)
  
} #spline.pred

#PLOT #2: dbh x WD predictions by plot
dbh.seq <- 10:120

#p.preds is wd x dbh x plot array
p.preds <- array(dim=c(3,length(dbh.seq),dim(p.wd)[3]))
for (iWD in 1:3) { #loop across WD values (min, median, max)
  wd.pos <- c(1,round(length(wd.seq)/2,0),length(wd.seq))[iWD]
  
  for (iPlot in 1:dim(p.wd)[3]) { #loop across plots
    p.preds[iWD,,iPlot] <- spline.pred(dbh.seq,p.wd[1:length(knots),wd.pos,iPlot],p.wd[7,wd.pos,iPlot])
    
  } #for iPlot
} #for iWD

#make set of 3 plots
layout(matrix(1:3,nrow=1))
matplot(dbh.seq,p.preds[1,,],type="l",lty=1,col=rgb(0,0,0,0.2),xlim=c(10,120),ylim=c(0,2),log="",xlab="DBH",ylab="Annual Mortality")
lines(dbh.seq,rowMeans(p.preds[1,,]),col="red",lwd=3)
title(main="WD 0.3")
matplot(dbh.seq,p.preds[2,,],type="l",lty=1,col=rgb(0,0,0,0.2),xlim=c(10,120),ylim=c(0,2),log="",xlab="DBH",ylab="Annual Mortality")
lines(dbh.seq,rowMeans(p.preds[2,,]),col="red",lwd=3)
title(main="WD 0.5")
matplot(dbh.seq,p.preds[3,,],type="l",lty=1,col=rgb(0,0,0,0.2),xlim=c(10,120),ylim=c(0,2),log="",xlab="DBH",ylab="Annual Mortality")
lines(dbh.seq,rowMeans(p.preds[3,,]),col="red",lwd=3)
title(main="WD 0.7")

#PLOT #3: predictions one plot at a time (for WD=0.5)
#same plots, but now one plot at a time (layout 12 at a time)
layout(matrix(1:12,nrow=3,byrow=T))
par(mar=c(3,2,1,0))
for(i in 1:189) plot(dbh.seq,p.preds[2,,i],ylim=c(0,2),type="l",main=i,xlab="",ylab="")

#obs vs pred plots
gdata$sp.num <- match(gdata$binomial,unique(gdata$binomial))
gdata$plot.num <- match(gdata$PlotCode,unique(gdata$PlotCode))

#read in predictions from file
gdata$pred <- read.csv(paste(dropbox,"RainforGrowthPars_meanpreds.csv",sep="/"),header=T)[,1]

#group obs and pred by DBH bin
layout(1)
dbin <- 10*round(gdata$dbh0/10,0)
obs.dbin <- tapply(gdata$dbhgrowth,dbin,FUN=mean)
pred.dbin <- tapply(gdata$pred,dbin,FUN=mean)

d.ord <- order(unique(dbin))

#PLOT #4: mean obs vs pred by DBH bin
plot(unique(dbin)[d.ord],obs.dbin,type="l")
lines(unique(dbin)[d.ord],pred.dbin,col="red")

#same idea, now for WD bins
wdbin <- round(gdata$WD,1)
obs.wdbin <- tapply(gdata$dbhgrowth,wdbin,FUN=mean)
pred.wdbin <- tapply(gdata$pred,wdbin,FUN=mean)

wd.ord <- order(unique(wdbin))

#PLOT #5: mean obs vs pred by WD bin
plot(unique(wdbin)[wd.ord],obs.wdbin,type="l")
lines(unique(wdbin)[wd.ord],pred.wdbin,col="red")

#PLOT #6: obs vs pred by plot
plot(tapply(gdata$pred,gdata$PlotCode,FUN="mean"),tapply(gdata$dbhgrowth,gdata$PlotCode,FUN="mean"))
abline(0,1)

##################################################################

