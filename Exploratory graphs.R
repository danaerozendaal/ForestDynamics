##Exploratory graphs
alldata<-read.table("All data for analysis.txt",h=T)

alldata$select<-ifelse((alldata$area.sub.T1==0.04 | alldata$area.sub.st==0.04),0,
                       ifelse(!(is.na(alldata$perc.trees.x.std) | alldata$perc.trees.x.std==0),1,
                              ifelse(!is.na(alldata$perc.trees.x.),2,3)))

##Plots xy-coordinates (xy_field, xy_standard)    ###PROBLEMS WITH OUTLIERS
xydata<-alldata[!is.na(alldata$perc.trees.x) & alldata$perc.trees.x>=90 
                & !is.na(alldata$x) & !is.na(alldata$y),]
xydata.std<-alldata[!is.na(alldata$perc.trees.x.std) & alldata$perc.trees.x.std>=90 
                    & !is.na(alldata$x_standard) & !is.na(alldata$y_standard),]

pdf("Scatter plots xy coordinates.pdf", width = 10, height = 8)

par(mfrow=c(2,2))

for(plot in unique(xydata$PlotCode)){
  
  data<-xydata[xydata$PlotCode==plot,]
  plot(data$x,data$y,xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y),max(data$y)),
        main=paste(unique(data$region)," PlotCode= ",plot," area(ha)= ",round(mean(data$PlotArea.use),2),sep=""))
  
} 

for(plot in unique(xydata.std$PlotCode)){

    data<-xydata.std[xydata.std$PlotCode==plot,]
    plot(data$x_standard,data$y_standard,xlim=c(min(data$x_standard),max(data$x_standard)),
         ylim=c(min(data$y_standard),max(data$y_standard)),
         main=paste(unique(data$region)," PlotCode= ",plot," area(ha)= ", round(mean(data$PlotArea.use),2),sep=""))
    
} 

dev.off()

##SAME, but only for plots without suitable subplots (including any number of trees)
x.std.data<-alldata[(alldata$select==1 | alldata$select==2) & !is.na(alldata$x_standard) & !is.na(alldata$y_standard),]
x.data<-alldata[(alldata$select==1 | alldata$select==2) & !is.na(alldata$x) & !is.na(alldata$y),]

pdf("Scatter plots xy coordinates 2.pdf", width = 10, height = 8)

par(mfrow=c(2,2))

for(plot in unique(x.std.data$PlotCode)){
  
  data<-x.std.data[x.std.data$PlotCode==plot,]
  plot(data$x_standard,data$y_standard,xlim=c(min(data$x_standard),max(data$x_standard)),
       ylim=c(min(data$y_standard),max(data$y_standard)),
       main=paste(unique(data$region)," ",plot," tot.area=",round(mean(data$PlotArea.use),2)," ha",
                    " sub.fld=",round(mean(data$area.sub.T1),3)," sub.std=",round(mean(data$area.sub.st),3),
                    " %trees.x=",round(mean(data$perc.trees.x.std),0),sep=""),cex.main=0.8)
  
} 

for(plot in unique(x.data$PlotCode)){
  
  data<-x.data[x.data$PlotCode==plot,]
  plot(data$x,data$y,xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y),max(data$y)),
       main=paste(unique(data$region)," ",plot," tot.area=",round(mean(data$PlotArea.use),2)," ha",
                  " sub.fld=",round(mean(data$area.sub.T1),3)," sub.std=",round(mean(data$area.sub.st),3),
                  " %trees.x=",round(mean(data$perc.trees.x.std),0),sep=""),cex.main=0.8)
  
} 

dev.off()

####################
###SAME FOR SUBPLOTS
#All plots with 0.04 ha plots, to check subplot sizes.

#Check that all trees are placed in a subplot
T1.data<-alldata[alldata$area.sub.T1==0.04,]
check<-T1.data[is.na(T1.data$SubPlotT1),] #36 trees (in MDC-02, MDC-03, MDC-04) have no subplot listed (and no std).
std.data<-alldata[alldata$area.sub.st==0.04,]
check<-std.data[is.na(std.data$Subplot_Standard),] #1 tree (in RST-01) has no std subplot, but has T1 listed.

x.std.data<-alldata[!is.na(alldata$x_standard) & !is.na(alldata$y_standard) & alldata$area.sub.st==0.04,]
x.data<-alldata[(!is.na(alldata$x) & !is.na(alldata$y) & alldata$area.sub.T1==0.04 & 
                  !(alldata$PlotCode %in% unique(x.std.data$PlotCode)))|
                (!is.na(alldata$x) & !is.na(alldata$y) & alldata$area.sub.st==0.04 & 
                     !(alldata$PlotCode %in% unique(x.std.data$PlotCode))),]

pdf("Scatter subplots xy coordinates 0.04 ha.pdf", width = 10, height = 8)

par(mfrow=c(2,2))

#Select plots and plot coordinates per subplot

for(plot in unique(x.std.data$PlotCode)){
  
  plot.data<-x.std.data[x.std.data$PlotCode==plot,] 
  
  for(subplot in unique(plot.data$Subplot_Standard)){
    
    data<-plot.data[plot.data$Subplot_Standard==subplot & !is.na(plot.data$x_standard) & 
                        !is.na(plot.data$y_standard),]
    plot(data$x_standard,data$y_standard,xlim=c(min(data$x_standard),max(data$x_standard)),
                    ylim=c(min(data$y_standard),max(data$y_standard)),
                    main=paste(unique(data$region)," ",plot," subplot.std=",subplot," no.subplots.std=",
                    mean(data$no.subplots.standard),
                    " %trees.x=",round(mean(data$perc.trees.x),0),sep=""),cex.main=0.8)
    
  }
}

for(plot in unique(x.data$PlotCode)){
  
  plot.data<-x.data[x.data$PlotCode==plot,] 
  
  for(subplot in unique(plot.data$SubPlotT1)){
    
    data<-plot.data[plot.data$SubPlotT1==subplot & !is.na(plot.data$x) & !is.na(plot.data$y),]
    plot(data$x,data$y,xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y),max(data$y)),
         main=paste(unique(data$region)," ",plot," subplot.std=",subplot," no.subplots.std=",
                    mean(data$no.subplots.T1),
                    " sub.fld=",round(mean(data$area.sub.T1),3),
                    " %trees.x=",round(mean(data$perc.trees.x),0),sep=""),cex.main=0.8)
    
  }
} 


for(plot in c("FMH-03")){
  
  plot.data<-x.data[x.data$PlotCode==plot,] 
  
  for(subplot in unique(plot.data$Subplot_Standard)){
    
    data<-plot.data[plot.data$Subplot_Standard==subplot & !is.na(plot.data$x) & !is.na(plot.data$y),]
    plot(data$x,data$y,xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y),max(data$y)),
         main=paste(unique(data$region)," ",plot," subplot.std=",subplot," no.subplots.std=",
                    mean(data$no.subplots.standard),
                    " sub.fld=",round(mean(data$area.sub.st),3),
                    " %trees.x=",round(mean(data$perc.trees.x),0),sep=""),cex.main=0.8)
    
  }
} 


dev.off()

#####################
#Same, but for a selection of plots with subplots that need to be checked.
x.std.data<-alldata[!is.na(alldata$x_standard) & !is.na(alldata$y_standard),]
x.data<-alldata[!is.na(alldata$x) & !is.na(alldata$y) & !is.na(alldata$SubPlotT1),]

####################################################################################################
##Quality xy coords subplots
trees.x.xy<-aggregate(alldata[,c("TreeID","x","y")],
                    list(alldata$PlotCode,alldata$SubPlotT1),function(x) length(na.omit(x)))
names(trees.x.xy)<-c("PlotCode","SubPlotT1","count.TreeID","count.x","count.y")
trees.x.std.xy<-aggregate(alldata[,c("TreeID","x_standard","y_standard")],
                      list(alldata$PlotCode,alldata$Subplot_Standard),function(x) length(na.omit(x)))
names(trees.x.std.xy)<-c("PlotCode","Subplot_Standard","count.TreeID","count.x.std","count.y.std")

trees.x.xy$perc.trees.sub.x<-(trees.x.xy$count.x/trees.x.xy$count.TreeID)*100
trees.x.xy2<-trees.x.xy[,c("PlotCode","SubPlotT1","perc.trees.sub.x")]

trees.x.std.xy$perc.trees.sub.x<-(trees.x.std.xy$count.x.std/trees.x.std.xy$count.TreeID)*100
trees.x.std.xy2<-trees.x.std.xy[,c("PlotCode","Subplot_Standard","perc.trees.sub.x")]

x.std.data2<-merge(x.std.data,trees.x.std.xy2,all.x=T)
x.data2<-merge(x.data,trees.x.xy2,all.x=T)
####################################################################################################

#Plots to include:
sel.plots.x<-c("AGJ-01","ALM-01","ASN-02","ASN-04","BAC-01","BAC-02",
               "BAC-03","BAC-04","BAC-05","BAC-06","CLA-03","CLA-04",
               "ELD-01","ELD-02","ELD-03","ELD-04","MNU-01","MNU-03",
               "MNU-04","MNU-06","MNU-08","OGI-07","RIO-01","RIO-02",
               "TAM-06","TIP-01","TIP-02","TIP-03","VCR-01")

sel.plots.x.std<-c("ELD-01","ELD-02")

pdf("Scatter subplots xy coordinates rest.pdf", width = 10, height = 8)

par(mfrow=c(2,2))

#Select plots and plot coordinates per subplot

for(plot in sel.plots.x){
  
  plot.data<-x.data2[x.data2$PlotCode==plot,] 
  
 for(subplot in unique(plot.data$SubPlotT1)){
  
  data<-plot.data[plot.data$SubPlotT1==subplot & !is.na(plot.data$x) & !is.na(plot.data$y),]
  plot(data$x,data$y,xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y),max(data$y)),
       main=paste(unique(data$region)," ",plot," subplot.T1=",subplot," no.subplots=",mean(data$no.subplots.T1),
                  " sub.fld=",round(mean(data$area.sub.T1),3),
                  " %trees.x=",round(mean(data$perc.trees.sub.x),0),sep=""),cex.main=0.8)

 }
}

for(plot in sel.plots.x.std){
  
  plot.data<-x.std.data2[x.std.data2$PlotCode==plot,] 
  
  for(subplot in unique(plot.data$Subplot_Standard)){
    
    data<-plot.data[plot.data$Subplot_Standard==subplot & !is.na(plot.data$x_standard) & 
                        !is.na(plot.data$y_standard),]
    plot(data$x_standard,data$y_standard,xlim=c(min(data$x_standard),max(data$x_standard)),
         ylim=c(min(data$y_standard),max(data$y_standard)),
         main=paste(unique(data$region)," ",plot," subplot.std=",subplot," no.subplots=",
                    mean(data$no.subplots.standard),
                    " sub.std=",round(mean(data$area.sub.st),3),
                    " %trees.x=",round(mean(data$perc.trees.sub.x),0),sep=""),cex.main=0.8)
    
  }
} 

dev.off()

######################################################################################################
######################################################################################################
##Look at subplot BA per subplot category
#subplot.type: standard, similar, large, small

alldata<-read.table("All data for analysis.txt",h=T)

agg.alldata<-aggregate(alldata[,c("subplotBA.m2ha","subplot.no.trees")],list(alldata$subplot.type,
                       alldata$subplot),mean,na.rm=T)
names(agg.alldata)<-c("subplot.type","subplot","subplotBA.m2ha","subplot.no.trees")

standard<-agg.alldata[agg.alldata$subplot.type=="standard",]
similar<-agg.alldata[agg.alldata$subplot.type=="similar",]
large<-agg.alldata[agg.alldata$subplot.type=="large",]
small<-agg.alldata[agg.alldata$subplot.type=="small",]

hist(standard$subplotBA.m2ha)
summary(standard$subplotBA.m2ha)  #1.6-93.5 m2/ha mean=27.5
summary(standard$subplot.no.trees) #4-83 mean=25.2

hist(similar$subplotBA.m2ha)
summary(similar$subplotBA.m2ha)  #1.6-182.50 m2/ha mean=31.1
summary(similar$subplot.no.trees) #2-81 mean=32.3

hist(large$subplotBA.m2ha)
summary(large$subplotBA.m2ha)  #9.3-75.1 m2/ha mean=33.7 
summary(large$subplot.no.trees) #15-83 mean=50

hist(small$subplotBA.m2ha)
summary(small$subplotBA.m2ha)  #0-135.7 m2/ha m mean=26.45
summary(small$subplot.no.trees) #1-54 mean=10.9

##############################################
##Histograms
gdata<-alldata[!is.na(alldata$dbhgrowth),]
gdata$PlotCode<-factor(gdata$PlotCode)

pdf("Histograms dbhgrowth.pdf", width = 10, height = 8)

par(mfrow=c(3,3),oma=c(0,0,3,0))

for(plot in unique(gdata$PlotCode)){
    
    data<-gdata[gdata$PlotCode==plot,]
    hist(data$dbhgrowth,breaks=seq(-7,7,l=50),
         main=paste(unique(data$region)," PlotCode= ",plot," area(ha)= ",round(mean(data$PlotArea.use),2),sep=""))
    
} 

dev.off()

##Overall
#Look at growth, mortality, recruitment, biomass, etc.

