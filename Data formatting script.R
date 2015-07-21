#FORMAT RAINFOR DATA

library(reshape)

##########################################################################
##FORMAT TREE DATA########################################################
#Load raw data

dropbox <- "C:/Users/DMAR/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/rozendad/Dropbox/Rainfor data edited"
#dropbox <- "C:/Users/vande202/Dropbox/Rainfor data edited"

org.data<-read.table(paste(dropbox,"Treedata prep.txt",sep="/"),h=T)

#PlotCode MAK-10 not in plotdata (has MAK-01 instead...)? Exclude for now.
data<-org.data[org.data$PlotCode!="MAK-10",]

###Look at characteristics raw data

#Census numbers
census<-aggregate(data$DBH4,list(data$PlotCode,data$censusno),mean,na.rm=T)
names(census)<-c("PlotCode","censusno","avgDBH4")
census2<-census[order(census$PlotCode,census$censusno),]

min.census<-aggregate(census2$censusno,list(census2$PlotCode),min,na.rm=T)
names(min.census)<-c("PlotCode","min.censusno")

census.no<-merge(census2,min.census)
census.no$censusno.new<-ifelse(census.no$censusno==census.no$min.censusno,1,2)
census.no2<-census.no[,c("PlotCode","censusno","censusno.new")]

#Merge censusno.new into data, select columns, and reshape
data2<-merge(data,census.no2,all.x=T)
data3<-data2[,c("PlotCode","censusno","TreeID","TagNumber","FamilyAPGName","GenusName",
                "FullSpeciesName","plotviewid","PlotID","PlotViewPlotCensusID","Alive",
                "MeanDecimalDate","DBH4","SubPlotT1","Subplot_Standard","x","y",
                "x_standard","y_standard","truerecruit","censusno.new")]

#TreeID is unique across dataset, TagNumber not
length(unique(data3$TreeID))
length(unique(data3$TagNumber))

#Length is the same
length(unique(data3$PlotCode))
length(unique(data3$plotviewid))

#Add census count to data
census.count<-aggregate(data3$censusno.new, list(data3$TreeID), function(x) length(na.omit(x)))
names(census.count)<-c("TreeID","census.count")

data4<-merge(data3,census.count,all.x=T)
dim(data3)
dim(data4)

##Trees with two censuses
data4.2<-data4[data4$census.count==2,]

#Does not work correctly
#data4.2r<-reshape(data4.2, idvar="TreeID", varying=c("censusno","Alive","MeanDecimalDate","DBH4","truerecruit"),
#                  timevar="censusno.new",direction="wide",sep=".")
   
#Very ugly... format reshaped dataframe               
data4.2r<-reshape(data4.2, idvar="TreeID",timevar="censusno.new",direction="wide",sep=".")

data4.2fin<-data4.2r[,c("TreeID","PlotCode.1","TagNumber.1","FamilyAPGName.1",
                        "GenusName.1","FullSpeciesName.1","plotviewid.1","PlotID.1",
                        "PlotViewPlotCensusID.1","SubPlotT1.1","Subplot_Standard.1","x.1",
                        "y.1","x_standard.1","y_standard.1","truerecruit.1","census.count.1",
                        "Alive.1","Alive.2","DBH4.1","DBH4.2")]

names(data4.2fin)<-c("TreeID","PlotCode","TagNumber","FamilyAPGName",
                     "GenusName","FullSpeciesName","plotviewid","PlotID",
                     "PlotViewPlotCensusID","SubPlotT1","Subplot_Standard","x",
                     "y","x_standard","y_standard","truerecruit","census.count",
                     "Alive.c1","Alive.c2","DBH4.c1","DBH4.c2")

##ADD demography category
#1) Alive in census 1 and 2
#2) Alive in census 1, dead in census 2
#3) Alive in census 1, missing in census 2. Died before census 2.
#4) Missing in census 1, alive in census 2. Optionally a recruit.
#5) Missing in census 1, alive in census 2, and truerecruit = 1. Recruit.
min(data4.2fin$Alive.c1)  #All trees are alive in first census
data4.2fin$Status<-ifelse(data4.2fin$Alive.c2==1,1,2)

##Trees with one census
#Trees present in first census only
data4.1_1<-data4[data4$census.count==1 & data4$censusno.new==1,]
data4.1_1$Alive.c1<-data4.1_1$Alive
data4.1_1$Alive.c2<-0
data4.1_1$Status<-3
data4.1_1$DBH4.c1<-data4.1_1$DBH4
data4.1_1$DBH4.c2<-NA

#Trees present in second census only
data4.1_2<-data4[data4$census.count==1 & data4$censusno.new==2,]
data4.1_2$Alive.c2<-data4.1_2$Alive
data4.1_2$Alive.c1<-0
data4.1_2$Status<-ifelse(data4.1_2$truerecruit==0,4,5)
data4.1_2$DBH4.c1<-NA
data4.1_2$DBH4.c2<-data4.1_2$DBH4

data4.1tot<-rbind(data4.1_1,data4.1_2)
data4.1fin<-data4.1tot[,c("TreeID","PlotCode","TagNumber","FamilyAPGName",
                "GenusName","FullSpeciesName","plotviewid","PlotID",
                "PlotViewPlotCensusID","SubPlotT1","Subplot_Standard","x",
                "y","x_standard","y_standard","truerecruit","census.count",
                "Alive.c1","Alive.c2","DBH4.c1","DBH4.c2","Status")]

#Combine the two datasets again
treedata<-rbind(data4.1fin,data4.2fin)

##Look at subplots and presence of xy-coordinates
#Compare SubPlotT1 and Subplot_Standard
no.subplots<-aggregate(treedata[,c("SubPlotT1","Subplot_Standard","x","y",
                                   "x_standard","y_standard")],list(treedata$PlotCode),
                                   function(x) length(unique(na.omit(x))))
names(no.subplots)<-c("PlotCode","no.subplots.T1","no.subplots.standard","count.x",
                      "count.y","count.x.std","count.y.std")

##################################################################################
##ADD PLOT DATA AND FORMAT########################################################
#Load raw data

pdata<-read.table(paste(dropbox,"Plotdata prep.txt",sep="/"),h=T)
pdata$region<-pdata$ContinentName

#Add E and CWD for all plots for biomass estimates (Chave et al. 2014 GCB) 
#Deleted.

#Add in number of subplots and subplot area
pdata2<-merge(pdata,no.subplots)
pdata2$area.sub.T1<-pdata2$PlotArea.use/pdata2$no.subplots.T1
pdata2$area.sub.st<-pdata2$PlotArea.use/pdata2$no.subplots.standard

#Add calendar year
pdata2$CensusDate<-(pdata2$MeanDecimalDate.c1+pdata2$MeanDecimalDate.c2)/2

##################################################################################
##COMBINE TREE DATA AND PLOT DATA#################################################
#Discard PlotID column, they do not match between treedata and plotdata in some cases
pdata2$PlotID<-NULL
treedata2<-treedata
treedata2$PlotID<-NULL

totdata<-merge(treedata2,pdata2,all.x=T)

##################################################################################
##FORMAT TOTAL DATA###############################################################
#Number of plots with xy coords
xydata<-totdata[!(totdata$count.x==0 & totdata$count.y==0 & totdata$count.x.std==0 & 
                    totdata$count.y.std==0),]
length(unique(xydata$PlotCode))     #151 plots have xy coords

#Quality xy coords: missing data!
trees.xy<-aggregate(xydata[,c("TreeID","x","y","x_standard","y_standard")],
                      list(xydata$PlotCode),function(x) length(na.omit(x)))
names(trees.xy)<-c("PlotCode","count.TreeID","count.x",
                      "count.y","count.x.std","count.y.std")
trees.xy$perc.trees.x<-(trees.xy$count.x/trees.xy$count.TreeID)*100
trees.xy$perc.trees.x.std<-(trees.xy$count.x.std/trees.xy$count.TreeID)*100
trees.xy<-trees.xy[,c("PlotCode","count.TreeID","perc.trees.x","perc.trees.x.std")]
summary(trees.xy)

###GO BACK to plot data, add trees.xy
##Count number of suitable plots based on subplot size etc. and xy coords
pdata3<-merge(pdata2,trees.xy,by="PlotCode",all.x=T)

#Look at subplots / xy coords (at least 90% of trees mapped)
pdata3$xy.90<-ifelse(!is.na(pdata3$perc.trees.x) & !is.na(pdata3$perc.trees.x.std) & 
                  (pdata3$perc.trees.x>=90 | pdata3$perc.trees.x.std>=90),1,0)
xy.90<-pdata3[pdata3$xy.90==1,]       #105 plots
rest<-pdata3[pdata3$xy.90==0,]        #103 plots
length(rest[rest$area.sub.T1==0.04 | rest$area.sub.st==0.04,]) #31 plots

#Save plot metadata (pdata3)
write.table(pdata3,"Plotdata.txt",row.names=F,quote=F,sep="\t")

#Add in pdata3
totdata2<-merge(totdata,pdata3[,c("PlotCode","perc.trees.x","perc.trees.x.std")],by="PlotCode",all.x=T)

##Add in basal area and other variables
totdata2$dbh0<-totdata2$DBH4.c1/10
totdata2$dbh1<-totdata2$DBH4.c2/10
totdata2$BA0<-pi*((totdata2$dbh0/2)^2)
totdata2$BA1<-pi*((totdata2$dbh1/2)^2)

##Add in demographics
#dbhgrowth
totdata2$dbhgrowth<-ifelse(totdata2$Status==1,(totdata2$dbh1-totdata2$dbh0)/totdata2$IntervalLength,NA)
#1=dead. Correct for IntervalLength later.
totdata2$dead<-ifelse(totdata2$Status==2 | totdata2$Status==3,1,0) 
#1=recruit. Correct for IntervalLength and subplot area later.
totdata2$recruit<-ifelse(totdata2$Status==4 | totdata2$Status==5,1,0) 

##Change column names
totdata2$family<-totdata2$FamilyAPGName
totdata2$genus<-totdata2$GenusName
totdata2$binomial<-totdata2$FullSpeciesName

#######################
##Add in WD for spp per region, otherwise a genus or family average per region
WD.chave<-read.table("WD_Chave.txt",h=T)

#Check spelling of species, genus, and family names that are not present in 
#the WD databases
#Chave database species
if(any(!(totdata2$binomial %in% (WD.chave$binomial)))){
  print(unique(totdata2[!(totdata2$binomial %in% WD.chave$binomial),]$binomial)) 
}

#Chave database genus
if(any(!(totdata2$genus %in% (WD.chave$genus)))){
  print(unique(totdata2[!(totdata2$genus %in% WD.chave$genus),]$genus)) 
}

#Chave database family
if(any(!(totdata2$family %in% (WD.chave$family)))){
  print(unique(totdata2[!(totdata2$family %in% WD.chave$family),]$family)) 
}

##Select WD per continent
regions<-c("Africa","South.America")

vector1<-c()
for (reg in regions) {

    sel.data<-totdata2[totdata2$region==reg,]
    WD.chave2<-WD.chave[WD.chave$region==reg,]

    #Add columns with different wood density values
    #Make sure comparison column is a character
    sel.data$binomial<-as.character(sel.data$binomial)
    WD.chave2$binomial<-as.character(WD.chave2$binomial)
    sel.data$genus<-as.character(sel.data$genus)
    WD.chave2$genus<-as.character(WD.chave2$genus)
    sel.data$family<-as.character(sel.data$family)
    WD.chave2$family<-as.character(WD.chave2$family)

    for(j in 1:nrow(sel.data)){
      sel.data[j,"WDspp.chave"]<-mean(WD.chave2[WD.chave2$binomial==
                                      sel.data[j,"binomial"],"WDspp.chave"])
    }

    for(j in 1:nrow(sel.data)){
     sel.data[j,"WDgen.chave"]<-mean(WD.chave2[WD.chave2$genus==
                                      sel.data[j,"genus"],"WDgen.chave"])
    }
    
    for(j in 1:nrow(sel.data)){
      sel.data[j,"WDfam.chave"]<-mean(WD.chave2[WD.chave2$family==
                                                  sel.data[j,"family"],"WDfam.chave"])
    }

    #Select correct wood density value
    sel.data$WD2<-ifelse(!is.na(sel.data$WDspp.chave),sel.data$WDspp.chave,
                         ifelse(!is.na(sel.data$WDgen.chave),sel.data$WDgen.chave,
                         ifelse(!is.na(sel.data$WDfam.chave),sel.data$WDfam.chave,NA)))

    ###################################################
    #Add average WD of all (identified) stems per plot (following Brienen et al. 2015) 
                     
    plotWD<-aggregate(sel.data$WD2,list(sel.data$PlotCode),mean,na.rm=T)
    names(plotWD)<-c("PlotCode","meanWD")
    
    sel.data2<-merge(sel.data,plotWD,all.x=T)
                     
    vector1<-rbind(vector1,sel.data2)
}

totdata3<-vector1    
    
#Select correct WD value                     
totdata3$WD3<-ifelse(is.na(totdata3$WD2),totdata3$meanWD,totdata3$WD2)
totdata3$WD<-totdata3$WD3

#Replace with equation(s) from Feldpausch et al.
##Add biomass: new Chave equation (Chave et al. 2014 GCB)
#totdata3$biomass0<-exp(-1.803 - (0.976 * totdata3$E) + 
#                                  0.976 * log(totdata3$WD) + 2.673 * log(totdata3$dbh0) - 
#                                  0.0299 * ((log(totdata3$dbh0))^2))

##################################################################################
##FORMAT SUBPLOTS FOR CALCULATING "NEIGHBOURHOOD" INDICES

st.subplots<-c("ELD-01","ELD-02","TAM-06","TIP-01")

  
field.subplots<-c("AGJ-01","ALM-01","ASN-02","ASN-04","BES-01","BOG-01","BOG-02","CAX-01",
                  "CAX-02","CLA-03","CLA-04","HAB-02","HAB-03","HAB-04","HAB-05","HAB-06",
                  "HAB-07","LAS-02","LKM-05","MDC-02","MDC-03","MDC-04","MNU-01","MNU-03",
                  "MNU-04","MNU-05","MNU-06","MNU-08","MRB-02","MRB-03","NXV-06","NXV-07",
                  "NXV-08","OGI-04","OGI-05","OGI-06","OGI-07","PPB-01","PPB-02","PPB-03",
                  "PTB-01","PTB-02","RIO-01","RIO-02","TAM-03","VCR-01","VCR-02","YAN-01")

#New columns subplot, subplot area, and subplot ID
totdata3$subplot<-ifelse((totdata3$area.sub.st==0.04 | totdata3$PlotCode %in% st.subplots),
                          totdata3$Subplot_Standard,
                         ifelse(totdata3$PlotCode %in% field.subplots,totdata3$SubPlotT1,NA))

totdata3$subplot.area<-ifelse((totdata3$area.sub.st==0.04 | totdata3$PlotCode %in% st.subplots),
                              0.04,ifelse(totdata3$PlotCode %in% field.subplots,
                              totdata3$area.sub.T1,NA))

totdata3$subplotID<-paste(totdata3$PlotCode,totdata3$subplot,sep="-")

#Plots to exclude: see plot metadata (NOW NOT NEEDED)
exclude<-c("ACL-01","BAC-01","BAC-02","BAC-03","BAC-04","BAC-05","BAC-06","CAI-03",
           "CAI-04","CAI-05","CAI-06","ELD-03","ELD-04","MRB-01","MSH-01",
           "POR-01","POR-02")

#Subplots close to 0.04 ha (and/or to 20x20 m)
similar.plots<-c("AGJ-01","ALM-01","ASN-02","ASN-04","BEE-01","CAX-01","CAX-02",
                 "CLA-03","CLA-04","ELD-01","ELD-02","JEN-12","MNU-03","MNU-04",
                 "MNU-05","MNU-06","MNU-08","MRB-02","MRB-03","PPB-01","PPB-02",
                 "PPB-03","PTB-01","PTB-02","RIO-01","RIO-02","TAM-03","TAM-06")

#0.09 ha or larger subplots
larger.plots<-c("HAB-02","HAB-03","HAB-04","HAB-05","HAB-06","HAB-07","LKM-05",
                "MNU-01","OGI-04","OGI-05","OGI-06","OGI-07","TIP-02","TIP-03",
                "YAN-01")

#0.01 ha subplots
small.plots<-c("BES-01","BOG-01","BOG-02","NXV-06","NXV-07","NXV-08",
               "VCR-01","VCR-02")

#Add variable to select plots based on subplot type
totdata3$subplot.type<-ifelse(totdata3$PlotCode %in% exclude,NA,
                              ifelse(totdata3$PlotCode %in% similar.plots,"similar",
                              ifelse(totdata3$PlotCode %in% larger.plots,"large",       
                              ifelse(totdata3$PlotCode %in% small.plots,"small","standard"))))

totdata4<-totdata3[!is.na(totdata3$subplot.type),]
  
#Add no. of trees per subplot and subplot basal area
agg.subplot1<-aggregate(totdata4$TreeID,list(totdata4$PlotCode,totdata4$subplot),
                        function(x) length(na.omit(x)))
names(agg.subplot1)<-c("PlotCode","subplot","subplot.no.trees")

agg.subplot2<-aggregate(totdata4$BA0,list(totdata4$PlotCode,totdata4$subplot),sum,na.rm=T)
names(agg.subplot2)<-c("PlotCode","subplot","subplotBA")

#Add subplot data, correct units
m.totdata4<-merge(totdata4,agg.subplot1,all.x=T)
totdata5<-merge(m.totdata4,agg.subplot2,all.x=T)

totdata5$subplotBA.m2ha<-(totdata5$subplotBA/10000)/totdata5$subplot.area

#Select relevant columns
totdata6<-totdata5[,c("PlotCode","subplot","TreeID","SubPlotT1","Subplot_Standard",
                      "x","y","x_standard","y_standard","truerecruit","region",
                      "PlotArea.use","IntervalLength","no.subplots.T1","no.subplots.standard",
                      "count.x","count.y","count.x.std","count.y.std","area.sub.T1","area.sub.st",
                      "CensusDate","perc.trees.x","perc.trees.x.std","dbh0","dbh1","BA0","BA1",
                      "dbhgrowth","dead","recruit","family","genus","binomial","WD","subplot.area",
                      "subplotID","subplot.type","subplot.no.trees","subplotBA.m2ha")]

#Save data file
write.table(totdata6,"All data for analysis.txt",row.names=F,quote=F,sep="\t")

#######
#######
##TO DO:
#Add rainfall + other variables (Worldclim) to plot data

######
######
##To think about:
#Mean calendar dates: census intervals range from 1968 to 2010 (mean 1997, median 1998).
#Subplots and neighbourhood structures.
#What resolution? Species-level not possible. Trait values...?

#Amazon
length(unique(data[data$region=="South.America",]$PlotCode))  #121 plots
length(unique(data[data$region=="South.America",]$binomial))  #3323 spp!!
length(unique(data[data$region=="South.America",]$TreeID))    #92624 trees

#Africa
length(unique(data[data$region=="Africa",]$PlotCode))  #62 plots
length(unique(data[data$region=="Africa",]$binomial))  #734 spp
length(unique(data[data$region=="Africa",]$TreeID))    #26668 trees

