###############################################################################
###############################################################################
#Edit Chave database
WD.Chave<-read.table("C:/Users/rozendad/Dropbox/Rainfor data edited/Wood density data Chave.txt",h=T)

#Add species average per region
avg.spp<-aggregate(WD.Chave$WD,list(WD.Chave$region,WD.Chave$family,
                                    WD.Chave$binomial,WD.Chave$genus,WD.Chave$species),
                   mean,na.rm=T)
names(avg.spp)<-c("region","family","binomial","genus","species","WDspp.chave")

#Add overall species average
avg.spp2<-aggregate(WD.Chave$WD,list(WD.Chave$family,
                                     WD.Chave$binomial,WD.Chave$genus,WD.Chave$species),
                    mean,na.rm=T)
names(avg.spp2)<-c("family","binomial","genus","species","WDspp.chave2")

#Add genus average per region 
avg.gen<-aggregate(WD.Chave$WD,list(WD.Chave$region,WD.Chave$family,WD.Chave$genus),
                   mean,na.rm=T)
names(avg.gen)<-c("region","family","genus","WDgen.chave")

#Add overall genus average 
avg.gen2<-aggregate(WD.Chave$WD,list(WD.Chave$family,WD.Chave$genus),
                    mean,na.rm=T)
names(avg.gen2)<-c("family","genus","WDgen.chave2")

#Add family average per region 
avg.fam<-aggregate(WD.Chave$WD,list(WD.Chave$region,WD.Chave$family),
                   mean,na.rm=T)
names(avg.fam)<-c("region","family","WDfam.chave")

#Add overall family average 
avg.fam2<-aggregate(WD.Chave$WD,list(WD.Chave$family),
                    mean,na.rm=T)
names(avg.fam2)<-c("family","WDfam.chave2")

#Merge
WD.Chave2<-merge(avg.spp,avg.spp2,all.x=T)
for (df in list(avg.gen,avg.gen2,avg.fam,avg.fam2)){
  WD.Chave2<-merge(WD.Chave2,df,all.x=T)
} 

write.table(WD.Chave2,"WD_Chave.txt",row.names=F,quote=F,sep="\t")