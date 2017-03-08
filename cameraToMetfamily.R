cameraToMetfamily<-function(Intensities=NA,dataGroups=NA,sampleType=NA,injectionOrder=NA,mzAndRT,
                            imputeNAs=T)
{
  headersTMP<-c("m/z",
                "RT",
                "Annotation",
                "Alignment ID",
                "Average Rt(min)",
                "Average Mz",
                "Metabolite name",
                "Adduct ion name",
                "Fill %",
                "MS/MS included",
                "INCHIKEY",
                "SMILES",
                "LINK",
                "Dot product",
                "Reverse dot product",
                "Fragment presence %",
                "Spectrum reference file name")
  
  if(imputeNAs)
  {
    require(mixOmics)
    imputed.X = (mixOmics::nipals(t(Intensities), reconst = TRUE)$rec)
    
    #  only replace the imputation for the missing values
    nas = is.na(t(Intensities))
    imputed.X[!nas] = t(Intensities)[!nas]
    
    Intensities<-t(imputed.X)
  }
  dataTmp<-matrix("0",nrow = nrow(Intensities)+4,ncol=ncol(dataAbudance)+length(headersTMP))
  
  dataTmp[1,]<-c(rep("",length(headersTMP)-1),"Class",dataGroups)
  dataTmp[2,]<-c(rep("",length(headersTMP)-1),"Type",sampleType)
  dataTmp[3,]<-c(rep("",length(headersTMP)-1),"Injection order",injectionOrder)
  dataTmp[4,]<-c(headersTMP,names(dataAbudance))
  dataTmp[5:nrow(dataTmp),1]<-mzAndRT[,"mz"]# index from the header
  dataTmp[5:nrow(dataTmp),2]<-mzAndRT[,"rt"]# index from the header
  dataTmp[5:nrow(dataTmp),5]<-mzAndRT[,"rt"]# index from the header
  dataTmp[5:nrow(dataTmp),6]<-mzAndRT[,"mz"]# index from the header
  dataTmp[5:nrow(dataTmp),seq(length(headersTMP)+1,ncol((dataTmp)))]<-as.matrix(Intensities)
  
  return(dataTmp)
}
