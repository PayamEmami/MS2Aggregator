Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
groupMS2<-function(mappedMS2s=NA,unmappedMS2s=NA,cameraObject=NA,includedUnmapped=T,processUnmapped=T,mergingMethod="max",pickPeak=F,smooth=F,...)
{
  output<-list()
  for(x in names(mappedMS2s)) 
  {
    
    if(length(mappedMS2s[[x]])<2)
    {
      output[[x]]<-mappedMS2s[[x]][[1]]
      output[[x]]@centroided<-F
      
    }else
    {
    
      
      
      tempSpectrum<-new("Spectrum2")
      tempSpectrum@precursorMz<-as.numeric(cameraObject@groupInfo[as.numeric(x),"mz"])
      tempSpectrum@rt<-as.numeric(cameraObject@groupInfo[as.numeric(x),"rt"])
      tempSpectrum@precursorCharge<-Mode(sapply(mappedMS2s[[x]],function(x){x@precursorCharge}))
      tempSpectrum@collisionEnergy<-Mode(sapply(mappedMS2s[[x]],function(x){x@collisionEnergy}))
      tempSpectrum@precursorIntensity<-mean(sapply(mappedMS2s[[x]],function(x){x@precursorIntensity}))
      #tempSpectrum@msLevel<-2
     # tempSpectrum@acquisitionNum<-sapply(mappedMS2s[[x]],function(x){x@acquisitionNum})
     # tempSpectrum@scanIndex<-mode(sapply(mappedMS2s[[x]],function(x){x@scanIndex}))
      tempSpectrum@centroided<-F
      
      combinedMS<-c()
      i<-1
      for(j in mappedMS2s[[x]])
      {
        combinedMS<- rbind(combinedMS,data.frame(mz=j@mz,int=j@intensity,sample=i))
        i<-i+1
      }
      if(mergingMethod=="max")
      {
        groupIndex<-xcms:::mzClustGeneric(combinedMS[,c("mz","sample")],verbose=F)
        
        newMS2<-data.frame(matrix(nrow = length(groupIndex$idx),ncol = 2))
        names(newMS2)<-c("mz","intensity")
        for(k in 1:length(groupIndex$idx))
        {
          y<-groupIndex$idx[[k]]
          if(length(y)<2)
          {
            newMS2[k,]<-combinedMS[y,c("mz","int")]
          }else
          {
            #k<-41
            # k<-38
            selectedSpec<-combinedMS[groupIndex$idx[[k]],c("mz","int")]
            selectedSpec<-selectedSpec[order(selectedSpec$mz),]
            
            newMS2[k,]<-selectedSpec[which.max(selectedSpec$int),]
            # x<-selectedSpec[,c("mz","int")][,1]
            #y<-selectedSpec[,c("mz","int")][,2]
            
            
            
          }
        }
        newMS2<-newMS2[order(newMS2$mz),]
        tempSpectrum@mz<-newMS2$mz
        tempSpectrum@intensity<-newMS2$intensity
      }else
      {
        combinedMS<-combinedMS[order(combinedMS$mz),]
        tempSpectrum@mz<-combinedMS$mz
        tempSpectrum@intensity<-combinedMS$int
      }
      
      output[[x]]<-tempSpectrum
    }
  }
  if(pickPeak)
  {
    for(x in names(output))
    {
      output[[x]]<- MSnbase::pickPeaks(output[[x]])
      
    }
  }
  if(smooth)
  {
    for(x in names(output))
    {
      output[[x]]<- MSnbase::smooth(output[[x]])
      
    }
  }
  unmap<-list()
  if(includedUnmapped==T)
  {
    for(i in c(1:length(unmappedMS2s)))
    {
      tmp<-unmappedMS2s[[i]]
      tmp@centroided<-F
      if(pickPeak & processUnmapped)
      {
        tmp<-MSnbase::pickPeaks(tmp)
      }
      
      if(smooth & processUnmapped)
      {
        tmp<-MSnbase::smooth(tmp)
      }
      unmap<-c(unmap,tmp)
      #output<-c(output,tmp)
    }
   
  }
  
  return(list(mapped=output,unmapped=unmap))
}

