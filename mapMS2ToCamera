require(intervals)
ppmCal<-function(run,ppm)
{
  return((run*ppm)/1000000)
}
IntervalMerge<-function(cameraObject,MSMSdata, PlusTime,MinusTime,ppm,listOfMS2Mapped=list(),listOfUnMapped=list()){
  
  listofPrecursorsmz<-c()
  for(i in seq(1,length(MSMSdata)))
  {
    listofPrecursorsmz<-c(listofPrecursorsmz,MSMSdata[[i]]@precursorMz)
  }
  
  listofPrecursorsrt<-c()
  for(i in seq(1,length(MSMSdata)))
  {
    listofPrecursorsrt<-c(listofPrecursorsrt,MSMSdata[[i]]@rt)
  }
  
  CameramzColumnIndex<-which(colnames(cameraObject@groupInfo)=="mz")
  
  MassRun1<-Intervals_full(cbind(listofPrecursorsmz,listofPrecursorsmz))
  
  MassRun2<-Intervals_full(cbind(cameraObject@groupInfo[,CameramzColumnIndex]-
                                   ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm),
                                 cameraObject@groupInfo[,CameramzColumnIndex]+
                                   ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm)))
  
  Mass_iii <- interval_overlap(MassRun1,MassRun2)
  
  CamerartLowColumnIndex<-which(colnames(cameraObject@groupInfo)=="rtmin")
  CamerartHighColumnIndex<-which(colnames(cameraObject@groupInfo)=="rtmax")
  
  TimeRun1<-Intervals_full(cbind(listofPrecursorsrt,listofPrecursorsrt))
  
  TimeRun2<-Intervals_full(cbind(cameraObject@groupInfo[,CamerartLowColumnIndex]-MinusTime,
                                 cameraObject@groupInfo[,CamerartHighColumnIndex]+PlusTime))
  Time_ii <- interval_overlap(TimeRun1,TimeRun2)
  
  imatch = mapply(intersect,Time_ii,Mass_iii)
  
  for (i in 1:length(imatch)) {
    for(j in imatch[[i]])
    {
      
      listOfMS2Mapped[[as.character(j)]]<-
        c(listOfMS2Mapped[[as.character(j)]],MSMSdata[[i]])
    }
  }
  for (i in 1:length(imatch)) {
    
    if(length(imatch[[i]])==0)
    {
      listOfUnMapped<-c(listOfUnMapped,MSMSdata[[i]])
    }
  }
  return(list(mapped=listOfMS2Mapped,unmapped=listOfUnMapped))
}
