writeMSMStoMSP<-function(MSMS=NA,fileName)
{
  #MSMS<-ss
  outPut<-""
  for(MS2n in c(1:length(MSMS)))
  {
    MS2<-unlist(MSMS[[MS2n]])
    pattern<-"NAME: @NAME\nRETENTIONTIME: @RETENTIONTIME\nPRECURSORMZ: @PRECURSORMZ\nMETABOLITENAME: @METABOLITENAME\nADDUCTIONNAME: @ADDUCTIONNAME\nNum Peaks: @Num Peaks\n@MSs"
    pattern<-gsub("@NAME","Unknown",pattern)
    pattern<-gsub(pattern = "@RETENTIONTIME",replacement = MS2@rt,x = pattern,fixed = T)
    pattern<-gsub(pattern = "@PRECURSORMZ",replacement = MS2@precursorMz,x = pattern,fixed = T)
    pattern<-gsub(pattern = "@METABOLITENAME",replacement = "",x = pattern,fixed = T)
    pattern<-gsub(pattern = "@ADDUCTIONNAME",replacement = "",x = pattern,fixed = T)   
    pattern<-gsub(pattern = "@Num Peaks",replacement = MS2@peaksCount,x = pattern,fixed = T)
    pattern<-gsub(pattern = "@MSs",replacement =    paste(apply( cbind(MS2@mz,MS2@intensity) , 1 , paste , collapse = "\t" ),collapse = "\n")
,x = pattern,fixed = T)
   
    outPut<-paste(outPut,pattern,sep = "\n\n")
    }
  cat(outPut,file = fileName)
}
