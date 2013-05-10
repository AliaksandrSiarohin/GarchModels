defaultDirectory<-"D:\\Data\\Garch Source\\SNP500";
setwd (defaultDirectory);
logRet<-c();
date<-c();

data<-loadPricesFromFile("prices.csv");
date<-as.Date(data[1][[1]],"%d-%m-%Y");
date<-date[2:length(date)];
logRet<-transformPrices(data[2][[1]]);




createWorkingSubDirecory<-function(dirName)
{
  dir.create(file.path(defaultDirectory, dirName), showWarnings = FALSE)
  setwd(file.path(defaultDirectory, dirName));
}

#load prices from .csv file
loadPricesFromFile<-function(fileName=
  "S.csv")
{
  result<-read.csv(paste(defaultDirectory,fileName,sep="\\"));
  result$Date=as.Date(result$Date,format="%d.%m.%Y");
  return (result);
}

#Save result of garch model optimosing to file
saveParamsToFile<-function(params,fileName="optimRes.txt")
{
  write(params,fileName,sep=',');
}

#Save modeled garch data to file
saveDataToFile<-function(data,fileName="optimRes.txt")
{
  write.csv(data,fileName);
}

#save plot to file
savePlot<-function(fileName)
{
  dev.copy(function(...){return (png(width=1920,height=1080,...))}
           ,fileName);
  dev.off ();
}