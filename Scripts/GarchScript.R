
plot(date,logRet,type='l');
savePlot("dataPlot.png")

ylim<-c(min(logRet),max(logRet));
plotModeledData <- function(modelSTDFunction,filePrefix,numberOfPlots,color)
{
  
  createWorkingSubDirecory(filePrefix);
  saveParamsToFile(garchParams,paste(filePrefix,"GP.txt"));
  for(i in 1:numberOfPlots)
  {     
    
    modelData <- garchModel(length(logRet),modelSTDFunction,
                            garchParams, length(initParams));
   
    saveDataToFile(modelData,paste(filePrefix,"ModelData",i,".txt",sep=""));
    
    plot(modelData,type='l',col=color,ylim=ylim);
    savePlot(paste(filePrefix,"Plot",i,".png",sep=""));
    
    plot(date,logRet,type='l');
    lines(date,modelData,type='l',col=color);
    savePlot(paste(filePrefix,"AndDataPlot",i,".png",sep=""));
  }
}

#Normal Garch
initParams<-c(1,1);
garchParams <- fitGarchModel(logRet,stdNormal,initParams,control=list(maxit=5000))$par;


plotModeledData(function(n,params) { rnorm(n,0,1); },filePrefix="Normal",
               color='red',numberOfPlots=50);

#CTS Garch
initParams=c(1.5,1,1,1,1,1);
garchParams <- fitGarchModel(logRet,stdCTS,initParams,control=list(maxit=5000))$par;


plotModeledData(modelSTDFunctionCTS,filePrefix="CTS",
                color='blue',numberOfPlots=50);

#MTS Garch
initParams=c(1.5,1,1,1,1);
garchParams <- fitGarchModel(logRet,stdMTS,initParams,control=list(maxit=15000))$par;


plotModeledData(modelSTDFunctionMTS,filePrefix="MTS",
                color='green',numberOfPlots=50);

#KR Garch
initParams=c(1.5,1,1,1,1,1,1,1);
garchParams <- fitGarchModel(logRet,stdKR,initParams,control=list(maxit=5000))$par;

plotModeledData(modelSTDFunctionKR,filePrefix="KR",
                color='purple',numberOfPlots=50);
