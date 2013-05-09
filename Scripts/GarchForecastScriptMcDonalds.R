SPLIT_FIT_TEST_SET<-0.6;

plot(date,logRet,type='l');
savePlot("dataPlot.png");

ylim<-c(min(logRet),max(logRet));
plotForecastedData <- function(modelSTDFunction,filePrefix,numberOfPlots,color)
{  
  createWorkingSubDirecory(filePrefix); 
  lengthOfFitSet=round(length(logRet) * SPLIT_FIT_TEST_SET);
  lengthOfTestSet=length(logRet)-lengthOfFitSet;
  for(i in 1:numberOfPlots)
  {     
    
    modelData <-garchForecast(logRet=logRet[1:lengthOfFitSet],
                              sampleLength=lengthOfTestSet,
                              modelFunction=modelSTDFunction,
                              params=garchParams, 
                              numberOfStdDensityParams=length(initParams));
    plot(date,logRet,type='l');
    lines(date[(lengthOfFitSet+1):length(logRet)],modelData,type='l',col=color);
    savePlot(paste(filePrefix,"AndDataPlotForecast",i,".png",sep=""));
  }
}

#Normal Garch
initParams<-c(1,1);
garchParams <- c(2.228486,0.1155909,2.074909e-06,0.06697215,0.927272);


plotForecastedData(function(n,params) { rnorm(n,0,1); },filePrefix="Normal",
                color='red',numberOfPlots=50);

#CTS Garch
initParams=c(1.5,1,1,1,1,1);
garchParams <- c(0.486877,1.095581,0.9759551,1.715976,0.8991975,
                 1.121384,1.568464e-06,0.05627618,0.9386918)


plotForecastedData(modelSTDFunctionCTS,filePrefix="CTS",
                color='blue',numberOfPlots=50);

#MTS Garch
initParams=c(1.5,1,1,1,1);
garchParams <- c(0.5977937,1.360131,1.108706,0.9535164,0.9114298
                 5.841102e-06,0.1199947,0.8602524);


plotForecastedData(modelSTDFunctionMTS,filePrefix="MTS",
                color='green',numberOfPlots=50);

#KR Garch
initParams=c(1.5,1,1,1,1,1,1,1);
garchParams <- c(0.8774809,0.9139359,0.7632551,1.191551,0.8819413,
                 0.8758534,1.372726,1.759472,1.884582e-06,0.07497817,
                 0.9247544)

plotForecastedData(modelSTDFunctionKR,filePrefix="KR",
                color='purple',numberOfPlots=50);
