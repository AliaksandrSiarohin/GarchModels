require("elliptic");
require("hypergeo");
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
garchParams <- c(1.211927,1.974351,6.164198e-06,0.08484688,0.9024538);


plotForecastedData(function(n,params) { rnorm(n,0,1); },filePrefix="Normal",
                color='red',numberOfPlots=50);

#CTS Garch
initParams=c(1.5,1,1,1,1,1);
garchParams <- c(1.284153,1.009822,1.00666,1.243156,1.139199,
                 1.248145,6.251535e-06,0.1621888,0.8286405)


plotForecastedData(modelSTDFunctionCTS,filePrefix="CTS",
                color='blue',numberOfPlots=50);

#MTS Garch
initParams=c(1.5,1,1,1,1);
garchParams <- c(1.147742,-0.8670031,1.07572,0.5348109,0.473992,
                 2.832777e-06,0.1272681,0.8597438)


plotForecastedData(modelSTDFunctionMTS,filePrefix="MTS",
                color='green',numberOfPlots=50);

#KR Garch
initParams=c(1.5,1,1,1,1,1,1,1);
garchParams <- c(1.461509,1.092585,0.8305387,0.9899253,0.9606343,
                 1.071105,1.08335,5.667672e-05,5.667672e-06,0.1103979,
                 0.8663954)

plotForecastedData(modelSTDFunctionKR,filePrefix="KR",
                color='purple',numberOfPlots=50);
