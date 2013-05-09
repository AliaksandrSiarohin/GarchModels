


#Моделирование случайной величины с заданной численно плотностью распределения
model <- function(density,n,mean,sd,numberOfCuts=1000,deviationFromMeanBound=8,lower=NULL,upper=NULL)
{  
  if(!is.null(upper))
  {
    upperBound=upper;
  }
  else
  {
    upperBound=deviationFromMeanBound*sd+mean;
  }
  if(!is.null(lower))
  {
    lowerBound=lower;
  }
  else
  {
    lowerBound=-deviationFromMeanBound*sd+mean;
  }
  step=(upperBound-lowerBound)/numberOfCuts;
  
  #Вычисляем верояность попадания на отрезок размера step от lowerBound до upperBound 
  cutBeginings<-seq(from=lowerBound,to=upperBound-step,by=step);
  
  probabilytis<-c();
  
  for(currentCut in 1:length(cutBeginings))
  {
    prob=density(cutBeginings[currentCut]+step/2);
    if (prob <= 0)
      prob=0;     
    probabilytis=c(probabilytis,prob);   
  }
  
  probabilytis=probabilytis/sum(probabilytis);
  
  randomCutsBeginings=
    sample(x=cutBeginings,size=n,replace=TRUE,prob=probabilytis);
  
  #Апроксимируем случайную величину равномерной на каждом отрезке
  
  randomResult=c();
  for(currentCut in 1:length(randomCutsBeginings))
  {
    randomVariable=runif(n=1,
                         min=randomCutsBeginings[currentCut],
                         max=randomCutsBeginings[currentCut]+step);
    randomResult=c(randomResult,randomVariable)  
  }
  
  randomResult
}

#Вычисление точности моделирования с помощбю критерия хи квадрат Пирсона
accyracy <- function(density, sample, mean,sd,numberOfCuts=1000,deviationFromMeanBound=8,lower=NULL,upper=NULL)
{  
  if(!is.null(upper))
  {
    upperBound=upper;
  }
  else
  {
    upperBound=deviationFromMeanBound*sd+mean;
  }
  if(!is.null(lower))
  {
    lowerBound=lower;
  }
  else
  {
    lowerBound=-deviationFromMeanBound*sd+mean;
  }
  step=(upperBound-lowerBound)/numberOfCuts;
  
  #Находим Вероятность попадания эмпирических величин на каждый из отрезков
  
  
  pEmp=c();
  for(currentCut in 1:numberOfCuts)
  {
    pEmp=c(pEmp,0);
  }
  for(curentX in 1:length(sample))
  {
    cutIndex=(sample[curentX]-lowerBound)/step;
    cutIndex=cutIndex+1;
    if(sample[curentX]>=lowerBound && sample[curentX]<upperBound)
      pEmp[cutIndex]=pEmp[cutIndex]+1;
  }
  pEmp=pEmp / length(sample);
 
  #Находим Вероятность попадания случайных величин соответствующих этому распределению 
  #на каждый из отрезков, интеглалл на каждом вычисляем по формуле Симпсона
  
  pH0=c();
  for(currentCut in 1:numberOfCuts)
  {
    cutBegin=lowerBound+step*(currentCut-1);    
    cutEnd=cutBegin+step;   
    integrall=(step/3) * (density(cutBegin)+4*density((cutBegin+cutEnd)/2)+density(cutEnd));
    pH0=c(pH0, integrall);    
  }     
  xiSquared=sum((((pH0-pEmp)^2))/pH0,na.rm=TRUE)*length(sample);
  
  
  return (pchisq(xiSquared,length(sample)-1));
}

#Вычисление прараметров модели при помощи метода максимального правдоподобия
findParams <- function(sample,initParams,densityOfXAndvectorOfParams,...)
{
  print("start")
  fn<-function (params)
  {   
    res=0;
    for(i in 1:length(sample))
    {
      d=try(densityOfXAndvectorOfParams(sample[i],params));      
      if(!is.numeric(d) || d<=0)
      {
        res=res-100000000000;
      }
      else
      {
        res=res+log(d);
      }
    }  
   return (-res);
  }
  return (optim(par=initParams,fn=fn,...))
}






