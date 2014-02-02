require('moments');
stats <- function(X) {
  
  c(mean(X), var(X),  skewness(X), kurtosis(X));
}

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
    prob=integrate(density,cutBeginings[currentCut],cutBeginings[currentCut]+step)$value;
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

#Проверяем гипотезу согласия Xи квадрат Пирсона
X_hipothes <- function(density, sample, deviationFromMeanBound=10,lower=NULL,upper=NULL,
                       x_per_interval=5,numberOfCuts=1000)
{  
  mean= mean(sample);
  sd=sqrt(var(sample));  
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
  counts = pEmp;
  pEmp=pEmp / length(sample);
  
  #Находим Вероятность попадания случайных величин соответствующих этому распределению 
  #на каждый из отрезков, интеглалл на каждом вычисляем по формуле Симпсона
  
  pH0=c();
  for(currentCut in 1:numberOfCuts)
  {
    cutBegin=lowerBound+step*(currentCut-1);    
    cutEnd=cutBegin+step;   
    integrall=integrate(density,cutBegin,cutEnd)$value;
    pH0=c(pH0, integrall);    
  }  
  new_pH0=c();
  new_pEmp=c();
  new_counts=c();
  pre_count = 0; 
  pre_h = 0;
  pre_emp = 0;
  for(i in 1:length(counts))
  {
    pre_count = pre_count+counts[i];
    pre_h = pre_h+pH0[i];
    pre_emp = pre_emp+pEmp[i];
    if ( pre_count >= x_per_interval )
    {
      new_pH0=c(new_pH0,pre_h);
      new_pEmp=c(new_pEmp,pre_emp);      
      pre_count = 0;
      pre_emp=0;
      pre_h=0;
      
    }    
    
  }
  
  pH0=new_pH0;
  pEmp=new_pEmp;  
  numberOfCuts=length(pH0);
  xiSquared=sum((((pH0-pEmp)^2))/pH0,na.rm=TRUE)*length(sample); 
  
  return (pchisq(xiSquared,numberOfCuts-1));
}

optim_X_hipo <-function(alpha, interval=c(10,1000),...)
{
  return (optimize(function(n_c){X_hipothes(numberOfCuts=n_c,...)},interval)$objective)
}

#Вычисление прараметров модели при помощи метода максимального правдоподобия
findParams <- function(sample,initParams,densityOfXAndvectorOfParams,...)
{  
  fn<-function (params)
  {   
    res=0;    
    print(params);
    for(i in 1:length(sample))
    {
      d=try(suppressWarnings(densityOfXAndvectorOfParams(sample[i],params)));   
      
      if(is.na(d) || !is.numeric(d) || d<=0 || d == Inf || d ==-Inf)
      {
        res=res-100000000000;
      }
      else
      {
        res=res+log(d);
      }
    }  
    print(res);
    return (-res);
  }
  return (optim(par=initParams,fn=fn,...))
}


confidence_interval <- function(theta, alpha) {
  theta_hat=mean(theta);
  delta = qnorm(1-alpha)*sqrt(var(theta))/sqrt(length(theta)-1);
  return (c(mean(theta)-delta,mean(theta)+delta));
}


