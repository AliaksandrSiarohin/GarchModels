BiG_NUMBER=1000000000000000000000000;

#ÐÑÑÐ¸ÑÐ»ÐµÐ½Ð¸Ðµ Ð¿Ð»Ð¾ÑÐ½Ð¾ÑÑÐ¸ y=a+b*x 
linearDensityTransform<-function(density,a,b)
{
  result<-function(x)
  {
    return (density((x-a)/b)/abs(b));
  }
  return (result);
}

#ÐÑÐµÐ¾Ð±ÑÐ°Ð·Ð¾Ð²Ð°Ð½Ð¸Ðµ ÑÐµÐ½ Ðº log returns
transformPrices<-function(prices)
{
  logRet=c();
  for(i in 2:length(prices))
  {
    logRet=c(logRet,log(prices[i]/prices[i-1]));
  }
  return (logRet);
}

computeNextSigma<-function(alpha0,alpha,beta,p,q,sigma,h)
{
  nextSigma=alpha0;
  n=length(sigma)+1;
  C=0;
  for(i in 1:p)
  {
    if((n-i) > 0)
    {
      nextSigma=nextSigma+alpha[i] * (h[n-i] - C)  * (h[n-i] - C);
    }          
  }
  for(i in 1:q)
  {
    if((n-i) > 0)
    {
      nextSigma=nextSigma+beta[i] * sigma[n-i] * sigma[n-i];
    }          
  }
  return (sqrt(nextSigma));
}


#ÐÐ¾Ð´Ð¾Ð±ÑÐ°ÑÑ Ð¿Ð°ÑÐ°Ð¼ÐµÑÑÑ Garch Ð¼Ð¾Ð´ÐµÐ»Ð¸ Ð²Ð¸Ð´Ð° h_n= sigma_n * epsilon_n Ð³Ð´Ðµ 
#sigma_n ^ 2 = alpha_0 + sum(alpha_i * h_(n-i) ^ 2) + sum(beta_i * sigma_(n-i) ^ 2) 
#p and q model order
#epsilon_n std distribution
fitGarchModel<-function(logReturns,stdDensity,stdDensityInitParams,p=1,q=1,...)
{
  initParams=c(stdDensityInitParams,var(logReturns),rep(0.1,length.out=p),rep(0.86,length.out=q));
 
  joinDensity <- function(params)
  {
    eps_params<-params[1:length(stdDensityInitParams)];
    alpha0=params[length(stdDensityInitParams)+1];
    alpha=params[(length(stdDensityInitParams)+2) : (length(stdDensityInitParams)+1+p)];
    beta=params[(length(stdDensityInitParams)+p+2) : (length(params)-1)];  
    C=0;
    
    if(alpha0<0)
      return (BiG_NUMBER);
    if(length(alpha[alpha<0])!=0)
      return (BiG_NUMBER);
    if(length(beta[beta<0])!=0)
      return (BiG_NUMBER);
    
    sigma<-c();
    epsilon<-c();
    res<-0;
    
    for(n in 1:length(logReturns))
    {
      sigma[n]=computeNextSigma(alpha0,alpha,beta,p,q,sigma,logReturns);     
    
      
      epsilon[n]=(logReturns[n]-C)/sigma[n];      
      
      density<-function(x){return (stdDensity(x,eps_params));};
      h_density=linearDensityTransform(density,C,sigma[n]);
      
      d=try(h_density(logReturns[n]));      
      if(n!=1)
      {
        if(is.nan(d) || !is.numeric(d) || d<=0)
        {
          return (BiG_NUMBER);
        }
        else
        {
          res=res+log(d);
        }
      }
      
    } 
    return (-res);
  }
  
  return (optim(par=initParams,fn=joinDensity,
            #    lower=c(rep(-Inf,length(stdDensityInitParams)),0,rep(0.0001,length.out=p),rep(0.0001,length.out=q),-Inf),
         # upper=c(rep(Inf,length(stdDensityInitParams)),1,rep(0.9999,length.out=p),rep(0.9999,length.out=q),Inf),
                ))
}
  
garchModel<-function(sampleLength,modelFunction,params,numberOfStdDensityParams,p=1,q=1)
{
  eps_params<-params[1:numberOfStdDensityParams];
  alpha0=params[numberOfStdDensityParams+1];
  alpha=params[(numberOfStdDensityParams+2) : (numberOfStdDensityParams+1+p)];
  beta=params[(numberOfStdDensityParams+p+2) : (length(params)-1)];
  C=0;
  
  h<-c();
  sigma<-c();
  
  epsilon<-modelFunction(sampleLength,eps_params);
  for(n in 1:sampleLength)
  {
    sigma[n]=computeNextSigma(alpha0,alpha,beta,p,q,sigma,h);     
    h[n]=C+sigma[n]*epsilon[n];  
  }
  return (h);
  
}

garchForecast<-function(logRet,sampleLength,modelFunction,params,numberOfStdDensityParams,p=1,q=1)
{
  eps_params<-params[1:numberOfStdDensityParams];
  alpha0=params[numberOfStdDensityParams+1];
  alpha=params[(numberOfStdDensityParams+2) : (numberOfStdDensityParams+1+p)];
  beta=params[(numberOfStdDensityParams+p+2) : (length(params)-1)];
  C=0;
  
  h<-c();
  sigma<-c();  
  
  for(n in 1:length(logRet))
  {   
    sigma[n]=computeNextSigma(alpha0,alpha,beta,p,q,sigma,h);         
    h[n]=logRet[n];  
  }
  
  epsilon<-modelFunction(sampleLength,eps_params);
  for(n in (length(logRet) + 1) :(length(logRet)+sampleLength))
  {
    sigma[n]=computeNextSigma(alpha0,alpha,beta,p,q,sigma,h);
    h[n]=C+sigma[n]*epsilon[n - length(logRet)];
  }
  return (h[(length(logRet) + 1) :(length(logRet)+sampleLength)]);
}



stdNormal<-function(x,params)
{
  mean=params[1];
  sd=abs(params[2]);
  density<-linearDensityTransform(function(x){ dnorm(x,params[1],abs(params[2])) },-mean/sd,1/sd);
  return (density(x));
}

variation <- function(alpha0, alpha1, beta1) {
  alpha0/(1-alpha1-beta1);
}

kurtosis <-function(alpha0, alpha1, beta1, kur_e) {
  kur_e*(alpha0^2+alpha0^2(alpha1+beta1))/((1-alpha1-beta1)*(1-beta1^2-alpah1^2*kur_e-2beta1*alpha1))  
}

kur_e <- function(params, var, kur) {
  kur(params)/ (var(params)^2);
}
  
  
  
  