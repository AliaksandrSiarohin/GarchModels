require("elliptic");

I<-complex(real=0,imaginary=1);
aDefault  <- 0.5
mDefault  <- 0
l1Default <- 1
l2Default <- 1
C1Default <- 0.5
C2Default <- 3


#Характеристическая функция
cCTS <- function(u,
                     a=aDefault,
                     m=mDefault,
                     C1=C1Default,
                     C2=C2Default,
                     l1=l1Default,
                     l2=l2Default)
{ 
  
    exp(I*u*m + C1*gamma(-a)*((l1 - I*u)^a - l1^a) + 
      C2*gamma(-a)*((l2 + I*u)^a - l2^a));
  
}

#Плотность вычисленная из харарктеристической функции
dCTS <- function(x,
                 a=aDefault,
                 m=mDefault,
                 C1=C1Default,
                 C2=C2Default,
                 l1=l1Default,
                 l2=l2Default)
{ 
  param <- function(u)
  {
   cCTS(u,a,m,C1,C2,l1,l2)*exp(-I*x*u)/(2*pi);   
  }
  Re(myintegrate(param,lower=-Inf,upper=Inf,rel.tol=.Machine$double.eps^0.5));
}

cumCTS <- function (n,
                    a=aDefault,
                    m=mDefault,
                    C1=C1Default,
                    C2=C2Default,
                    l1=l1Default,
                    l2=l2Default)
{
    
    if (n==1)
    {
      return (m+gamma(1-a)*(C1*l1^(a-1)-C2*l2^(a-1)));
    }
    else
    {
      return (gamma(n-a)*(C1*l1^(a-n)+(-1)^n*C2*l2^(a-n)));
    }
    
}


#Математическое ожидание
meanCTS <- function (a=aDefault,
                     m=mDefault,
                     C1=C1Default,
                     C2=C2Default,
                     l1=l1Default,
                     l2=l2Default)
{
   cumCTS(1,a,m,C1,C2,l1,l2);
}

#Дисперсия

varCTS <-function (a=aDefault,
                   m=mDefault,
                   C1=C1Default,
                   C2=C2Default,
                   l1=l1Default,
                   l2=l2Default)
{
  cumCTS(2,a,m,C1,C2,l1,l2);
}

#Ассиметрия 
skewnessCTS <-function (a=aDefault,
                        m=mDefault,
                        C1=C1Default,
                        C2=C2Default,
                        l1=l1Default,
                        l2=l2Default)
{
  cumCTS(3,a,m,C1,C2,l1,l2)/ ((cumCTS(2,a,m,C1,C2,l1,l2))^1.5)
}

#Куртосис
kurtosisCTS <-function (a=aDefault,
                        m=mDefault,
                        C1=C1Default,
                        C2=C2Default,
                        l1=l1Default,
                        l2=l2Default)
{
  cumCTS(4,a,m,C1,C2,l1,l2)/ ((cumCTS(2,a,m,C1,C2,l1,l2))^2)
}


#Моделирование Случайных Величин CTS
rCTS<- function (n,
                 a=aDefault,
                 m=mDefault,
                 C1=C1Default,
                 C2=C2Default,
                 l1=l1Default,
                 l2=l2Default)
{
  density<- function(x)
  {
    dCTS(x,a,m,C1,C2,l1,l2);
  }
  model(density,n,meanCTS(a,m,C1,C2,l1,l2),sqrt(varCTS(a,m,C1,C2,l1,l2)));
}

accyracyCTS <- function(sample,
                        a=aDefault,
                        m=mDefault,
                        C1=C1Default,
                        C2=C2Default,
                        l1=l1Default,
                        l2=l2Default)
{
  density<- function(x)
  {    
    dCTS(x,a,m,C1,C2,l1,l2);
  }
  accyracy(density,sample,meanCTS(a,m,C1,C2,l1,l2),sqrt(varCTS(a,m,C1,C2,l1,l2)));
}

stdCTS <- function(x,params)
{
  mean=meanCTS(params[1],params[2],params[3],params[4],params[5],params[6]);
  sd=sqrt(abs(varCTS(params[1],params[2],params[3],params[4],params[5],params[6])));
  density<-linearDensityTransform(function(x){ 
    dCTS(x,params[1],params[2],params[3],params[4],params[5],params[6]
         ) },-mean/sd,1/sd);
  return (density(x));
}


printTableCTS <- function()
{
  aVal  <- c(0.5,0.5,1.9,0.5,0.5,0.5,0.5)
  mVal  <- c(0,  10, 0,  0,  0,  0,  0  )
  l1Val <- c(1,  1,  1,  10, 1,  1,  1  )
  l2Val <- c(1,  1,  1,  1,  10, 1,  1  )
  C1Val  <-c(0.5,0.5,0.5,0.5,0.5,0.1,0.5)
  C2Val  <-c(3,  3,  3,  3,  3,  3,  30 )
  
  for(i in 1:length(aVal))
  {
    samples<-rCTS(1000,aVal[i],mVal[i],C1Val[i],C2Val[i],l1Val[i],l2Val[i]);
    print(paste("mean : ",as.character(meanCTS(aVal[i],mVal[i],C1Val[i],C2Val[i],l1Val[i],l2Val[i])),
                mean(samples)));
    print(paste("variation : ",as.character(varCTS(aVal[i],mVal[i],C1Val[i],C2Val[i],l1Val[i],l2Val[i])),
                var(samples)));
    print(paste("skewness : ",as.character(skewnessCTS(aVal[i],mVal[i],C1Val[i],C2Val[i],l1Val[i],l2Val[i])),
                skewness(samples)));
    print(paste("kurtosis : ",as.character(kurtosisCTS(aVal[i],mVal[i],C1Val[i],C2Val[i],l1Val[i],l2Val[i])),
                kurtosis(samples)));
    
    print((paste("accyracy : ",accyracyCTS(samples,aVal[i],mVal[i],C1Val[i],C2Val[i],l1Val[i],l2Val[i]))));
  }
}


printParamsTableCTS <- function()
{
  n<-c(100,500,1000,5000);
  initParams=c(aDefault,mDefault,C1Default,C2Default,l1Default,l2Default);
  print("True params :");
  print(initParams);  
  for(i in 1:length(n))
  {
    samples<-rCTS(n[i]);
    
    params<-findParams(samples,initParams,function(x,params){
      dCTS(x,params[1],params[2],params[3],params[4],params[5],params[6])});
    print(paste("Finded params for ",n[i]," samples"));
    print(params);
  }
}


modelSTDFunctionCTS<-function(n, params)
{
  density<- function(x)
  {
    stdCTS(x,params);
  }
  return (model(density,n,0,1));  
}

