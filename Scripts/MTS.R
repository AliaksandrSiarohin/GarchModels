require("elliptic");
require("hypergeo");


I<-complex(real=0,imaginary=1);
aDefault  <- 0.3
mDefault  <- 10
l1Default <- 2
l2Default <- 2
CDefault  <- 2


cr <- function(u,a,m,C,l1,l2)
{
  (sqrt(pi)*C*gamma(-a/2)/(2^((a+3)/2)))*((l1*l1+u*u)^(a/2)-l1^a+(l2*l2+u*u)^(a/2)-l2^a);
}
cl <- function(u,a,m,C,l1,l2)
{
  (I*u*C*gamma(-a/2)/(2^((a+3)/2)))*((l1^(a-1))*hypergeo(1,(1-a)/2,3/2,-u*u/(l1*l1))-
    (l2^(a-1))*hypergeo(1,(1-a)/2,3/2,-u*u/(l2*l2)));
}


#Характеристическая функция
cMTS <- function(u,
                  a=aDefault,
                  m=mDefault,
                  C=CDefault,
                  l1=l1Default,
                  l2=l2Default)
{
  
  exp(I*u*m + cr(u,a,m,C,l1,l2)+cl(u,a,m,C,l1,l2));
}

#Плотность вычисленная из харарктеристической функции
dMTS<-function(x,
               a=aDefault,
               m=mDefault,
               C=CDefault,
               l1=l1Default,
               l2=l2Default,...)
{  
  
  I<-complex(real=0,imaginary=1);
  param <- function(u)
  {
    cMTS(u,a,m,C,l1,l2)*exp(-I*x*u)/(2*pi);
  }
  Re(myintegrate(param,lower=-Inf,upper=Inf,rel.tol=.Machine$double.eps^0.5,...));
}

cumMTS <- function (n,
                    a=aDefault,
                    m=mDefault,
                    C=CDefault,
                    l1=l1Default,
                    l2=l2Default)
{
  
  if (n %% 2 == 1)
  {
    cum=0;
    if (n==1)
    {
      cum=m;
    }
    cum=cum+(2^(n-(a+3)/2)) * 
      factorial((n-1)/2) * 
      C *
      gamma((n-a)/2) *
      (l1^(a-n)-l2^(a-n));
    return (cum);
  }
  else
  {
    cum=2^(-(a+3)/2) *
        sqrt(pi) *
        ( factorial(n) / factorial(n/2) ) *
        C *
        gamma((n-a)/2) *
        (l1^(a-n)+l2^(a-n))
    return (cum);
  }
  
}


#Математическое ожидание
meanMTS <- function (a=aDefault,
                     m=mDefault,
                     C=CDefault,
                     l1=l1Default,
                     l2=l2Default)
{
  cumMTS(1,a,m,C,l1,l2);
}

#Дисперсия

varMTS <-function (a=aDefault,
                   m=mDefault,
                   C=CDefault,
                   l1=l1Default,
                   l2=l2Default)
{
  cumMTS(2,a,m,C,l1,l2);
}

#Ассиметрия 
skewnessMTS <-function (a=aDefault,
                        m=mDefault,
                        C=CDefault,
                        l1=l1Default,
                        l2=l2Default)
{
  cumMTS(3,a,m,C,l1,l2)/ ((cumMTS(2,a,m,C,l1,l2))^1.5)
}

#Куртосис
kurtosisMTS <-function (a=aDefault,
                        m=mDefault,
                        C=CDefault,
                        l1=l1Default,
                        l2=l2Default)
{
  cumMTS(4,a,m,C,l1,l2)/ ((cumMTS(2,a,m,C,l1,l2))^2)
}


#Моделирование Случайных Величин MTS
rMTS<- function (n,
                 a=aDefault,
                 m=mDefault,
                 C=CDefault,
                 l1=l1Default,
                 l2=l2Default)
{
  density<- function(x)
  {
    dMTS(x,a,m,C,l1,l2);
  }
  model(density,n,meanMTS(a,m,C,l1,l2),sqrt(varMTS(a,m,C,l1,l2)));
}

accyracyMTS <- function(sample,
                        a=aDefault,
                        m=mDefault,
                        C=CDefault,
                        l1=l1Default,
                        l2=l2Default)
{
  density<- function(x)
  {
    dMTS(x,a,m,C,l1,l2);
  }
  accyracy(density,sample,meanMTS(a,m,C,l1,l2),sqrt(varMTS(a,m,C,l1,l2)));
}

printTableMTS <- function()
{
  aVal  <- c(0.3,1.1,0.3,0.3,0.3)
  mVal  <- c(10,0,0,0,0)
  l1Val <- c(2,2,0.1,2,2)
  l2Val <- c(2,2,2,20,2)
  CVal  <- c(2,2,2,2,10)
  for (i in 1:length(aVal))
  {
    aDefault  <- aVal[i]
    mDefault  <- mVal[i]
    l1Default <- l1Val[i]
    l2Default <- l2Val[i]
    CDefault  <- CVal[i]
    samples<-rMTS(1000,aVal[i],mVal[i],CVal[i],l1Val[i],l2Val[i]);
    
    print(paste("mean : ",as.character(meanMTS(aVal[i],mVal[i],CVal[i],l1Val[i],l2Val[i])),mean(samples)));
    print(paste("variation : ",as.character(varMTS(aVal[i],mVal[i],CVal[i],l1Val[i],l2Val[i])),var(samples)));
    print(paste("skewness : ",as.character(skewnessMTS(aVal[i],mVal[i],CVal[i],l1Val[i],l2Val[i])),skewness(samples)));
    print(paste("kurtosis : ",as.character(kurtosisMTS(aVal[i],mVal[i],CVal[i],l1Val[i],l2Val[i])),kurtosis(samples)));
    print((paste("accyracy : ",accyracyMTS(samples,aVal[i],mVal[i],CVal[i],l1Val[i],l2Val[i]))));
  }
}

printParamsTableMTS <- function()
{
  n<-c(5000);
  initParams=c(aDefault,mDefault,CDefault,l1Default,l2Default);
  print("True params :");
  print(initParams);  
  for(i in 1:length(n))
  {
    samples<-rMTS(n[i]);
    
    params<-findParams(samples,initParams,function(x,params){
      dMTS(x,params[1],params[2],params[3],params[4],params[5])},control = list(maxit=100));
    print(paste("Finded params for ",n[i]," samples"));
    print(params);
  }
}

stdMTS <- function(x,params)
{
  mean=meanMTS(params[1],params[2],params[3],params[4],params[5]);
  sd=sqrt(abs(varMTS(params[1],params[2],params[3],params[4],params[5])));
  density<-linearDensityTransform(function(x){ 
    dMTS(x,params[1],params[2],params[3],params[4],params[5]) },-mean/sd,1/sd);
  return (density(x));
}


modelSTDFunctionMTS<-function(n, params)
{
  density<- function(x)
  {
    stdMTS(x,params);
  }
  return (model(density,n,0,1));  
}






