require("elliptic");
require("hypergeo");

I<-complex(real=0,imaginary=1);

aDefault  <- 0.5
mDefault  <- 0
k1Default <- 1
k2Default <- 1
r1Default <- 1
r2Default <- 1
p1Default <- 1
p2Default <- 1


H <- function(u,a,k,r,p)
{
  (k*gamma(-a)/p)*(hypergeo(p,-a,1+p,I*r*u)-1);
}

#Характеристическая функция
cKR <- function(u,
                a=aDefault,
                m=mDefault,
                k1=k1Default,
                k2=k2Default,
                r1=r1Default,
                r2=r2Default,
                p1=p1Default,
                p2=p2Default)
{
  exprTail=I*u*a*gamma(-a)*((k1*r1)/(p1+1)-(k2*r2)/(p2+1))
  exp(I*u*m + H(u,a,k1,r1,p1)+H(-u,a,k2,r2,p2)+ exprTail);
}

#Плотность вычисленная из харарктеристической функции
dKR<-function(x, 
              a=aDefault,
              m=mDefault,
              k1=k1Default,
              k2=k2Default,
              r1=r1Default,
              r2=r2Default,
              p1=p1Default,
              p2=p2Default)
{
  param <- function(u)
  {
    cKR(u,a,m,k1,k2,r1,r2,p1,p2)*exp(-I*x*u)/(2*pi);
  }
  Re(myintegrate(param,lower=-Inf,upper=Inf,rel.tol=.Machine$double.eps^0.5));
}

cumKR <- function (n,
                   a=aDefault,
                   m=mDefault,
                   k1=k1Default,
                   k2=k2Default,
                   r1=r1Default,
                   r2=r2Default,
                   p1=p1Default,
                   p2=p2Default)
{
  
  if(n==1)
  {
    return (m);
  }
  else
  {
    return (gamma(n-a) *
      ((k1*r1^n)/(p1+n)+((-1)^n)*(k2*r2^n)/(p2+n)));
  }
  
}


#Математическое ожидание
meanKR <- function (a=aDefault,
                    m=mDefault,
                    k1=k1Default,
                    k2=k2Default,
                    r1=r1Default,
                    r2=r2Default,
                    p1=p1Default,
                    p2=p2Default)
{
  cumKR(1,a,m,k1,k2,r1,r2,p1,p2);
}

#Дисперсия

varKR <-function (a=aDefault,
                  m=mDefault,
                  k1=k1Default,
                  k2=k2Default,
                  r1=r1Default,
                  r2=r2Default,
                  p1=p1Default,
                  p2=p2Default)
{
  cumKR(2,a,m,k1,k2,r1,r2,p1,p2);
}

#Ассиметрия 
skewnessKR <-function (a=aDefault,
                       m=mDefault,
                       k1=k1Default,
                       k2=k2Default,
                       r1=r1Default,
                       r2=r2Default,
                       p1=p1Default,
                       p2=p2Default)
{
  cumKR(3,a,m,k1,k2,r1,r2,p1,p2)/ ((cumKR(2,a,m,k1,k2,r1,r2,p1,p2))^1.5)
}

#Куртосис
kurtosisKR <-function (a=aDefault,
                       m=mDefault,
                       k1=k1Default,
                       k2=k2Default,
                       r1=r1Default,
                       r2=r2Default,
                       p1=p1Default,
                       p2=p2Default)
{
  cumKR(4,a,m,k1,k2,r1,r2,p1,p2)/ ((cumKR(2,a,m,k1,k2,r1,r2,p1,p2))^2)
}


#Моделирование Случайных Величин KR
rKR<- function (n,
                a=aDefault,
                m=mDefault,
                k1=k1Default,
                k2=k2Default,
                r1=r1Default,
                r2=r2Default,
                p1=p1Default,
                p2=p2Default)
{
  density<- function(x)
  {
    dKR(x,a,m,k1,k2,r1,r2,p1,p2);
  }
  model(density,n,meanKR(a,m,k1,k2,r1,r2,p1,p2),sqrt(varKR(a,m,k1,k2,r1,r2,p1,p2)));
}

accyracyKR <- function(sample,
                       a=aDefault,
                       m=mDefault,
                       k1=k1Default,
                       k2=k2Default,
                       r1=r1Default,
                       r2=r2Default,
                       p1=p1Default,
                       p2=p2Default)
{
  density<- function(x)
  {    
    dKR(x,a,m,k1,k2,r1,r2,p1,p2);
  }
  accyracy(density,sample,meanKR(a,m,k1,k2,r1,r2,p1,p2),sqrt(varKR(a,m,k1,k2,r1,r2,p1,p2)));
}

printTableKR <- function()
{
  aVal  <- c(0.5,1.9,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
  mVal  <- c(0,  0,  10, 0,  0,  0,  0,  0,  0  )
  k1Val <- c(1,  1,  1,  10, 1,  1,  1,  1,  1  )
  k2Val <- c(1,  1,  1,  1,  10, 1,  1,  1,  1  )
  r1Val <- c(1,  1,  1,  1,  1,  10, 1,  1,  1  )
  r2Val <- c(1,  1,  1,  1,  1,  1,  10, 1,  1  )
  p1Val <- c(1,  1,  1,  1,  1,  1,  1,  10, 1  )
  p2Val <- c(1,  1,  1,  1,  1,  1,  1,  1,  10 )
  for (i in 1:length(aVal))
  {   
    samples<-rKR(1000,aVal[i],mVal[i],k1Val[i],k2Val[i],r1Val[i],r2Val[i],p1Val[i],p2Val[i]);
    print(paste("mean : ",as.character(meanKR(aVal[i],mVal[i],k1Val[i],k2Val[i],r1Val[i],r2Val[i],p1Val[i],p2Val[i]))
                ,mean(samples)));
    print(paste("variation : ",as.character(varKR(aVal[i],mVal[i],k1Val[i],k2Val[i],r1Val[i],r2Val[i],p1Val[i],p2Val[i])),
                var(samples)));
    print(paste("skewness : ",as.character(skewnessKR(aVal[i],mVal[i],k1Val[i],k2Val[i],r1Val[i],r2Val[i],p1Val[i],p2Val[i])),
                skewness(samples)));
    print(paste("kurtosis : ",as.character(kurtosisKR(aVal[i],mVal[i],k1Val[i],k2Val[i],r1Val[i],r2Val[i],p1Val[i],p2Val[i])),
                kurtosis(samples)));
    print(paste("accyracy : ",accyracyKR(samples,aVal[i],mVal[i],k1Val[i],k2Val[i],r1Val[i],r2Val[i],p1Val[i],p2Val[i])));
    print("----------------------------------------------");
  }
}


printParamsTableKR <- function()
{
  n<-c(100,500,1000,5000);
  initParams=c(aDefault,mDefault,k1Default,k2Default,r1Default,r2Default,p1Default,p2Default);
  print("True params :");
  print(initParams);  
  for(i in 1:length(n))
  {
    samples<-rKR(n[i]);
    
    params<-findParams(samples,initParams,function(x,params){
      dKR(x,params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8])},control = list(maxit=100));
    print(paste("Finded params for ",n[i]," samples"));
    print(params);
  }
}


stdKR <- function(x,params)
{
  mean=meanKR(params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8]);
  sd=sqrt(abs(varKR(params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8])));
  density<-linearDensityTransform(function(x){ 
    dKR(x,params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8]) },-mean/sd,1/sd);
  return (density(x));
}

modelSTDFunctionKR<-function(n, params)
{
  density<- function(x)
  {
    stdKR(x,params);
  }
  return (model(density,n,0,1));  
}