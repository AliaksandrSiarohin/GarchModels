I<-complex(real=0,imaginary=1);
cr <- function(u,a,m,C,l1,l2)
{
  (sqrt(pi)*C*gamma(-a/2)/(2^((a+3)/2)))*((l1*l1+u*u)^(a/2)-l1^a+(l2*l2+u*u)^(a/2)-l2^a);
}
cl <- function(u,a,m,C,l1,l2)
{
  (I*u*C*gamma(-a/2)/(2^((a+3)/2)))*((l1^(a-1))*hypergeo(1,(1-a)/2,3/2,-u*u/(l1*l1))-
    (l2^(a-1))*hypergeo(1,(1-a)/2,3/2,-u*u/(l2*l2)));
}
charFunc <- function(u,a,m,C,l1,l2)
{
 
  exp(I*u*m + cr(u,a,m,C,l1,l2)+cl(u,a,m,C,l1,l2));
}
density <-function(x,a,m,C,l1,l2)
{  
  
  I<-complex(real=0,imaginary=1);
  param <- function(u)
  {
    charFunc(u,a,m,C,l1,l2)*exp(-I*x*u)/(2*pi);
  }
  Re(myintegrate(param,lower=-Inf,upper=Inf));
}
getFx <- function(x,a,m,C,l1,l2)
{
  fx<-c();
  
  for(i in 1:length(x))
  {
    fx<-c(fx,density(x[i],a,m,C,l1,l2));
  }
  print (fx);
}
x <- seq(-4, 3,length=100)
a <- c(0.5,0.5,0.5,0.5)
m <- c(1,1,1,1)
l1<- c(1,1,1,1)
l2<- c(1,1,1,1)
C<-c(0.5,1,6,10)

colors <- c("red", "blue", "darkgreen", "gold")

plot(x,
     getFx(x,a[1],m[1],C[1],l1[1],l2[1]),
     lwd=2, 
     col=colors[1],
     type="l",
     xlab="x value",
     
     ylab="Density",
     main="Comprassion Of MTS Distributions")

for (i in 2:length(m))
{
  lines(x, getFx(x,a[i],m[i],C[i],l1[i],l2[i]),type="l",
        lwd=2, col=colors[i])
}
labels=c()
for (i in 1:length(m))
{
  labels <- c(labels,paste(sep="",
                           "a = ",as.character(a[i]),
                           ", m = ",as.character(m[i]),
                           ", l+ = ",as.character(l1[i]),
                           ", l- = ",as.character(l2[i]),
                           ", C = ",as.character(C[i])));
}
legend("topleft", inset=.05, title="Distributions",
       legend=labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
