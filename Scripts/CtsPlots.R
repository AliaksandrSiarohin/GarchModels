
charFunc <- function(u,a,m,C1,C2,l1,l2)
{
  I<-complex(real=0,imaginary=1);
  exp(I*u*m + C1*gamma(-a)*((l1 - I*u)^a - l1^a) + 
    C2*gamma(-a)*((l2 + I*u)^a - l2^a));
}
density <-function(x,a,m,C1,C2,l1,l2)
{
  

  I<-complex(real=0,imaginary=1);
  param <- function(u)
  {
    charFunc(u,a,m,C1,C2,l1,l2)*exp(-I*x*u)/(2*pi);
  }
  Re(myintegrate(param,lower=-Inf,upper=Inf));
}
getFx <- function(x,a,m,C1,C2,l1,l2)
{
  fx<-c();
  
  for(i in 1:length(x))
  {
    fx<-c(fx,density(x[i],a,m,C1,C2,l1,l2));
  }
  print (fx);
}
x <- seq(-10, 30,length=500)
a <- c(0.5,0.5,0.5,0.5)
m <- c(0,1,2,3)
l1<- c(1,1,1,1)
l2<- c(1,1,1,1)
C1<-c(0.5,2.0,3.0,0.5)
C2<-c(3.0,3.0,2.0,4.0)

colors <- c("red", "blue", "darkgreen", "gold")

plot(x,
     getFx(x,a[1],m[1],C1[1],C2[1],l1[1],l2[1]),
     lwd=2, 
     col=colors[1],
     type="l",
     xlab="x value",
     ylab="Density",
     main="Comprassion Of CTS Distributions")

for (i in 2:length(m))
{
  lines(x, getFx(x,a[i],m[i],C1[i],C2[i],l1[i],l2[i]),type="l",
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
                           ", C1 = ",as.character(C1[i]),
                           ", C2 = ",as.character(C2[i])));
}
legend("topright", inset=.05, title="Distributions",
       legend=labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
