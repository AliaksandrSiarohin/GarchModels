I<-complex(real=0,imaginary=1);

H <- function(u,a,k,r,p)
{
  (k*gamma(-a)/p)*(hypergeo(p,-a,1+p,I*r*u)-1);
}
charFunc <- function(u,a,m,k1,k2,r1,r2,p1,p2)
{
  exprTail=I*u*a*gamma(-a)*((k1*r1)/(p1+1)-(k2*r2)/(p2+1))
  exp(I*u*m + H(u,a,k1,r1,p1)+H(-u,a,k2,r2,p2)+ exprTail);
}
density <-function(x,a,m,k1,k2,r1,r2,p1,p2)
{  
  
  I<-complex(real=0,imaginary=1);
  param <- function(u)
  {
    charFunc(u,a,m,k1,k2,r1,r2,p1,p2)*exp(-I*x*u)/(2*pi);
  }
  Re(myintegrate(param,lower=-Inf,upper=Inf));
}
getFx <- function(x,a,m,k1,k2,r1,r2,p1,p2)
{
  fx<-c();
  
  for(i in 1:length(x))
  {
    fx<-c(fx,density(x[i],a,m,k1,k2,r1,r2,p1,p2));
  }
  print (fx);
}
x <- seq(-6, 4,length=100)
a <- c(0.5,0.5,0.5,0.5)
m <- c(1.5,1.5,1.5,1.5)
k1<- c(1,1,1,1)
k2<- c(1,1,1,1)
p1<- c(1,1,1,1)
p2<- c(1,1,1,1)
r1<- c(-0.3,1,0.1,2)
r2<- c(-0.3,0.1,1,2)

colors <- c("red", "blue", "darkgreen", "gold")

plot(x,
     getFx(x,a[1],m[1],k1[1],k2[1],r1[1],r2[1],p1[1],p2[1]),
     lwd=2, 
     col=colors[1],
     type="l",
     xlab="x value",
     
     ylab="Density",
     main="Comprassion Of KR Distributions")

for (i in 2:length(m))
{
  lines(x, getFx(x,a[i],m[i],k1[i],k2[i],r1[i],r2[i],p1[i],p2[i]),type="l",
        lwd=2, col=colors[i])
}
labels=c()
for (i in 1:length(m))
{
  labels <- c(labels,paste(sep="",
                           "a = ",as.character(a[i]),
                           ", m = ",as.character(m[i]),
                           ", k+ = ",as.character(k1[i]),
                           ", k- = ",as.character(k2[i]),
                           ", r+ = ",as.character(r1[i]),
                           ", r- = ",as.character(r2[i]),
                           ", p+ = ",as.character(p1[i]),
                           ", p- = ",as.character(p2[i])));
}
legend("topleft", inset=.05, title="Distributions",
       legend=labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
