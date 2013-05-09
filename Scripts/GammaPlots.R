

x <- seq(0, 10,length=100)


theta <- c(2,2,2,2)
k<-c(1,2,3,9)
colors <- c("red", "blue", "darkgreen", "gold")

plot(x,
     dgamma(x=x,shape=k[1],scale=theta[1]),
     lwd=2, 
     col=colors[1],
     type="l",
     xlab="x value",
     ylab="Density",
     main="Comprassion Of Gamma Distributions")

for (i in 2:length(mean))
{
  lines(x, dgamma(x=x,shape=k[i],scale=theta[i]),type="l",
        lwd=2, col=colors[i])
}
labels=c()
for (i in 1:length(mean))
{
  labels <- c(labels,paste(sep="","theta - ",as.character(theta[i]),
                           ", k - ",as.character(k[i])))
}
legend("topright", inset=.05, title="Distributions",
       legend=labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
