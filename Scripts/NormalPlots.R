

x <- seq(-10, 10, length=100)


mean <- c(0,1,0,0)
sd <- c(1,2,3,4)
colors <- c("red", "blue", "darkgreen", "gold")

plot(x,
     dnorm(x=x,mean=mean[1],sd=sd[1]),
     lwd=2, 
     col=colors[1],
     type="l",
     xlab="x value",
     ylab="Density",
     main="Comprassion Of Normal Distribution")

for (i in 2:length(mean))
{
  lines(x, dnorm(x=x,mean=mean[i],sd=sd[i]), lwd=2, col=colors[i])
}
labels=c()
for (i in 1:length(mean))
{
  labels <- c(labels,paste("mean - ",as.character(mean[i]),
                           "sd - ",as.character(sd[i])))
}
legend("topright", inset=.05, title="Distributions",
       legend=labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
