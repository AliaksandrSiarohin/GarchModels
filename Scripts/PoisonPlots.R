

x <- seq(0, 15, by=1)


lambda <- c(0.1,1,5,10)

colors <- c("red", "blue", "darkgreen", "gold")
matrix=c()
for (i in 1:length(mean))
{
   matrix=c(matrix,dpois(x=x,lambda=lambda[i]))
}


barplot(matrix(matrix,nrow=4,byrow=TRUE), main="Comprassion of Poison distribution",
        ylab= "probability",beside=TRUE, col=colors,names.arg=x)

labels=c()
for (i in 1:length(mean))
{
  labels <- c(labels,paste("lambda - ",as.character(lambda[i])))
}
legend("topright", inset=.05, title="Distributions",
       legend=labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
