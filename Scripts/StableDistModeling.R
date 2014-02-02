require(stabledist)

model_function <- function(x, params) {
  dstable(x, 2, params[1], params[2], params[3])
}

#p=findParams(sample=logRet, initParams=c(0.5,0.5,0.5,0.5),
#           densityOfXAndvectorOfParams=model_function, 
#           method="L-BFGS-B", lower = c(0,-1,0,-Inf) + 0.01, upper = c(2,1,Inf,Inf) - 0.01)

p=findParams(sample=logRet, initParams=c(0.5,0.5,0.5),
           densityOfXAndvectorOfParams=model_function, 
           method="L-BFGS-B", lower = c(-1,0,-Inf) + 0.0001, upper = c(1,Inf,Inf) - 0.0001)



saveParamsToFile(p$par,"params_stable_alpha_2")


params <-c(1.800588,-0.128442,0.01,0.0002265085);
plot(logRet, type ='l')
x=rstable(length(logRet),alpha=params[1], 
          beta = params[2], gamma=params[3], delta = params[4])

lines(x, col = 'red')
