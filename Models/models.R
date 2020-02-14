
######
# predictive models whose power has to be tested on synthetic data
######



# basic auto-regressive model
#  -> use package MTS

library(MTS)

#d = data.frame(eurusd=diff(log(test$EURUSD$mid/test$EURUSD$mid[1])),eurgbp=diff(log(test$EURGBP$mid/test$EURGBP$mid[1])))
#s1=read.csv('data/filtered/EURGBP_201506_600.csv',header=FALSE)[,2];x1=log(s1/s1[1])
#s2=read.csv('data/filtered/EURUSD_201506_600.csv',header=FALSE)[,2];x2=log(s2/s2[1])

# Tests
#m = VARMA(d[1:100,],p=10,q=0)
# forecasts then given by VARMApred
#VARMApred(m,h=1)



predictionMSE <- function(x1,x2,Tw,p){
  dx1 = diff(x1)
  dx2 = diff(x2)
  d = data.frame(dx1,dx2)
  
  pred = matrix(NA,length(dx1)-Tw,2)
  
  for(t in 1:nrow(pred)){
    m = tryCatch(VARMA(d[t:(t+Tw-1),],p=p,q=0))
    pred[t,] = tryCatch(VARMApred(m)$pred)
  }
  
  # mse can be directly computed on returns as prediction done step by step
  mse=c(0,0)
  
  expected = d[(Tw+1):length(dx1),]
  
  res=list()
  
  res[["expected"]] = expected
  res[["pred"]] = pred
  
  return(res)
}


# test
#x1=sample(gaussianFilter(log(test$EURUSD$mid/test$EURUSD$mid[1]),180),180)
#x2=sample(gaussianFilter(log(test$EURGBP$mid/test$EURGBP$mid[1]),180),180)

#mse = predictionMSE(x1,x2,144,2)

#plot(abs(mse$expected[,1]-mse$pred[,1]),type='l',col='red');points(abs(mse$expected[,1]),type='l')






