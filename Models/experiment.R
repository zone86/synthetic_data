

#########
## Experiment
#########

resdir = 'res/'
assets = c("EURUSD","EURGBP")

filterings = c(6,9,12)
correlations = seq(from=-0.95,to=0.95,by=0.05)
omega0 = 144

source('functions.R')

library(doParallel)
cl <- makeCluster(40)
registerDoParallel(cl)


for(filtering in filterings){
    # load data
    x1=lpDataFromDir('data/filtered','EURUSD')
    x2=lpDataFromDir('data/filtered','EURGBP')
  
    # refilter
    xf1 = gaussianFilter(x1,filtering);
    xf2 = gaussianFilter(x2,filtering);
    # generate synthetic data for each level of correlation
    synthlength = length(xf1)-(4*omega0)
    synth = matrix(0,2*synthlength,length(correlations))
    for(i in 1:length(correlations)){
      s=synthAssets(xf1,xf2,correlations[i],omega0);
      synth[,i]=c(s[,1],s[,2])
    }
    
    rhoeff=c()
    for(j in 1:length(correlations)){
      x1s=synth[1:synthlength,j];
      x2s=synth[(synthlength+1):(2*synthlength),j];
      rhoeff = append(rhoeff,cor(diff(x1s),diff(x2s)))
    }
    
    # compute perfs in //
    #   estimated computing time : for t \in [1:100] at filtering=6 : tau = 160s
    #      length(x1s)=6000 -> tautot = tau*60 ~ 3h !
    #   filtering = 12 -> tautot ~ 1,5h
    #
    res <- foreach(j=1:ncol(synth)) %dopar% {
      source('functions.R');source('models.R')
      x1s=sample(synth[1:synthlength,j],filtering/2);x2s=sample(synth[(synthlength+1):(2*synthlength),j],filtering/2);
      #t=system.time(predictionMSE(x1s[1:200],x2s[1:200],288*2/filtering,2));t
      m=predictionMSE(x1s,x2s,288*2/filtering,2)
      error=(m$expected-m$pred)^2
      res=c(error[,1],error[,2])
      res
    }
    
    # get results into data frame
    vals_mat = matrix(0,length(res),length(res[[1]]))
    for(a in 1:length(res)){vals_mat[a,]=res[[a]]}
    v = data.frame(vals_mat);
    
    write.table(v,file=paste0(resdir,"all_",filtering,".csv"),sep=";",row.names=FALSE,col.names=FALSE)
    write.table(rhoeff,file=paste0(resdir,"rhoeff_",filtering,".csv"),sep=";",row.names=FALSE,col.names=FALSE)
    
}


stopCluster(cl)
