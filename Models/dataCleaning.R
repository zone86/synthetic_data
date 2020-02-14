
#setwd(paste0(Sys.getenv('CS_HOME'),'/FinancialNetwork/SyntheticAsset/Models'))

# clean and filter data for each month

formatDates <- function(d){as.numeric(format(as.POSIXct(substr(d,1,15),format="%Y%m%d %H%M%S",tz="UTC"),format="%s"))}

gaussianFilter <- function(x,sigma){
  # kernel : cut at +- 2 sigma
  k = exp(- (((-2*sigma):(2*sigma))^2) / (sigma/2)^2)
  k = k/sum(k)
  return(convolve(x,k,type="filter"))
}


extractFromSupport<-function(reducedsupport,support,ts){
  # start is position of first of reduced support in support
  start = which(support==reducedsupport[1])[1]
  pend =  which(support==reducedsupport[length(reducedsupport)])
  end = pend[length(pend)]
  res=c(ts[start]) 
  suppind = 1
  
  show(paste0("start : ",start," ; end : ",end))

  for(k in start:(end-1)){
    # check if jump in reduced support compared to support
    # if no jump 
   #if(!is.numeric(reducedsupport[suppind + 1])){show(paste0("suppind : ",suppind," ; ",reducedsupport[suppind + 1]," ; k : ",k," ; ",support[k + 1]))}
   #if(!is.numeric(support[k + 1])){show(paste0("suppind : ",suppind," ; ",reducedsupport[suppind + 1]," ; k : ",k," ; ",support[k + 1]))}

   #show(paste0("suppind : ",suppind," ; ",reducedsupport[suppind + 1]," ; k : ",k," ; ",support[k + 1]))

   if(suppind<length(reducedsupport)){  
     if(reducedsupport[suppind + 1]==support[k+1]){
        # if current index different, at the end of a jump, must adjust res ts
        if(support[k]!=reducedsupport[suppind]){
          res = res + (ts[k+1] - res[length(res)])
        }
        res=append(res,ts[k+1])
        #show(paste0("supind : ",suppind," , k : ",k," , res = ",length(res))) #DEBUG
        suppind=suppind+1
      }
    }
  } 
  return(res)
}


cleanData <- function(dir,month,assets){
  show(paste0("Cleaning data : ",length(assets)," assets, month : ",month))
  # load the data
  show("-> loading data")
  d = list()
  for(asset in assets){d[[asset]]=read.csv(paste0(dir,"/",asset,'_',month,'.csv'),header=TRUE)}
  # convert times
  show("-> converting timestamps")
  for(asset in assets){
       d[[asset]][,1]=sapply(d[[asset]][,1],formatDates)
       d[[asset]] = d[[asset]][!is.na(d[[asset]][,1]),]
  }
  # compute mids
  show("-> computing mids")
  for(asset in assets){d[[asset]]$mids=(d[[asset]][,2]+d[[asset]][,3])/2}
  # adjust TS
  show("-> computing common support")
  support = d[[assets[1]]][,1]
  for(i in 2:length(assets)){support=intersect(support,d[[assets[i]]][,1])}
  support = unique(support)
  support = support[!is.na(support)]
  # compute clean ts
  show("-> cleaning ts")
  res=list()
  for(asset in assets){res[[asset]]=data.frame(ts=support,mid=extractFromSupport(support,d[[asset]][,1],d[[asset]]$mids))}
  
  return(res)
}

##
# filter clean data
filterData<-function(data,sigma){
  show("Filtering data...")
  res = list()
  for(k in 1:length(data)){
      x = data[[names(data)[k]]]$mid
      filtersupport = ((2*sigma)+1):(length(x)-2*sigma)
      f = gaussianFilter(x,sigma)
      fts = data[[names(data)[k]]]$ts[filtersupport]
      # now subsample in sigma / 3
      sample = seq(from=1,to=length(f),by=floor(sigma/3))
      res[[names(data)[k]]] = data.frame(ts=fts[sample],fmid=f[sample])
  }
  return(res)
}


# test
#test <- cleanData('data/test',"201506",c("EURGBP","EURUSD"))

assets = c("EURUSD","EURGBP")
months = paste0("2015",c("08","09","10","11","07"))
dir = "data/raw"

filtering = 600

for(month in months){
   data = cleanData(dir,month,assets)
   # export raw data
   for(asset in assets){
     write.table(data[[asset]],file=paste0("data/clean/",asset,"_",month,".csv"),row.names=FALSE,col.names=FALSE,sep=",")
   }
   
   # filter, sample and export
   fdata = filterData(data,filtering)
   # export filtered
   for(asset in assets){
     write.table(fdata[[asset]],file=paste0("data/filtered/",asset,"_",month,"_",filtering,".csv"),row.names=FALSE,col.names=FALSE,sep=",")
   }
   
}






