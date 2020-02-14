

# functions


##
#  file names assumed as asset+"_"+sortable_date+end
lpDataFromDir <-function(dir,asset){
  files <- list.files(dir)
  assetfiles = c()
  for(file in files){if(substr(file,1,nchar(asset))==asset){assetfiles=append(assetfiles,substr(file,nchar(asset)+2,nchar(file)))}}
  s=c()
  for(file_suffix in sort(assetfiles)){
    s1=read.csv(paste0(dir,'/',asset,'_',file_suffix),header=FALSE)[,2];
    s = s - s[length(s)]+s1[1]
    s=append(s,s1)
  }
  x = log(s/s[1])
  return(x)
}



gaussianFilter <- function(x,sigma){
  k = exp(- (((-2*sigma):(2*sigma))^2) / (sigma/2)^2)
  k = k/sum(k)
  return(convolve(x,k,type="filter"))
}

##
# sampling by providing sampling step
sample <- function(x,step){
  s = seq(from=1,to=length(x),by=max(floor(step),1))
  return(x[s])
}

brownian <- function(sigma,steps){
  cumsum(rnorm(steps,0,sigma))
}


laggedCorrs <- function(x1,x2,taumin,taumax,taustep){
  taus = seq(from=taumin,to=taumax,by=taustep)
  corrs=c();corrsmin=c();corrsmax=c()
  for(tau in taus){
    # determine common support bounds
    if(tau <= 0){
      x2e = x2[(1-tau):length(x2)];x1e= x1[1:(length(x1)+tau)]
    }else{
      x2e = x2[1:(length(x2)-tau)];x1e= x1[(1+tau):length(x1)]
    }
    t = cor.test(x1e,x2e)
    corrs=append(corrs,t$estimate);corrsmin=append(corrsmin,t$conf.int[1]);corrsmax=append(corrsmax,t$conf.int[2])
  }
  return(data.frame(taus,corrs,corrsmin,corrsmax))
}


##
#
synthAssets<- function(x1,x2,rho,omega0){
  # base components
  x10 = gaussianFilter(x1,omega0);
  x20 = gaussianFilter(x2,omega0)
  #x10 = x10 - mean(diff(x10));x20 = x20 - mean(diff(x20))
  
  # reduced support for assets
  x1r = x1[(2*omega0+1):(length(x1)-2*omega0)];
  x2r = x2[(2*omega0+1):(length(x1)-2*omega0)]
  
  # estimate sd for black-scholes
  sigma1 = sd(diff(x1r));
  sigma2 = sd(diff(x2r))
  
  if(sigma1<sigma2){
    xs1 = brownian(sigma1,length(x1))
    xs1par = brownian(sigma1,length(x1))
    xs2 = rho*xs1 + sqrt(1 - ((sigma1/sigma2)*rho)^2)*xs1par
  }else{
    xs2 = brownian(sigma2,length(x2))
    xs2par = brownian(sigma2,length(x2))
    xs1 = rho*xs2 + sqrt(1 - ((sigma2/sigma1)*rho)^2)*xs2par
  }
  
  #return(data.frame(xs1=xs1[(2*omega0+1):(length(xs1)-2*omega0)],xs2=xs2[(2*omega0+1):(length(xs1)-2*omega0)]))
  return(data.frame(xs1=x10+(xs1[(2*omega0+1):(length(xs1)-2*omega0)]-gaussianFilter(xs1,omega0)),xs2=x20+(xs2[(2*omega0+1):(length(xs2)-2*omega0)]-gaussianFilter(xs2,omega0))))
}









