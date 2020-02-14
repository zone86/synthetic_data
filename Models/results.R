
# Results of experiments

setwd(paste0(Sys.getenv('CS_HOME'),'/FinancialNetwork/SyntheticAsset/Models'))


### Tests

#res <- read.csv('res/201506_3.csv',sep=";",header=FALSE)
#res=t(res)

# rough plot
#plot(colMeans(res))
#colMeans(res)
#apply(res,2,sd)

# 95% CI on means
#apply(res,2,function(c){2*1.96*sd(c)/sqrt(length(c))})


##############
##
resdir = 'res/'
filterings = c(6,9,12)
correlations = seq(from=-0.95,to=0.95,by=0.05)
assets = c('EURUSD','EURGBP')

x1=lpDataFromDir('data/filtered','EURUSD')
x2=lpDataFromDir('data/filtered','EURGBP')

rho=c();rhoeff=c();filt=c();asset=c();perf=c();perfmin=c();perfmax=c()

for(filtering in filterings){
  res <- t(read.csv(paste0(resdir,'all_',filtering,'.csv'),sep=";",header=FALSE))
  rescor <-c(read.csv(paste0(resdir,'rhoeff_',filtering,'.csv'),sep=";",header=FALSE))$V1
  # noramlize by sigma2 -> need to recompute it
  
  xf1 = gaussianFilter(x1,filtering);
  xf2 = gaussianFilter(x2,filtering);
  sigma1 = sd(diff(xf1));sigma2 = sd(diff(xf2));sigmas=list();sigmas[[assets[1]]]=sigma1;sigmas[[assets[2]]]=sigma2;
  for(j in 1:length(correlations)){
    d = data.frame(res[1:floor(nrow(res)/2),j],res[(floor(nrow(res)/2)+1):nrow(res),j]);colnames(d)<-assets
    for(a in assets){
      rho=append(rho,correlations[j]);rhoeff=append(rhoeff,rescor[j]);asset=append(asset,a);filt=append(filt,filtering)
      p=d[[a]]/sigmas[[a]]^2
      perf=append(perf,mean(p))
      cisize = 1.96*sd(p)/sqrt(length(p))
      # CI on mean
      perfmin=append(perfmin,mean(p)-cisize);perfmax=append(perfmax,mean(p)+cisize)
    }
  }  
}

d=data.frame(rho,rhoeff,filt,asset,perf,perfmin,perfmax)
#filtering=12

for(filtering in filterings){
  g=ggplot(d[d$filt==filtering,],aes(x=rhoeff,y=perf,colour=asset))
  g+geom_point()+geom_errorbar(aes(ymin=perfmin,ymax=perfmax))+ggtitle(bquote(omega[1]*"="*.(filtering/6)*"h"))+stat_smooth()+
    xlab(expression('Effective correlation '*rho[e]))+ylab(expression('Performance '*pi))+
    scale_color_discrete(name="Asset")+stdtheme
  ggsave(file=paste0('../Results/Prediction/predictionPerformance-rhoeff_filt',filtering,'.png'),width=22,height=16,units='cm')
}

######
# lag corr between fundamental components -> cause of different behavior ?


filterings=c(6,9,12,144)
tau=c();filt=c();rhol=c();rholmin=c();rholmax=c()

for(filtering in filterings){
  l = laggedCorrs(x1 = diff(gaussianFilter(x1,filtering)),x2=diff(gaussianFilter(x2,filtering)),taumin = -48,taumax=48,taustep=1)
  filt=append(filt,rep(filtering,nrow(l)))
  tau=append(tau,l$taus);rhol=append(rhol,l$corrs);rholmin=append(rholmin,l$corrsmin);rholmax=append(rholmax,l$corrsmax)
}

d=data.frame(tau=tau/6,omega1=factor(paste0(as.character(filt/6),"h"),levels=c("1h","1.5h","2h","24h")),rhol,rholmin,rholmax)
#cols=c("#00BFC4","#F8766D","#C77CFF","#7CAE00")
#cols=rainbow(24)[c(1,3,8,20)]
#xlim(-48,48)

g=ggplot(d,aes(x=tau,y=rhol,colour=omega1))
g+geom_point(shape=19)+
  geom_errorbar(aes(ymin=rholmin,ymax=rholmax))+
  geom_line(aes(x=tau,y=rhol,colour=omega1,group=omega1))+
  #scale_color_manual(values=cols,breaks=c("1","1.5","2","24"),name="omega1 (h)")+
  scale_color_discrete(name=expression(omega[1]))+
  geom_vline(xintercept=c(-3,-5.95,3,5.95)/6,col=rep(cols[1],4),linetype=rep(2,4))+
  geom_vline(xintercept=c(-4.5,-9,4.5,9)/6,col=rep(cols[2],4),linetype=rep(2,4))+
  geom_vline(xintercept=c(-6.05,-12,6.05,12)/6,col=rep(cols[3],4),linetype=rep(2,4))+
  coord_cartesian(xlim = c(-8, 8))+xlab(expression('Temporal lag '*tau*' (h)'))+
  ylab(expression('Lagged correlation '*rho[tau]))+stdtheme
ggsave(file='../Results/Prediction/laggedCorrelations.png',width=22,height=16,units='cm')
  
#unique(ggplot_build(g+geom_point(shape=19))$data[[1]]$colour)


