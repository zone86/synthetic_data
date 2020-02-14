

######
# computation of effective correlations with confidence intervals, for different scales


setwd(paste0(Sys.getenv('CS_HOME'),'/FinancialNetwork/SyntheticAsset/Models'))
source('functions.R')

#month="201506"
#s1=read.csv(paste0('data/filtered/EURGBP_',month,'_600.csv'),header=FALSE)[,2];x1=log(s1/s1[1])
#s2=read.csv(paste0('data/filtered/EURUSD_',month,'_600.csv'),header=FALSE)[,2];x2=log(s2/s2[1])
x1=lpDataFromDir('data/filtered','EURUSD')
x2=lpDataFromDir('data/filtered','EURGBP')

omega0 = 288
filtering=3

plot(diff(gaussianFilter(x1,omega0)),type='l');points(diff(gaussianFilter(x2,omega0)),type='l',col='red')
cor(diff(gaussianFilter(x1,omega0)),diff( gaussianFilter(x2,omega0)))
cor(diff(x1[(2*omega0+1):(length(x1)-2*omega0)]-gaussianFilter(x1,omega0)),diff(x2[(2*omega0+1):(length(x2)-2*omega0)]-gaussianFilter(x2,omega0)))

s=synthAssets(x1,x2,correlations[3],omega0*3/filtering)
cor(diff(s$xs1),diff(s$xs2))
plot(gaussianFilter(x1,omega0),type='l');
points(s$xs1,type='l',col='red')
points(x1[(2*omega0+1):(length(x1)-2*omega0)],type='l',col='purple')



#########################

plot(x1,type='l')
sd(diff(x1))
sd(diff(x2))
cor(diff(x1),diff(x2))
cov(diff(x1),diff(x2))
x10 = gaussianFilter(x1,omega0);#x10=x10-mean(diff(x10));
x1r = x1[(2*omega0+1):(length(x1)-2*omega0)]-mean(diff(x10));
x20 = gaussianFilter(x2,omega0);#x20=x20-mean(diff(x20));
x2r = x2[(2*omega0+1):(length(x1)-2*omega0)]-mean(diff(x20));
plot(x1r,type='l');points(x10,type='l',col='red')
plot(x1r - x10,type='l')

(sd(diff(x10))/sd(diff(x1)))^2

y1 = x1r - x10
y2 = x2r - x20
rhoY=cor(diff(y1),diff(y2))
rhoT=cor(diff(x10),diff(x20))
rho=cor(diff(x1),diff(x2))
e1 = (sd(diff(x10))/sd(diff(x1)));e2=(sd(diff(x20))/sd(diff(x2)))
(e1*e2*rhoT + rhoY)*(1 - 0.5*(e1^2)+ e2^2)
(e1*e2*rhoT + rhoY)/sqrt((1+e1^2)*(1+e2^2))
# does not work, since Y \indep T does not hold directly on real data.


abs(mean(diff(x10))/max(diff(x10)))
abs(mean(diff(y1))/max(diff(y1)))

# retest synth asset
x1f=gaussianFilter(x = x1,12);x2f=gaussianFilter(x = x2,12)
correlations=seq(from=-0.8,to=0.8,by=0.05)
rhoeff=c()
for(rho in correlations){
  s=synthAssets(x1=x1f,x2=x2f,rho =rho,omega0 = 144)
  rhoeff=append(rhoeff,cor.test(diff(s$xs1),diff(s$xs2))$estimate)
}
plot(correlations,correlations,type='l');points(correlations,rhoeff,pch="+",col='red')

#######################

omega0 = 144
filterings = c(1,3,6,12,24)
correlations = seq(from=-0.95,to=0.95,by=0.05)
bootstrap_num = 1

effcorrs = c();effcorrsmin=c();effcorrsmax=c()
filt=c();expcorr=c();thcorr=c();e1mean=c();e2mean=c()

# first compute vol/cor at omega0 - estimated on all ts, a bit different than after second filtering,
# but should be negletible
x10=gaussianFilter(x1,omega0);x20=gaussianFilter(x2,omega0)
sigma10 = sd(diff(x10));sigma20 = sd(diff(x20))
rho0 = cor(diff(x10),diff(x20))

for(f in 1:length(filterings)){
  filtering=filterings[f]
  
  show(paste0('filtering = ',filtering))
  
  #x1 = sample(gaussianFilter(x1,filtering),filtering);x2 = sample(gaussianFilter(x2,filtering),filtering)
  # Do not sample - ok here as long as no pred. model
  xf1=gaussianFilter(x = x1,filtering);xf2=gaussianFilter(x = x2,filtering)
  
  synthlength = length(xf1)-(4*omega0)
  synth1 = matrix(0,bootstrap_num*synthlength,length(correlations))
  synth2 = matrix(0,bootstrap_num*synthlength,length(correlations))
  for(i in 1:length(correlations)){
    show(paste0('  rho = ',correlations[i]))
    for(k in 0:(bootstrap_num-1)){
       s=synthAssets(xf1,xf2,correlations[i],omega0);
       synth1[k*synthlength+(1:synthlength),i]=s$xs1
       synth2[k*synthlength+(1:synthlength),i]=s$xs2
    }
  }
  for(j in 1:length(correlations)){
    xs1 = synth1[,j];xs2 = synth2[,j]
    ctest = cor.test(diff(xs1),diff(xs2));
    sigmaY1 = sd(diff(xs1-x10[(2*filtering+1):(length(x10)-2*filtering)]));sigmaY2 = sd(diff(xs2-x20[(2*filtering+1):(length(x10)-2*filtering)]));e1=sigma10/sigmaY1;e2=sigma20/sigmaY2
    rhoY = correlations[j]
    #show(cor(diff(synth1[,j]),diff(synth2[,j])))
    effcorrs=append(effcorrs,ctest$estimate);
    effcorrsmin=append(effcorrsmin,ctest$conf.int[1])
    effcorrsmax=append(effcorrsmax,ctest$conf.int[2])
    filt=append(filt,filtering);
    expcorr=append(expcorr,correlations[j])
    thcorr=append(thcorr,(e1*e2*rho0 + rhoY)/sqrt((1+e1^2)*(1+e2^2)))
    e1mean=append(e1mean,e1);e2mean=append(e2mean,e2);
  } 
}

reseffcorrs = data.frame(
  rho=expcorr,
  rho_eff=effcorrs,
  rhomin=effcorrsmin,
  rhomax=effcorrsmax,
  rhoth=thcorr,
  filtering=filt*600/3600
)

save(reseffcorrs,file='res/effcorrs.RData')
library(RColorBrewer)

discrfilt = ifelse(filt==1,"10min",ifelse(filt==3,"30min",ifelse(filt==6,"1h",ifelse(filt==12,"2h","4h"))))
reseffcorrs$discrfilt = factor(discrfilt,levels = c("10min","30min","1h","2h","4h"))

g=ggplot(reseffcorrs,aes(x=rho,y=rho_eff,colour=discrfilt))
g+geom_point(shape=19)+geom_line(linetype=2,size=0.4)+
  geom_errorbar(aes(ymin=rhomin,ymax=rhomax),width=0.03)+
  #scale_color_gradientn(colours=rainbow(length(filterings)),name=expression(omega[1]*" (h)"))+
  #scale_color_gradientn(colours=brewer.pal(length(filterings),"Paired"),name=expression(omega[1]*" (h)"))+
  scale_color_discrete(name=expression(omega[1]))+
  geom_line(aes(x=rho,y=rhoth,group=filtering,colour=discrfilt),size=1)+
  geom_abline(slope=1,intercept = 0,colour="black",linetype=2,size=1)+
  geom_vline(xintercept=2*rho0*mean(e1mean)*mean(e2mean)/(mean(e1mean)^2+mean(e2mean)^2),colour='red')+
  xlab(expression("Targeted correlation "*rho))+ylab(expression("Effective correlation "*rho[e]))+
  #scale_x_log10()+scale_y_log10()+ # must add 1 - makes no sense anyway
  stdtheme
ggsave(file='../Results/Cors/effectiveCorrs.png',width=24,height=20,units='cm')






