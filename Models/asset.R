
# synthetic correlation between assets

setwd(paste0(Sys.getenv('CS_HOME'),'/FinancialNetwork/SyntheticAsset/Models'))

library(ggplot2)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))
source('functions.R')

# tests asset viz
eurusd <- read.csv('data/raw/EURUSD_201511.csv',header=FALSE)
eurgbp <- read.csv('data/raw/EURGBP_201511.csv',header=FALSE)

# package smoother -- NO, shitty complexity
library(smoother)


# should proceed to temporal adjustement here : different lengths !
# -> dirty part of data normalization and completion

formatDates <- function(d){as.numeric(format(as.POSIXct(substr(d,1,15),format="%Y%m%d %H%M%S",tz="UTC"),format="%s"))}
sapply(eurusd[,1],formatDates)

# midmarket signal
s1 = (eurusd[,2]+eurusd[,3])/2
x1 = log(s1/s1[1])
dx1 = diff(x1)
s2 = (eurgbp[,2]+eurgbp[,3])/2
x2 = log(s2/s2[1])
dx2 = diff(x2)


# plots 
x=x1

#plot(x[1801:12600],type='l',main='EUR/USD, 1st November 2015',ylab='log(s/s0)',xlab='time (s)')
#points(convolve(x[1801:12600],k,type="filter"),col='red',type='l')
#filtering test
#f1=smth.gaussian(x,window=1800)
#f2=smth.gaussian(x,window=450)

f1=smth.gaussian(x[1:14400],window=1800,tails = T)
f2=smth.gaussian(x[1:14400],window=300,tails = T)

#points(f1[1801:12600],col='red',type='l')
#points(f2[1801:12600],col='purple',type='l')

l = length(1801:12600)
d = data.frame(t=rep(1:l,3),
               x=c(x[1801:12600],f1[1801:12600],f2[1801:12600]),
               filtering = factor(c(rep("raw",l),rep("30min",l),rep("5min",l)),levels=c("raw","5min","30min"))
               )
# redo as ggplot
#pal<-function(n){return(c("black","red","green"))}

g = ggplot(d[d$filtering!="raw",],aes(x=t,y=x,group=filtering,color=filtering))
g+geom_line(data=d[d$filtering=="raw",],size=0.3,color='black')+
  geom_line(size=1)+#ggtitle('EUR/USD, 1st November 2015')+
  xlab('Time (s)')+ylab(expression("log"*(s/s[0])))+
  scale_color_discrete(name=expression(omega))+stdtheme
ggsave(file='../Results/Stylized/filteringExample_20151101_5-30min.png',width=24,height=20,units='cm')


### same plot with 10min more readable ?

f2=smth.gaussian(x[1:14400],window=600,tails = T)
d = data.frame(t=rep(1:l,3),x=c(x[1801:12600],f1[1801:12600],f2[1801:12600]),filtering = factor(c(rep("raw",l),rep("30min",l),rep("10min",l)),levels=c("raw","10min","30min")))


g = ggplot(d[d$filtering!="raw",],aes(x=t,y=x,group=filtering,color=filtering))
g+geom_line(data=d[d$filtering=="raw",],size=0.3,color='black')+
  geom_line(size=1.5)+
  #ggtitle('EUR/USD, 1st November 2015')+
  xlab('Time (s)')+ylab(expression("log"*(s/s[0])))+xlim(0,l)+
  scale_color_discrete(name=expression(omega))+stdtheme
ggsave(file='../Results/Stylized/filteringExample_20151101_10-30min.png',width=28,height=18,units='cm')



###########

# stylized facts ? -> check fat tail
hist(dx1,breaks=1000)


### Filtering : smoothing function is totally shitty in complexity ; must use convolve



# corr synth data -> correlation between log-returns at a given scale.

# natural corrs ?
# check lag
df1 = gaussianFilter(dx1,1800)
df2 = gaussianFilter(dx2,1800)
corrs=c()
for(k in (0:240)*60+1){
  corrs=append(corrs,cor(f1[7201:(25*3600)],f2[(1:(23*3600))+k]))
}


# noise at 30min
plot(df1,type='l')
points(df2,type='l',col='red')
plot(dx1,type='l')

# Q : how correlate at given scale ?
#   -> estimate sigma of Black-Scholes after sampling.
#  pb : correlates returns or signal ? if returns, integration errors ?


# check if function works well
hist(diff(brownian(0.001,100000)),breaks=1000)
hist(df1[df1<quantile(df1,0.999)&df1>quantile(df1,0.001)],breaks=1000)

# very rough estimation : sd(f1) is sigma for the blackshcoles cumsum(f1)
sigma = sd(df1)
f1 = cumsum(df1)
# f1 can be sampled at 10min as filtered at 30min, should be no aliasing
f1 = f1[seq(from=1,to=length(f1),by=600)]



#quick test
plot(diff(f1),type='l')
corrs=c()
for(k in 1:1000){
  if(k%%100==0){show(k)}
  synth = synthAsset(f1,0.7)
  corrs = append(corrs,cor(diff(synth),diff(f1)))
  #points(diff(synth),type='l',col=k)
}

hist(corrs,breaks=50)



# test synth asset
omega0=72
s=synthAssets(x1,x2,0.01,omega0)
plot(x1[(2*omega0+1):(length(x1)-2*omega0)],type='l');points(s[,1],type='l',col='red')

# TODO : plot of examples of synth data




