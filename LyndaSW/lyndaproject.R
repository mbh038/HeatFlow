# Hukseflux heat Flux Meter
rhfm=0.00625 # W/m2

#let Vin=2*(T1-T2) - assuming U value of 2)
#lfit<-lm(formula=Vin~Vout,data=a)

# amplifier functions from assuming U=2.0
ampch1<-function(vout){
    (0.003081 * vout +1.152144) #cal 24-02-17 lynda & mike
}

# heat flux function

Qt<-function(R1,R2,Tm_init,C){
  Tm[1]=Tm_init
  for (i in 1:(length(t)-1)){
    dt=60*(t[i+1]-t[i])
    Q[i]=(Tint[i]-Tm[i])/R1
    Tm[i+1]=((Tint[i+1]/R1)+(Text[i+1]/R2)+C*Tm[i]/dt)/(1/R1 + 1/R2 + C/dt)
  }
  Q
}

### Window trial   
df<-read.table("trial001.csv",sep=",",stringsAsFactors=FALSE,header=TRUE)
id<-seq(1,nrow(df))
df<-cbind(id,df)
names(df)<-c("id","date","time","T1","T2","T3","T4","hp1","hp2")

ymax=max(df$T1,na.rm=TRUE)
ymin=min(df$T2,na.rm=TRUE)

library(rafalib)
mypar(1,1)
plot(df$T1,type="l",ylim=c(ymin,ymax)) # inside
lines(df$T2,col=2) # outside
lines(df$T3,col=3) # inside
lines(df$T4,col=4) # outside

mypar(2,1)
# heat plate plots
df$hp1preamp1<-ampch1(df$hp1)
df$hp1flux1<-df$hp1preamp1/0.060


index=1:(nrow(df))

plot(df$id[index],df$hp1[index],type="l",main="Vout from amp Channel 1")
#lines(df$id[index],df$hp2[index],col="blue",type="l")

plot(df$id[index],df$hp1flux1[index],type="l",main="Heat flux hp 1")
#lines(df$id[index],df$hp2flux2[index],col="blue",type="l")

Qexp<-df$hp1flux1[index]
#Qexp<- -Qexp[-1]+15
Tint<-df$T1[index]
Text<-df$T2[index]
t=df$id[index]

Tm=numeric()
Q=numeric()

mypar(1,1)
R1=0.25
R2=0.25
C=10000
Tm_init=13
tau=60

Q<-Qt(R1,R2,Tm_init,C)
Q<-append(Q,tail(Q,n=1))

plot(Qexp,type="l",ylim=c(min(min(Q),min(Qexp)),max(max(Q),max(Qexp))),xlab="Time (min)",ylab="Heat flux Q (W/m^2)",col="red")
lines(Q,col="blue")
legend("topright", c("Measured", "Predicted"), pch="o", col=c("red", "blue"))

LL <- function(R1,R2,Tm_init,C, mu, sigma) {
  R = Qexp-Qt(R1,R2,Tm_init,C)
  #
  R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
  #
  -sum(R)
}


library(stats4)
fit4p<-mle(LL, 
           start = list(R1=0.25,R2=0.25,Tm_init=15,C=10,sigma=2),
           fixed=list(mu=0),
           nobs = length(Q),
           lower = c(.01,.01,10,.1,0.1),
           upper = c(2,2.,20,1000,5),
           method= "L-BFGS-B"
)

Q<-Qt(coef(fit4p)[1],coef(fit4p)[2],coef(fit4p)[3],coef(fit4p)[4])
lines(Q,col="green")

print(summary(fit4p))

U=1/(coef(fit4p)[1]+coef(fit4p)[2]-rhfm)
print(U)
