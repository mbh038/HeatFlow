# Hukseflux heat Flux Meter
rhfm=0.00625 # W/m2


#let Vin=2*(T1-T2) - assuming U value of 2)
#lfit<-lm(formula=Vin~Vout,data=a)

# amplifier functions from assuming U=2.0
ampch1<-function(vout){
    (0.003081 * vout +1.152144) #cal 24-02-17 lynda & mike
}

ampch2<-function(vout){
  (0.0003245 * vout +1.066) #cal 24-02-17 lynda & mike
}

# heat flux function assuming Bidulph model
Qt<-function(R1,R2,Tm_init,C,Tint,Text){
  Tm[1]=Tm_init
  for (i in 1:(length(t)-1)){
    dt=60*(t[i+1]-t[i])
    Q[i]=(Tint[i]-Tm[i])/R1
    Tm[i+1]=((Tint[i+1]/R1)+(Text[i+1]/R2)+C*Tm[i]/dt)/(1/R1 + 1/R2 + C/dt)
  }
  Q
}

#heat flux model assuming zero heat capacity
Qins<-function(R3,Tint,Text){
  Q=(Tint-Text)/R3
  Q
}

### Window trial   
df<-read.table("trial002.csv",sep=",",stringsAsFactors=FALSE,header=TRUE)
id<-seq(1,nrow(df))
df<-cbind(id,df)
names(df)<-c("id","date","time","T1","T2","T3","T4","hp1","hp2")

ymax=max(max(df$T1,na.rm=TRUE),max(df$T3,na.rm=TRUE))
ymin=min(min(df$T2,na.rm=TRUE),min(df$T4,na.rm=TRUE))

library(rafalib)
mypar(2,1)
plot(df$T1,type="l",ylim=c(ymin,ymax)) # inside
lines(df$T2,col=2) # outside
lines(df$T3,col=3) # inside
lines(df$T4,col=4) # outside
legend("topright", c("Ch1 Tin", "Ch1 Tout","Ch2 Tin", "Ch2 Tout"), pch="-",col=c(1,2,3,4),cex=0.7)

ymax=max(df$T1-df$T2,df$T3-df$T4,na.rm=TRUE)
ymin=min(df$T1-df$T2,df$T3-df$T4,na.rm=TRUE)

plot(df$T1-df$T2,type="l",ylim=c(ymin,ymax)) # Ch1
lines(df$T3-df$T4,col=2) # Ch2
legend("topright", c("Ch1", "Ch2"), pch="-",col=c(1,2))

#find amplifier calibration constants, assuming U=2
df$hp1Vin=2*(df$T1-df$T2)*0.06
hp1fit<-lm(formula=df$hp1Vin~df$hp1,data=df)
df$hp2Vin=2*(df$T3-df$T4)*0.06
hp2fit<-lm(formula=df$hp2Vin~df$hp2,data=df)



mypar(2,1)
# heat plate plots
df$hp1preamp1<-ampch1(df$hp1)
df$hp1flux1<-df$hp1preamp1/0.060
df$hp2preamp2<-ampch2(df$hp2)
df$hp2flux2<-df$hp2preamp2/0.060


#index=1:(nrow(df))
index=4500:5000

plot(df$id[index],df$hp1[index],type="l",main="Vout from amp Channel 1")
lines(df$id[index],df$hp2[index],col="red",type="l")
legend("topright", c("Ch1", "Ch2"), pch="-",col=c("black", "red"))

plot(df$id[index],df$hp1flux1[index],type="l",main="Heat flux hp 1")
lines(df$id[index],df$hp2flux2[index],col="red",type="l")
legend("topright", c("Ch1", "Ch2"), pch="-", col=c("black", "red"))

Qexp1<-df$hp1flux1[index]
Qexp2<-df$hp2flux2[index]
#Qexp<- -Qexp[-1]+15
Tint1<-df$T1[index]
Text1<-df$T2[index]
Tint2<-df$T3[index]
Text2<-df$T4[index]
t=df$id[index]

Tm=numeric()
Q1=numeric()
Q2=numeric()

mypar(1,1)
R1=0.2
R2=0.2
C=10000
Tm_init=13
tau=60

R3=0.4

Q1<-Qt(R1,R2,Tm_init,C,Tint1,Text1)
Qi1<-Qins(R3,Tint1,Text1)
Q1<-append(Q1,tail(Q1,n=1))
Qi1<-append(Qi1,tail(Qi1,n=1))

Q2<-Qt(R1,R2,Tm_init,C,Tint2,Text2)
Qi2<-Qins(R3,Tint2,Text2)
Q2<-append(Q2,tail(Q2,n=1))
Qi2<-append(Qi2,tail(Qi2,n=1))

mypar(2,1)
plot(Qexp1,type="l",ylim=c(min(min(Q1),min(Qexp1)),max(max(Q1),max(Qexp1))),xlab="Time (min)",ylab="Heat flux Q1 (W/m^2)",col="red")
lines(Q1,col="blue")
lines(Qi1,col="cyan")
legend("topright", c("Measured", "Predicted"), pch="o", col=c("red", "blue"))

plot(Qexp2,type="l",ylim=c(min(min(Q2),min(Qexp2)),max(max(Q2),max(Qexp2))),xlab="Time (min)",ylab="Heat flux Q2 (W/m^2)",col="red")
lines(Q2,col="blue")
lines(Qi2,col="cyan")
legend("topright", c("Measured", "Predicted"), pch="o", col=c("red", "blue"))

LL <- function(R1,R2,Tm_init,C, mu, sigma) {
  Rem = Qexp-Qt(R1,R2,Tm_init,C,Tint,Text)
  #
  Rem = suppressWarnings(dnorm(Rem, mu, sigma, log = TRUE))
  #
  -sum(Rem)
}

LLins <- function(R3, mu, sigma,Tint,Text) {
  R = Qexp-Qins(R3,Tint,Text)
  #
  R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
  #
  -sum(R)
}


library(stats4)
fit4p1<-mle(LL, 
           start = list(R1=0.25,R2=0.25,Tm_init=15,C=10,sigma=4),
           fixed=list(mu=0),
           nobs = length(Q),
           lower = c(.01,.01,10,.1,0.1),
           upper = c(2,2.,20,1000,5),
           method= "L-BFGS-B"
)

fitins1<-mle(LLins,
           start = list(R3=0.5,sigma=4),
           fixed=list(mu=0),
           nobs = length(Q),
           lower = c(.01,0.1),
           upper = c(2,5),
           method= "L-BFGS-B"
)

Q1<-Qt1(coef(fit4p1)[1],coef(fit4p1)[2],coef(fit4p1)[3],coef(fit4p1)[4],Tint1,Text1)
lines(Q1,col="grey")

Qi1<-Qins(coef(fitins)[1],Tint1,Text1)
lines(Qi1,col="green")

print(summary(fit4p1))
print(summary(fitins1))

U=1/(coef(fit4p1)[1]+coef(fit4p1)[2]-rhfm)
print(U)

U2=1/(coef(fitins1)[1]-rhfm)
print(U2)
