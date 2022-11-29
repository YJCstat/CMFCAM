#' Estimating the parameters of the frailty model
#' @param T1 observed mediator event time (vector)
#' @param T2 observed terminal event time (vector)
#' @param d2 1 for terminal event occured 0 for censored (vector)
#' @param tol maximum tolerance of change during the iteration
#' @param step maximum  number of the iteration
#' @param int_theta initial value (>0) for theta used for iteration 
#' @keywords Xu2010
#' @export
#' @examples data=meta.gen(500,theta_0=0.5,theta_1=0.5,L1=1,L2=1,L3=1,b01=0,b02=0,b03=0,cc=2,dd="uniform")
#' @examples 
#' @examples ans=Xu2010(data$X1,data$X2,data$D, int_theta=1 ,tol=0.01,step=50)
#' @examples ans

Xu2010=function(T1,T2,d2,int_theta,tol=tol,step){
d1=(T1<T2)*1
m=length(T1)
T2M=matrix(T2,m,m)
ST1=sort(T1,index.return=TRUE)
st1=ST1$x
T1M=matrix(T1,m,m)
ST1M=t(matrix(st1,m,m))
Risk1=ifelse(T1M>=ST1M,1,0)
sd1=d1[ST1$ix]
d1i=ifelse(T1M==ST1M,1,0)*matrix(d1,m,m)
ST23=sort(T2,index.return=TRUE)
sd23=( d2*(1-d1) )[ST23$ix]
st23=ST23$x
ST23M=t(matrix(st23,m,m))
Risk2=ifelse(T1M>=ST23M,1,0)
d23i=ifelse(T2M==ST23M,1,0)*matrix(d2*(1-d1),m,m)
ST24=sort(T2,index.return=TRUE)
sd24=(d1*d2)[ST24$ix]
st24=ST24$x
ST24M=t(matrix(st24,m,m))
Risk12=ifelse((T2M>=ST24M) & (ST24M>T1M),1,0)

d24i=ifelse(T2M==ST24M,1,0)*matrix(d1*d2,m,m)

pLA=function(yy,tt,LL){
pla=function(yi){
loc=sum(yi>=tt)
if(loc==0)ans=0
if(loc>0)ans=LL[loc]
ans
}
apply(matrix(yy),1,pla)
}

Nelson_Aalen=function(yy,dd,risk_na){
down_na=apply(risk_na,2,sum)
down_na=ifelse(down_na==0,1,down_na)
dL=dd/down_na
dL
}
dL1=Nelson_Aalen(st1,sd1,Risk1)
dL2=Nelson_Aalen(st23,sd23,Risk2)
dL3=Nelson_Aalen(st24,sd24,Risk12)


LA1=pLA(T1,st1,cumsum(dL1))
LA2=pLA(T1,st23,cumsum(dL2))
LA32=pLA(T2,st24,cumsum(dL3))
LA31=pLA(T1,st24,cumsum(dL3))
LA3=LA32-LA31
Ai=LA1+LA2+LA3
theta=int_theta
ETA=c(theta,dL1,dL2,dL3) #¼ȩw

Score=function(ETA){
theta=ETA[1]
dL1=ETA[2:(m+1)]
dL2=ETA[(m+2)  :(2*m+1)]
dL3=ETA[(2*m+2):(3*m+1)]

Bi=1/theta+(d1+d2)

LA1=pLA(T1,st1,cumsum(dL1))
LA2=pLA(T1,st23,cumsum(dL2))
LA32=pLA(T2,st24,cumsum(dL3))
LA31=pLA(T1,st24,cumsum(dL3))
LA3=LA32-LA31
Ai=LA1+LA2+LA3


down1=apply(Risk1*(1+theta*d1+theta*d2)/(1+theta*Ai),2,sum)
down2=apply(Risk2*(1+theta*d1+theta*d2)/(1+theta*Ai),2,sum)
down3=apply(Risk12*(1+theta*d1+theta*d2)/(1+theta*Ai),2,sum)


U1=function(x){
Bi=1/x+d1+d2
 
#candA=which(1+x*Ai>0)
#replaceA=max(Ai[candA])
#Ai[which(1+x*Ai<=0)]=replaceA

sum(
          d1*d2/(1+x)+
          1/(x^2)*log(1+x*Ai)-
          Bi*Ai/(1+x*Ai) 
    )
}
u1=U1(theta)
u2=sd1/dL1-down1
u3=sd23/dL2-down2
u4=sd24/dL3-down3
ans=c(u1,u2,u3,u4)
ans
}
ETA=c(theta,dL1,dL2,dL3) 

SS_1=Score(ETA)
SSS_1=abs(sum(SS_1[SS_1!="NaN"]))
SSS_0=10000

theta=theta_1=int_theta
kkk=0
ETA_1=c(int_theta,dL1,dL2,dL3)
ETA_0=ETA_1*0
dist_0=1000
dist_1=100
USp=history=US=Ans_ETA=NULL
dist_0=10000
SS_1=Score(ETA_1)
can_to_stop=c(1,1+(1:10),m+(1:10),2*m+(1:10))
dist_1= max(abs(ETA_0[can_to_stop]-ETA_1[can_to_stop]))
divk=0

usp=max(abs(SS_1[SS_1!="NaN"]))

hist_U=NULL
st_c=0
Ans_ETA=NULL
diverge="F"
#print(c(st_c,theta ,usp,tol , dist_1,tol, kkk))
while(st_c<10&
theta>10^-4&
  (usp>tol | dist_1>tol) &
kkk<step
){

if(dist_1-dist_0>0){
divk=divk+1
divk=max(divk,3)
theta=theta_1=int_theta*0.1^divk
kkk=0
ETA_1=c(int_theta,dL1,dL2,dL3)
ETA_0=ETA_1*0
dist_0=1000
dist_1=100
USp=history=US
dist_0=10000
SS_1=Score(ETA_1)
can_to_stop=c(1,1+(1:10),m+(1:10),2*m+(1:10))
}

SSS_0=SSS_1
Bi=1/theta+(d1+d2)
rr=30
rrr=15
#rr=min(50,rr+divk)

for(jj in 1:rr){
## U2
kkrr=0
#while(prod(1 + theta * Ai)>0 & kkrr<rrr){
for(i in 1:rrr){
kkrr=kkrr+1
LA1=pLA(T1,st1,cumsum(dL1))
Ai=LA1+LA2+LA3
down1=apply(Risk1*(1+theta*d1+theta*d2)/(1+theta*Ai),2,sum)
down1=ifelse(down1==0,1,down1)
dL1=sd1/down1
}
#
# U3
for(i in 1:rrr){
LA2=pLA(T1,st23,cumsum(dL2))
Ai=LA1+LA2+LA3

down2=apply(Risk2*(1+theta*d1+theta*d2)/(1+theta*Ai),2,sum)
down2=ifelse(down2==0,1,down2)
dL2=sd23/down2
}
for(i in 1:rrr){
LA32=pLA(T2,st24,cumsum(dL3))
LA31=pLA(T1,st24,cumsum(dL3))
LA3=LA32-LA31
Ai=LA1+LA2+LA3
down3=apply(Risk12*(1+theta*d1+theta*d2)/(1+theta*Ai),2,sum)
down3=ifelse(down3==0,1,down3)
dL3=sd24/down3
}
}

U1=function(x){
Bi=1/x+d1+d2
sum(
          d1*d2/(1+x)+
          1/(x^2)*log(1+x*Ai)-
          Bi*Ai/(1+x*Ai) 
    )
}

U1_2=function(x){
#x=theta
Bi=1/x+d1+d2

sum(
-d1*d2/(1+x)^2-2/x^3*log(1+x*Ai)+2*Ai/x^2/(1+x*Ai)+Bi*(Ai/(1+x*Ai))^2
)
}
U_st=1
if(theta>10^-5){
while(abs(U1(theta))>0.0000001&U_st==1){
if(theta_1>0) {theta_1=theta_1-1/U1_2(theta_1)*U1(theta_1)}
if(theta_1<0) {U_st=0}
if(theta_1>0) theta=theta_1 else theta=10^-6
}}

ETA_0=ETA_1
ETA_1=c(theta,dL1,dL2,dL3)
SS_1=Score(ETA_1)
SSS_1=sum(abs(SS_1[SS_1!="NaN"]))

kkk=kkk+1

dist_0=dist_1
dist_1=max(abs(ETA_0[can_to_stop]-ETA_1[can_to_stop]))

ETA=c(theta,dL1,dL2,dL3)

Ans_ETA=rbind(Ans_ETA,ETA )

hist_U=c(hist_U,max(abs(SS_1),na.rm=TRUE))

if(length(hist_U)>7){
st_c=sum(diff(hist_U[-5])>0)
}
#print(theta)
usp=max(abs(SS_1[SS_1!="NaN"]))
}
cho_ans=which(hist_U==min(hist_U))
ETA=Ans_ETA[cho_ans,]
if(theta==10^-6)
diverge="T"
a0=length(dL1)
a1=length(dL2)
a2=length(dL3)
x=theta=ETA[1]
dL1=ETA[(1+1):(1+a0)]
dL2=ETA[(1+a0+1):(1+a0+a1)]
dL3=ETA[(1+a0+a1+1):(1+a0+a1+a2)]

#U1_theta

Bi=1/x+d1+d2
LA1=pLA(T1,st1,cumsum(dL1))
LA2=pLA(T1,st23,cumsum(dL2))
LA32=pLA(T2,st24,cumsum(dL3))
LA31=pLA(T1,st24,cumsum(dL3))
LA3=LA32-LA31
Ai=LA1+LA2+LA3
SS=Score(ETA)
LPC=which(ETA!=0)
eee=(1/m)^2
ee=diag(eee,length(ETA))
UU=NULL
app_score=function(ii){( Score(ETA+ee[ii,])-SS )/eee}
uu=apply(matrix(LPC),1,app_score)
selectNA=apply(uu,1,sum)
UU=uu[selectNA !="NaN",]
UU=0.5*(UU+t(UU))
JM=solve(-UU)
sdd=diag(JM)^0.5
j1=sum(dL1>0)
j2=sum(dL2>0)
j3=sum(dL3>0)
theta_sd=sdd[1]
sdd=sdd[-1]
dLA1_sd=sdd[1:j1]
sdd=sdd[-(1:j1)]
if(j2>0){dLA2_sd=sdd[1:j2]
          sdd=sdd[-(1:j2)]
          dLA3_sd=sdd
}
if(j2==0){dLA2_sd=0
          sdd=sdd[-(1:j1)]
          dLA3_sd=sdd
}

sd.err=diag(JM)^0.5

a_1=sum(dL1!=0)
a_2=sum(dL2!=0)
a_3=sum(dL3!=0)
report=list(
theta=theta,
dL1=dL1,
dL2=dL2,
dL3=dL3,
st1=st1,st2=st23,
d1=sd1,d23=sd23,d24=sd24,
JM=JM,
theta_sd=sd.err[1] ,
diverge=diverge
)
return(report)
}

