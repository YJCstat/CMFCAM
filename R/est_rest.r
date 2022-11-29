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
#' @examples ans=Xu2010_rest(data$X1,data$X2,data$D, int_theta=1 ,tol=0.01,step=50)
#' @examples ans

Xu2010_rest=function(T1,T2,d2,int_theta,tol=0.01,step){
FIG="FALSE"
d1=(T1<T2)*1
m=length(T1)

T2M=matrix(T2,m,m)
#######for U2
ST1=sort(T1,index.return=TRUE)
st1=ST1$x
T1M=matrix(T1,m,m)
ST1M=t(matrix(st1,m,m))
Risk1=ifelse(T1M>=ST1M,1,0)

sd1=d1[ST1$ix]

###########for U3
ST23=sort(T2,index.return=TRUE)
sd23=d2[ST23$ix]

st23=ST23$x
ST23M=t(matrix(st23,m,m))
Risk2=ifelse(T2M>=ST23M,1,0)

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
##########################################
Score=function(ETA){
theta=ETA[1]
dL1=ETA[2:(m+1)]
dL2=ETA[(m+2)  :(2*m+1)]
Bi=1/theta+(d1+d2)

LA1=pLA(T1,st1,cumsum(dL1))
LA2=pLA(T2,st23,cumsum(dL2))
Ai=LA1+LA2


RRt10=Risk1*(1+theta*d1+theta*d2)/(1+theta*Ai)
RR010=Risk1*1
RRin1=ifelse(RRt10>0,RRt10,RR010)
down1=apply(RRin1,2,sum)
down1=ifelse(down1<=0,down10,down1)


RRt20=Risk2*(1+theta*d1+theta*d2)/(1+theta*Ai)
RR020=Risk2*1
RRin2=ifelse(RRt20>0,RRt20,RR020)
down2=apply(RRin2,2,sum)
down2=ifelse(down2<=0,down20,down2)
U1=function(x){
Bi=1/x+d1+d2
 
candA=which(1+x*Ai>0)
replaceA=max(Ai[candA])
Ai[which(1+x*Ai<=0)]=replaceA

sum(
          d1*d2/(1+x)+
          1/(x^2)*log(1+x*Ai)-
          Bi*Ai/(1+x*Ai) 
    )
}
u1=U1(theta)
u2=sd1/dL1-down1
u3=sd23/dL2-down2
c(u1,u2,u3)
}

################################################################
ETA=c(int_theta,dL1,dL2) 

SS_1=Score(ETA)
SSS_1=abs(sum(SS_1[SS_1!="NaN"]))
SSS_0=10000
#################################################################
#liter liter....liter
#################################################################

LA1=pLA(T1,st1,cumsum(dL1))
LA2=pLA(T2,st23,cumsum(dL2))
Ai=LA1+LA2
theta=int_theta
theta_0=1000
theta_1=int_theta
kkk=0
ETA_1=c(theta,dL1,dL2)
ETA_0=ETA_1*0
USp=NULL
SS_1=Score(ETA_1)
usp=max(abs(SS_1[SS_1!="NaN"]))
can_to_stop=c(1,1+(1:10),m+(1:10))
dist_1= max(abs(ETA_0[can_to_stop]-ETA_1[can_to_stop]))
Ans_ETA=NULL
#############################################################
#          start for while                                  #
#############################################################
while(
theta>0&
  (usp>tol | dist_1>tol) &
kkk<step
){
theta_0=theta_1
Bi=1/theta+(d1+d2)
###########################################################
rrr=50                                                   #
## U2                                                     # 
for(i in 1:rrr){                                          #
LA1=pLA(T1,st1,cumsum(dL1))                               #
Ai=LA1+LA2                                                #     
                                                          # 
RRt10=Risk1*(1+theta*d1+theta*d2)/(1+theta*Ai)            #
RR010=Risk1*1                                             #
RRin1=ifelse(RRt10>0,RRt10,RR010)                         #
                                                          #
down1=apply(RRin1,2,sum)                                  #
down1=ifelse(down1<=0,down10,down1)                       #
dL1=sd1/down1                                             #
#}                                                        # 
# U3                                                      #
#for(i in 1:rrr){                                         #
LA2=pLA(T2,st23,cumsum(dL2))                              #
Ai=LA1+LA2                                                #
                                                          # 
RRt20=Risk2*(1+theta*d1+theta*d2)/(1+theta*Ai)            #
RR020=Risk2*1                                             #
RRin2=ifelse(RRt20>0,RRt20,RR020)                         #
down2=apply(RRin2,2,sum)                                  #  
down2=ifelse(down2<=0,down20,down2)                       # 
dL2=sd23/down2                                            #
#}                                                        #
}                                                         #
###########################################################

U1=function(x){
Bi=1/x+d1+d2
Bi=1/x+d1+d2
candA=which(1+x*Ai>0)
replaceA=max(Ai[candA])
Ai[which(1+x*Ai<=0)]=replaceA
sum(
          d1*d2/(1+x)+
          1/(x^2)*log(1+x*Ai)-
          Bi*Ai/(1+x*Ai) 
    )
}
###################################
###################################

U1_2=function(x){
#x=theta
Bi=1/x+d1+d2
candA=which(1+x*Ai>0)
replaceA=max(Ai[candA])
Ai[which(1+x*Ai<=0)]=replaceA
sum(
-d1*d2/(1+x)^2-2/x^3*log(1+x*Ai)+2*Ai/x^2/(1+x*Ai)+Bi*(Ai/(1+x*Ai))^2
)
}

##############################
##############################

while(abs(U1(theta_1))>0.0000001){
theta_1=theta_1-0.1/U1_2(theta_1)*U1(theta_1)
theta=theta_1
}

ETA_0=ETA_1
ETA_1=c(theta,dL1,dL2)
SS_1=Score(ETA_1)
SSS_1=abs(sum(SS_1[SS_1!="NaN"]))
kkk=kkk+1

dist_0=dist_1
dist_1=max(abs(ETA_0[can_to_stop]-ETA_1[can_to_stop]))

#print(c(kkk,theta,abs(ETA_0[1]-ETA_1[1])))
#############################################
#plot(st1,cumsum(dL1),ylim=c(0,max(T2)))

ETA=c(theta,dL1,dL2)
usp=max(abs(SS_1[SS_1!="NaN"]))
USp=c(USp,usp)

renewTheta=c(theta,dL1,dL2)

Ans_ETA=rbind(Ans_ETA,c(theta,dL1,dL2) )
#print(c(abs(SS_1),theta,dist_1,dist_1-dist_0,divk,rr))
if(FIG=="TRUE"){
plot(Ans_ETA[,1],USp,type="l")
abline(h=0)
}
}
#
#
#
#
#ETA= Ans_ETA[which(US==min(US))[1],]
ETA=renewTheta

#################################################################
SS=Score(ETA)
LPC=which(ETA!=0)
eee=(1/m)^2
ee=diag(eee,length(ETA))
UU=NULL
app_score=function(ii){( Score(ETA+ee[ii,])-SS )/eee}
uu=apply(matrix(LPC),1,app_score)
UU=uu[apply(uu,1,sum)!="NaN",]
JM=solve(-UU)
sdd=diag(JM)^0.5
j1=sum(dL1>0)
theta_sd=sdd[1]
sdd=sdd[-1]
dLA1_sd=sdd[1:j1]
sdd=sdd[-(1:j1)]
dLA2_sd=sdd

######################################################


report=list(
theta=theta,
dL1=dL1,
dL2=dL2,
st1=st1,
st2=st23,
d1=sd1,d23=sd23,
JM=JM,
theta_sd=theta_sd,
dLA1_sd=dLA1_sd,
dLA2_sd=dLA2_sd)
return(report)
}




