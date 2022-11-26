#' Estimating the direct and indirect of the frailty model
#' @import Copula.surv
#' @param data data.frame(X1,X2,D,Z)
#' @param interpolation time can be vector or scalar
#' @param int_theta initial value of theta for iteration, nonnegative values vector of length 2
#' @param int_theta initial value of theta for iteration, nonnegative values vector of length 2
#' @param tol maximum tolerance of change during the iteration
#' @param step maximum  number of the iteration
#' @keywords causal inference, semicompeting risks, frailty model
#' @export
#' @examples data=meta.gen(500,theta_0=0.5,theta_1=0.5,L1=0.5,L2=0.5,L3=1,b01=1,b02=0,b03=0,cc=2,dd="uniform")
#' @examples P.time=seq(0,1,by=0.01)
#' @examples ans=Frailty(data,P.time,int_theta=c(0.5,0.5),tol=0.01,step=50)
#' @examples plot(P.time,ans$DE,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$DE+ans$DE_sd,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$DE-ans$DE_sd,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$IE,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples points(P.time,ans$IE+ans$IE_sd,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples points(P.time,ans$IE-ans$IE_sd,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples legend(0,0.45,c("direct effect","indirect effect"),col=1:2,lty=1)


Frailty=function(data,P.time,int_theta,tol,step){
data$S=data$Z
TOL=tol
typeS=names(table(data$S))
Z=data$S
Z[data$S==typeS[1]]=0
Z[data$S==typeS[2]]=1
x.obs_0=data$X1[Z==0]
y.obs_0=data$X2[Z==0]
dy_0=data$D[Z==0]
dx_0=ifelse(x.obs_0==y.obs_0,0,1)
x.obs_1=data$X1[Z==1]
y.obs_1=data$X2[Z==1]
dy_1=data$D[Z==1]
dx_1=ifelse(x.obs_1==y.obs_1,0,1)

z.obs=ifelse(data$X1<=data$X2,data$X1,data$X2)
dz=data$D
z.obs_0=z.obs[Z==0]
z.obs_1=z.obs[Z==1]
dz_0=dz[Z==0]
dz_1=dz[Z==1]

if(length(int_theta)==0){
theta_0=U2.Clayton(x.obs_0,y.obs_0,dx_0,dy_0)
theta_1=U2.Clayton(x.obs_1,y.obs_1,dx_1,dy_1)
TH0=theta_0[1]
TH1=theta_1[1]
}
if(length(int_theta)!=0){
TH0=int_theta[1]
TH1=int_theta[2]
}
Xu_0=Xu2010(x.obs_0,y.obs_0,dy_0,TH0,TOL,step)
Xu_1=Xu2010(x.obs_1,y.obs_1,dy_1,TH1,TOL,step)

theta_z1=Xu_1$theta
theta_z0=Xu_0$theta

dL1_z0=Xu_0$dL1
dL1_z1=Xu_1$dL1

dL2_z0=Xu_0$dL2
dL2_z1=Xu_1$dL2

dL3_z0=Xu_0$dL3
dL3_z1=Xu_1$dL3


Effect=function(theta_z0,dL1_z0,dL2_z0,dL3_z0,theta_z1,dL1_z1,dL2_z1,dL3_z1){

A1_z00=pLA(Xu_0$st2,Xu_0$st1,cumsum(dL1_z0))
A2_z00=pLA(Xu_0$st2,Xu_0$st2,cumsum(dL2_z0))
A1_z11=pLA(Xu_1$st2,Xu_1$st1,cumsum(dL1_z1))
A2_z11=pLA(Xu_1$st2,Xu_1$st2,cumsum(dL2_z1))

A1_z01=pLA(Xu_1$st2,Xu_0$st1,cumsum(dL1_z0))
A2_z01=pLA(Xu_1$st2,Xu_0$st2,cumsum(dL2_z0))
A1_z10=pLA(Xu_0$st2,Xu_1$st1,cumsum(dL1_z1))
A2_z10=pLA(Xu_0$st2,Xu_1$st2,cumsum(dL2_z1))

w_z00=(1+theta_z0* A1_z00  )^(-1/theta_z0)
w_z11=(1+theta_z1* A1_z11  )^(-1/theta_z1)
w_z01=(1+theta_z0* A1_z01  )^(-1/theta_z0)
w_z10=(1+theta_z1* A1_z10  )^(-1/theta_z1)


dL0_z00=dL2_z0/(1+theta_z0*(A1_z00+A2_z00) )
dL0_z11=dL2_z1/(1+theta_z1*(A1_z11+A2_z11) )
dL0_z01=dL2_z0/(1+theta_z0*(A1_z10+A2_z10) )
dL0_z10=dL2_z1/(1+theta_z1*(A1_z01+A2_z01) )

dL1_z00=(1+theta_z0)*dL3_z0/(1+theta_z0*((A1_z00+A2_z00)) )
dL1_z11=(1+theta_z1)*dL3_z1/(1+theta_z1*((A1_z11+A2_z11)) )
dL1_z01=(1+theta_z0)*dL3_z0/(1+theta_z0*((A1_z10+A2_z10)) )
dL1_z10=(1+theta_z1)*dL3_z1/(1+theta_z1*((A1_z01+A2_z01)) )

DE_10=cumsum(dL0_z11*w_z01+dL1_z11*(1-w_z01))
DE_00=cumsum(dL0_z00*w_z00+dL1_z00*(1-w_z00))

IE_11=cumsum(dL0_z11*w_z11+dL1_z11*(1-w_z11))

DE=pLA(P.time,Xu_1$st2,DE_10)-pLA(P.time,Xu_0$st2,DE_00)
IE= pLA(P.time,Xu_1$st2,IE_11)-pLA(P.time,Xu_1$st2,DE_10)
data.frame(DE=DE,IE=IE)

}
ee=(1/length(data$X1))/10
Diff_DE=Diff_IE=NULL

DEIE=Effect(theta_z0,dL1_z0,dL2_z0,dL3_z0,theta_z1,dL1_z1,dL2_z1,dL3_z1)
DE=DEIE[,1]
IE=DEIE[,2]

DEIE_ee=Effect(theta_z0+ee,dL1_z0,dL2_z0,dL3_z0,theta_z1,dL1_z1,dL2_z1,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]

DE_theZ0=(DE_ee-DE)/ee
IE_theZ0=(IE_ee-IE)/ee

Diff_DE=rbind(Diff_DE,DE_theZ0)
Diff_IE=rbind(Diff_IE,IE_theZ0)

for(i in which(dL1_z0!=0)){
dL1_z0_process=dL1_z0
dL1_z0_process[i]=dL1_z0[i]+ee
DEIE_ee=Effect(theta_z0,dL1_z0_process,dL2_z0,dL3_z0,theta_z1,dL1_z1,dL2_z1,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_dL1_z0=(DE_ee-DE)/ee
IE_dL1_z0=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_dL1_z0)
Diff_IE=rbind(Diff_IE,IE_dL1_z0)
}

for(i in which(dL2_z0!=0)){
dL2_z0_process=dL2_z0
dL2_z0_process[i]=dL2_z0[i]+ee
DEIE_ee=Effect(theta_z0,dL1_z0,dL2_z0_process,dL3_z0,theta_z1,dL1_z1,dL2_z1,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_dL2_z0=(DE_ee-DE)/ee
IE_dL2_z0=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_dL2_z0)
Diff_IE=rbind(Diff_IE,IE_dL2_z0)
}
for(i in which(dL3_z0!=0)){
dL3_z0_process=dL3_z0
dL3_z0_process[i]=dL3_z0[i]+ee
DEIE_ee=Effect(theta_z0,dL1_z0,dL2_z0,dL3_z0_process,theta_z1,dL1_z1,dL2_z1,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_dL3_z0=(DE_ee-DE)/ee
IE_dL3_z0=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_dL3_z0)
Diff_IE=rbind(Diff_IE,IE_dL3_z0)
}

DEIE_ee=Effect(theta_z0,dL1_z0,dL2_z0,dL3_z0,theta_z1+ee,dL1_z1,dL2_z1,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_theZ1=(DE_ee-DE)/ee
IE_theZ1=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_theZ1)
Diff_IE=rbind(Diff_IE,IE_theZ1)


for(i in which(dL1_z1!=0)){
dL1_z1_process=dL1_z1
dL1_z1_process[i]=dL1_z1[i]+ee
DEIE_ee=Effect(theta_z0,dL1_z0,dL2_z0,dL3_z0,theta_z1,dL1_z1_process,dL2_z1,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_dL1_z1=(DE_ee-DE)/ee
IE_dL1_z1=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_dL1_z1)
Diff_IE=rbind(Diff_IE,IE_dL1_z1)
}


for(i in which(dL2_z1!=0)){
dL2_z1_process=dL2_z1
dL2_z1_process[i]=dL2_z1[i]+ee
DEIE_ee=Effect(theta_z0,dL1_z0,dL2_z0,dL3_z0,theta_z1,dL1_z1,dL2_z1_process,dL3_z1)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_dL1_z2=(DE_ee-DE)/ee
IE_dL1_z2=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_dL1_z2)
Diff_IE=rbind(Diff_IE,IE_dL1_z2)
}


for(i in which(dL3_z1!=0)){
dL3_z1_process=dL3_z1
dL3_z1_process[i]=dL3_z1[i]+ee
DEIE_ee=Effect(theta_z0,dL1_z0,dL2_z0,dL3_z0,theta_z1,dL1_z1,dL2_z1,dL3_z1_process)
DE_ee=DEIE_ee[,1]
IE_ee=DEIE_ee[,2]
DE_dL1_z3=(DE_ee-DE)/ee
IE_dL1_z3=(IE_ee-IE)/ee
Diff_DE=rbind(Diff_DE,DE_dL1_z3)
Diff_IE=rbind(Diff_IE,IE_dL1_z3)
}


JMZ0=Xu_0$JM
JMZ1=Xu_1$JM
LJM=dim(JMZ1)[1]+dim(JMZ0)[1]
Big_JM=matrix(0,LJM,LJM)

Big_JM[1:dim(JMZ0)[1],1:dim(JMZ0)[1]]=JMZ0
Big_JM[-(1:dim(JMZ0)[1]),-(1:dim(JMZ0)[1])]=JMZ1

DE_sd=IE_sd=NULL
for(ii in 1:dim(Diff_DE)[2]){
def=Diff_DE[,ii]
ief=Diff_IE[,ii]
DE_sd=c(DE_sd,t(def)%*%Big_JM%*%def)
IE_sd=c(IE_sd,t(ief)%*%Big_JM%*%ief)
}
DE_sd=DE_sd^0.5
IE_sd=IE_sd^0.5

 #DE_sd[which(DE_sd=="NaN")]=DE_sd[which(DE_sd=="NaN")[1]-1]
 #IE_sd[which(IE_sd=="NaN")]=IE_sd[which(IE_sd=="NaN")[1]-1]
gc()


#} #cal_sd: calculate or not the sd.
Report=list(
theta1=Xu_0$theta,
theta2=Xu_1$theta,
DE=DE,
IE=IE,
DE_sd=DE_sd,
IE_sd=IE_sd,
diverge_0=Xu_0$diverge,
diverge_1=Xu_1$diverge
)

 
return(Report)
}