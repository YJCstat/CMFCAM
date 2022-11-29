
#' Estimating the direct and indirect of the Multistate model.
#' @import Copula.surv
#' @import survival
#' @import BSDA
#' @param data data.frame(X1,X2,D,Z)
#' @param P.time interpolation time can be vector or scalar
#' @keywords causal inference, semicompeting risks, frailty model
#' @export
#' @examples #install.packages("mstate")
#' @examples library(mstate)
#' @examples data=meta.gen(500,theta_0=0.5,theta_1=0.5,L1=0.5,L2=0.5,L3=1,b01=1,b02=0,b03=0,cc=2,dd="uniform")
#' @examples P.time=seq(0,1,by=0.01)
#' @examples ans=ms(data,P.time)
#' @examples plot(P.time,ans$DE,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$IE,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples legend(0,0.45,c("direct effect","indirect effect"),col=1:2,lty=1)
ms=function(data,P.time){

typeS=names(table(data$Z))
data$S=data$Z
Z=data$Z
Z[data$S==typeS[1]]=0
Z[data$S==typeS[2]]=1
x.obs_0=data$X1[Z==0]
y.obs_0=data$X2[Z==0]
dy_0=data$D[Z==0]
dx_0=ifelse(x.obs_0==y.obs_0,0,1)
#---------------------------------------------------
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
#########################################################

###################################################################
Ms.haz=function(data,mediated,primary,dx,dy){

tmat <- transMat(
                 x = list(c(2,3),c(3),c()),
                 names = c("Tx",mediated,primary)
                 )
msebmt <- msprep(
data   = data,
trans  = tmat, 
time   = c(NA,mediated, primary), 
status = c(NA,dx, dy) 
                )
c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
            data = msebmt,
            method = "breslow"
            )
msf0 <- msfit(object = c0, vartype = "greenwood", trans = tmat)
msf0
}
###################################################################
data_0=data.frame(x.obs_0,y.obs_0,dx_0,dy_0)
data_1=data.frame(x.obs_1,y.obs_1,dx_1,dy_1)

haz0=Ms.haz(data_0,"x.obs_0","y.obs_0","dx_0","dy_0")$Haz
haz1=Ms.haz(data_1,"x.obs_1","y.obs_1","dx_1","dy_1")$Haz
###################################################################
time0=unique(haz0[,1])
haz0_SD=haz0[haz0[,3]==2,2]
haz0_SM=haz0[haz0[,3]==1,2]
haz0_MD=haz0[haz0[,3]==3,2]

time1=unique(haz1[,1])
haz1_SD=haz1[haz1[,3]==2,2]
haz1_SM=haz1[haz1[,3]==1,2]
haz1_MD=haz1[haz1[,3]==3,2]
#########################
P.time=as.matrix(P.time)
predict_z=function(x,time,R){
L=sum(x>=time) 
if(L==0) ans=0
if(L>0) ans=R[L]
if(L>=length(R)) ans=R[length(R)]
ans
}
################################
HSD_0=function(x){predict_z(x,time0,haz0_SD)}
haz0_SD_P=apply(P.time,1,HSD_0)
################################
HSM_0=function(x){predict_z(x,time0,haz0_SM)}
haz0_SM_P=apply(P.time,1,HSM_0)
################################
HMD_0=function(x){predict_z(x,time0,haz0_MD)}
haz0_MD_P=apply(P.time,1,HMD_0)
################################
################################
HSD_1=function(x){predict_z(x,time1,haz1_SD)}
haz1_SD_P=apply(P.time,1,HSD_1)
################################
HSM_1=function(x){predict_z(x,time1,haz1_SM)}
haz1_SM_P=apply(P.time,1,HSM_1)
################################
HMD_1=function(x){predict_z(x,time1,haz1_MD)}
haz1_MD_P=apply(P.time,1,HMD_1)
################################

W1tz=function(zz){
if(zz==0){wn1=1-exp(-haz0_SM_P)}
if(zz==1){wn1=1-exp(-haz1_SM_P)}
wn1
} 

lambda_0=function(zz){
if(zz==0){dL0=diff(c(0,haz0_SD_P)) }
if(zz==1){dL0=diff(c(0,haz1_SD_P)) }
ans=dL0
ifelse(ans=="NaN",0,ans)
} 

lambda_1=function(zz){
if(zz==0){dL1=diff(c(0,haz0_MD_P)) }
if(zz==1){dL1=diff(c(0,haz1_MD_P)) }
ans=dL1
ifelse(ans=="NaN",0,ans)
} 
#############################################

CountLA=function(za,zb){
W1tz(zb)*lambda_1(za)+(1-W1tz(zb))*lambda_0(za)
}

DE=cumsum(CountLA(1,0)-CountLA(0,0))
IE=cumsum(CountLA(1,1)-CountLA(1,0))

Report=list(
DE=DE,
IE=IE,
LA11=haz0_SM_P,
LA12=haz0_SD_P,
LA13=haz0_MD_P,
LA21=haz1_SM_P,
LA22=haz1_SD_P,
LA23=haz1_MD_P
)

return(Report)
}
