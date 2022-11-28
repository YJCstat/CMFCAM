library("Copula.surv")
Copula.DEIE=function(data,P.time){
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
dz1=which(data$X1<data$X2)
dz[dz1]=1
z.obs_0=z.obs[Z==0]
z.obs_1=z.obs[Z==1]
dz_0=dz[Z==0]
dz_1=dz[Z==1]
###################################################################

theta_0=U2.Clayton(x.obs_0,y.obs_0,dx_0,dy_0)


theta_1=U2.Clayton(x.obs_1,y.obs_1,dx_1,dy_1)


###################################################################


Sz0=survfit(Surv(z.obs_0,dz_0) ~ 1, conf.type = "log-log")
Sz1=survfit(Surv(z.obs_1,dz_1) ~ 1, conf.type = "log-log")
Sy0=survfit(Surv(y.obs_0,dy_0) ~ 1, conf.type = "log-log")
Sy1=survfit(Surv(y.obs_1,dy_1) ~ 1, conf.type = "log-log")

#########################
P.time=as.matrix(P.time)


predict_z0=function(x){
L=sum(x>=Sz0$time) 
if(L==0) ans=1
if(L>0) ans=Sz0$surv[L]
if(L>=length(Sz0$surv)) ans=Sz0$surv[  max(which(Sz0$surv>0))]
ans
}

Sz0.y=apply(P.time,1,predict_z0)
################################

predict_z1=function(x){
L=sum(x>=Sz1$time) 
if(L==0) ans=1
if(L>0) ans=Sz1$surv[L]
if(L>=length(Sz1$surv)) ans=Sz1$surv[ max(which(Sz1$surv>0))]
ans
}

Sz1.y=apply(P.time,1,predict_z1)
################################

predict_y0=function(x){
L=sum(x>=Sy0$time) 
if(L==0) ans=1
if(L>0) ans=Sy0$surv[L]
if(L>=length(Sy0$surv)) ans=Sy0$surv[ max(which(Sy0$surv>0))]
ans
}

Sy0.y=apply(P.time,1,predict_y0)
################################

predict_y1=function(x){
L=sum(x>=Sy1$time) 
if(L==0) ans=1
if(L>0) ans=Sy1$surv[L]
if(L>=length(Sy1$surv)) ans=Sy1$surv[max(which(Sy1$surv>0))]
ans
}

Sy1.y=apply(P.time,1,predict_y1)
################################

Sx0.y=(Sz0.y^(-theta_0[1])-Sy0.y^(-theta_0[1])+1 )^(1/(-theta_0[1]))
Sx1.y=(Sz1.y^(-theta_1[1])-Sy1.y^(-theta_1[1])+1 )^(1/(-theta_1[1]))

#plot(P.time,Sx0.y,type="l")
#plot(P.time,Sx1.y,type="l")
#plot(P.time,Sz0.y^(1-theta_0[1]),type="l")
#plot(P.time,Sy0.y^(1-theta_0[1]),type="l")
#plot(P.time,Sz0.y^(1-theta_0[1])-Sy0.y^(1-theta_0[1]),type="l")


W1tz=function(zz){
if(zz==0){ theta_z=theta_0[1];Sc1=Sx0.y;Sc2=Sy0.y }
if(zz==1){ theta_z=theta_1[1];Sc1=Sx1.y;Sc2=Sy1.y }
1-(Sc1^(-theta_z)+Sc2^(-theta_z)-1)^(-1/theta_z) /Sc2
} 

lambda_0=function(zz){
if(zz==0){ theta_z=theta_0[1];Sc1=Sx0.y;Sc2=Sy0.y }
if(zz==1){ theta_z=theta_1[1];Sc1=Sx1.y;Sc2=Sy1.y }

LAc2=-log(Sc2)
dLac2=diff(c(0,LAc2))
ans=dLac2*Sc2^(-theta_z)/(Sc1^(-theta_z)+Sc2^(-theta_z)-1)
ifelse(ans=="NaN",0,ans)
} 

lambda_1=function(zz){
if(zz==0){ theta_z=theta_0[1];Sc1=Sx0.y;Sc2=Sy0.y }
if(zz==1){ theta_z=theta_1[1];Sc1=Sx1.y;Sc2=Sy1.y }

LAc2=-log(Sc2)
dLac2=diff(c(0,LAc2))

RP=(
Sc2-Sc2^(-theta_z)*(Sc1^(-theta_z)+Sc2^(-theta_z)-1)^(-1/theta_z-1)
)/(
Sc2-(Sc1^(-theta_z)+Sc2^(-theta_z)-1)^(-1/theta_z)
)
RP=ifelse(RP=="NaN",0,RP)
RP=ifelse(abs(RP)==Inf,0,RP)

dLac2*RP

}


CountLA=function(za,zb){
W1tz(zb)*lambda_1(za)+(1-W1tz(zb))*lambda_0(za)
#W1tz(1)
#lambda_1(1)
#lambda_0(1)
}

DE=cumsum(CountLA(1,0)-CountLA(0,0))
IE=cumsum(CountLA(1,1)-CountLA(1,0))


Report=list(
DE=DE,
IE=IE,
Sx1z0=Sx0.y,
Sx1z1=Sx1.y,
Sy1z0=Sy0.y,
Sy1z1=Sy1.y,
theta=c(theta_0[1],theta_1[1])
)
 
return(Report)
}
