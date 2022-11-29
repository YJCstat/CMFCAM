#' generating the data with 1/2 bivariate exposure from the frailty model by  three exponential distributions
#' @param n sample size (even)
#' @param theta_0 theta for Z=1
#' @param theta_1 theta for Z=2
#' @param L1 lambda for genarate T1
#' @param L2 lambda for generate T2 without given T1
#' @param L3 lambda for generate T2 given T1
#' @param b01 effect from Z to T1
#' @param b02 effect from Z to T2
#' @param b03 effect from T1 to T2
#' @param cc parameter for generating the censoring time, the regulator censoring rate 
#' @param dd  set "uniform" for U(0,cc); set weibull for weibull(shape=5,scale=cc)
#' @param output X1,X2,Z and D are observed mediated, terminal event times exposure and censoring index (1/0 for failure and censored)
#' @keywords meta.gen
#' @export
#' @examples 
#' @examples meta.gen(500,theta_0=1,theta_1=0.5,L1=1,L2=1,L3=1,b01=0.5,b02=0,b03=1,cc=2,dd="uniform")



meta.gen=function(n,theta_0,theta_1,L1,L2,L3,b01,b02,b03,cc=2,dd="uniform"){
#n=500 ;theta_0=0;theta_1=0.5 ; b01=0;b02=0.5;b03=0.25;L1=2.5;L2=2;L3=2.5;cc=0.75;dd="weibull"
  Z_1=rep(1,times=n/2)
 U1_1=runif(n/2)
 U2_1=runif(n/2)
U12_1=runif(n/2)
 T1_1=T2_1=rep(NA,n/2)
if(theta_0==0) rr_1=rep(1,n/2) else rr_1=rgamma(n/2,1/theta_0,1/theta_0)

T1_1=-log(U1_1)/(L1*rr_1*exp(b01*Z_1))
T2_1=-log(U2_1)/(L2*rr_1*exp(b02*Z_1))
n1_1=which(T2_1>T1_1)
T2u_1=-log(U2_1)-T1_1*rr_1*( L2*exp(b02*Z_1)-L3*exp(b03*Z_1) )
T2_1[n1_1]=( T2u_1/(L3*rr_1*exp(b03*Z_1)) )[n1_1]
T1_1[-n1_1]=10000000
T10_1=T1_1
if(dd=="weibull") C_1=rweibull(n/2,5,cc)
if(dd=="uniform") C_1=runif(n/2,0,cc)
Cm_1=ifelse(C_1<=T2_1,C_1,T2_1)
D1_1=ifelse(T1_1>Cm_1,0,1)
T1_1=ifelse(T1_1>Cm_1,Cm_1,T1_1)
D2_1=ifelse(T2_1>C_1,0,1)
T2_1=ifelse(T2_1>C_1,C_1,T2_1)
################################################
  Z_2=rep(2,times=n/2)
 U1_2=runif(n/2)
 U2_2=runif(n/2)
U12_2=runif(n/2)
 T1_2=T2_2=rep(NA,n/2)

if(theta_1==0) rr_2=rep(1,n/2) else rr_2=rgamma(n/2,1/theta_1,1/theta_1)

T1_2=-log(U1_2)/(L1*rr_2*exp(b01*Z_2))
T2_2=-log(U2_2)/(L2*rr_2*exp(b02*Z_2))
n1_2=which(T2_2>T1_2)
T2u_2=-log(U2_2)-T1_2*rr_2*( L2*exp(b02*Z_2)-L3*exp(b03*Z_2) )
T2_2[n1_2]=( T2u_2/(L3*rr_2*exp(b03*Z_2)) )[n1_2]
T1_2[-n1_2]=10000000
T10_2=T1_2
if(dd=="weibull") C_2=rweibull(n/2,5,cc)
if(dd=="uniform") C_2=runif(n/2,0,cc)
Cm_2=ifelse(C_2<=T2_2,C_2,T2_2)
D1_2=ifelse(T1_2>Cm_2,0,1)
T1_2=ifelse(T1_2>Cm_2,Cm_2,T1_2)
D2_2=ifelse(T2_2>C_2,0,1)
T2_2=ifelse(T2_2>C_2,C_2,T2_2)

################################################
X1=c(T1_1,T1_2)
X2=c(T2_1,T2_2)
D=c(D2_1,D2_2)
Z=c(Z_1,Z_2)
data.frame(X1,X2,Z,D)
}
