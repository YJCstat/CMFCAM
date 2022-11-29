#' Estimating the direct and indirect of the Multistate model. Obtaining the variance by bootstrapping
#' @import Copula.surv
#' @import survival
#' @import BSDA
#' @import mstate
#' @param data data.frame(X1,X2,D,Z)
#' @param P.time interpolation time can be vector or scalar
#' @keywords causal inference, semicompeting risks, frailty model
#' @export
#' @examples data=meta.gen(500,theta_0=0.5,theta_1=0.5,L1=0.5,L2=0.5,L3=1,b01=1,b02=0,b03=0,cc=2,dd="uniform")
#' @examples P.time=seq(0,1,by=0.01)
#' @examples ans=M_state(data,P.time)
#' @examples plot(P.time,ans$DE,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$DE+ans$DE_sd,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$DE-ans$DE_sd,type="l",ylim=c(-0.5,0.5))
#' @examples points(P.time,ans$IE,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples points(P.time,ans$IE+ans$IE_sd,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples points(P.time,ans$IE-ans$IE_sd,type="l",ylim=c(-0.5,0.5),col=2)
#' @examples legend(0,0.45,c("direct effect","indirect effect"),col=1:2,lty=1)


M_state=function(data,P.time){
Result=ms(data,P.time)
ResultBS=list()
n=length(data$X1)
DEBS=IEBS=NULL
for(bb in 1:500){
dataBB=data[sample(1:n,n,replace = TRUE),]

ResultBS=ms(dataBB,P.time)
IEBS=rbind(IEBS,ResultBS$IE)
DEBS=rbind(DEBS,ResultBS$DE)

}
Ans=list(
DE=Result$DE,
DE_sd=apply(DEBS,2,sd,na.rm=TRUE),
IE=Result$IE,
IE_sd=apply(IEBS,2,sd,na.rm=TRUE)
)
return(Ans)
}