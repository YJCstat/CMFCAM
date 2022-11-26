#' predicting the value for one to one vector
#' @param yy,tt,LL
#' @param yy any time you want to interpolate
#' @param tt time or x value with same length of LL, must >0
#' @param LL value as f(tt)
#' @keywords pLA
#' @export
#' @examples x=seq(0,0.5,by=0.01)
#' LL=pnorm(x,0,1)
#' pLA(c(0,0.1,0.25,0.01,3),x,LL)


pLA=function(yy,tt,LL){
pla=function(yi){
loc=sum(yi>=tt)
if(loc==0)ans=0
if(loc>0)ans=LL[loc]
ans
}
apply(matrix(yy),1,pla)
}
