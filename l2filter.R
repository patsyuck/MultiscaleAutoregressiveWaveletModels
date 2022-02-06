l2filter<-function(x,n,r,pmax){
if (length(x)==1)
  lx<-x
else
  lx<-length(x)
floor((lx-(pmax-1)*2^r-1-pmax*(r+1))/(2^n-1)+1)
}