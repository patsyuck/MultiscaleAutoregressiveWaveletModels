l.filter<-function(x,n,pmax){
if (length(x)==1)
  lx<-x
else
  lx<-length(x)
floor((lx-(pmax-1)*2^n-1-pmax*(n+1))/(2^n-1)+1)
}