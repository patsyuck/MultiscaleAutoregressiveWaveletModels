d.train<-function(x,l,n,pmax){
if (length(x)==1)
  lx<-x
else
  lx<-length(x)
return(lx-(l-1)*(2^n-1)-(pmax-1)*2^n-1)
}