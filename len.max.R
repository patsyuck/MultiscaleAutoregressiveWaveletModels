len.max<-function(lx,l,n,r,pmax){
return(lx-(l-1)*(2^n-1)-(pmax-1)*2^r-1)
}