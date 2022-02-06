future<-function(dec,p,n,r,cf,fut,lz){
a<-matrix(nrow=n+1,ncol=lz)# box of the enlarged approximation matrix 
d<-matrix(nrow=n,ncol=lz)# box of the enlarged detalization matrix
a[,1:dec$days]<-dec$appr
d[,1:dec$days]<-dec$detal
for(i in 1:fut)
{
  day<-dec$days+i-1
  dc<-vector(mode="numeric",length=sum(p))
  u<-1
  for(j in 1:r)
    if(p[j]>0)
      for(k in 1:p[j])
      {
	dc[u]<-d[j,day-(k-1)*2^j]
	u<-u+1
      }
  for (k in 1:p[r+1])
  {
    dc[u]<-a[r+1,day-(k-1)*2^r]
    u<-u+1
  }
  a[1,day+1]<-sum(cf*dc)
  for(j in 1:n)
  {
    a[j+1,day+1]<-0
    for(k in 1:dec$lf)
      a[j+1,day+1]<-a[j+1,day+1]+dec$f[k]*a[j,day+1-(k-1)*2^(j-1)]
    d[j,day+1]<-a[j,day+1]-a[j+1,day+1]
  }
}# we have the enlarged approximation and detalization matrices
return(progn=a[1,(dec$days+1):lz])
}