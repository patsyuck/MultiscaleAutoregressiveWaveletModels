r2matrix<-function(dec,r,p,top,day,len){
if (length(p)==1)
{
  b<-matrix(nrow=len,ncol=(r+1)*p+1)# box of matrix for modeling
  for (i in 1:len)
  {
    b[i,(r+1)*p+1]<-dec$appr[top+1,day+i]
    for (j in 1:r)
      for (k in 1:p)
	b[i,(j-1)*p+k]<-dec$detal[j,day+i-1-(k-1)*2^j]
    for (k in 1:p)
      b[i,r*p+k]<-dec$appr[r+1,day+i-1-(k-1)*2^r]
  }
}
else
{
  b<-matrix(nrow=len,ncol=sum(p)+1)# box of matrix for modeling
  for (i in 1:len)
  {
    b[i,sum(p)+1]<-dec$appr[top+1,day+i]
    s<-0
    for (j in 1:r)
      if (p[j]>0)
      {
	for (k in 1:p[j])
	  b[i,s+k]<-dec$detal[j,day+i-1-(k-1)*2^j]
	s<-s+p[j]
      }
    if (p[r+1]>0)
      for (k in 1:p[r+1])
	b[i,s+k]<-dec$appr[r+1,day+i-1-(k-1)*2^r]
  }
}
return(as.data.frame(b))# data for modelling
}