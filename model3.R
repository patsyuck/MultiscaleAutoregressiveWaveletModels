model3<-function(x,filter="haar",n=5,r=3,pmax=3,len=20){
lmax<-l2filter(x,n,r,pmax)
dec<-holesDWT(x,filter,n,r,pmax,len,lmax)# multiresolution decomposition
day<-dec$days-len
b<-r.matrix(dec,r,pmax,day,len)# matrix for modeling
frm<-paste("b[,",(r+1)*pmax+1,"]~0",sep="")
ARSmax<-0
vect<-vector(mode="numeric",length=r+1)
for (k in 1:pmax)
{
  p<-vector(mode="numeric",length=r+1)
  p[r+1]<-k
  while (sum(p[1:r])<(r*pmax))
  {
    for(i in 1:r)
      if (p[i]<pmax)
      {
	p[i]<-p[i]+1
	break
      }
      else
	p[i]<-0
    form<-frm
    for (i in 1:(r+1))
      if (p[i]>0)
	for (j in 1:p[i])
	  form<-paste(form,"+b[,",(i-1)*pmax+j,"]",sep="")
    form<-as.formula(form)# formula for modeling
    mdl<-lm(form,data=b)# linear regression model
    ARS<-summary(mdl)$adj.r.squared
    if (ARS>ARSmax)
    {
      ARSmax<-ARS
      vect<-p
    }
  }
}
return(vect)
}