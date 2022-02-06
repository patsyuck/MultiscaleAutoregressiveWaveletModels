adj.R.sq<-function(b,frm,r,pmax){
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