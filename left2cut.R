left2cut<-function(b,frm,form,r,pmax,p,pold,alpha,iter){
while (sum(pold==p)<r+1)
{
  form<-as.formula(form)# formula for modeling
  mdl<-lm(form,data=b)# linear regression model
  c<-summary(mdl)$coefficients
  pold<-p
  z<-0
  for (i in 1:(r+1))
    if (p[i]>0)
    {
      t<-p[i]
      while (c[z+t,4]>alpha)
      {
	t<-t-1
	if (t==0)
	  break
      }
      z<-z+p[i]
      p[i]<-t
    }
  if (p[r+1]==0)
    p[r+1]<-order(c[(sum(pold[1:r])+1):sum(pold),4])[1]
  form<-frm
  for (i in 1:(r+1))
    if (p[i]>0)
      for (j in 1:p[i])
	form<-paste(form,"+b[,",(i-1)*pmax+j,"]",sep="")
  iter<-iter+1
}
return(list(model=mdl,vect=p,iter=iter))
}