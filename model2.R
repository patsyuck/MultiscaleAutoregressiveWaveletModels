model2<-function(x,filter="haar",n=5,r=3,pmax=5,len=20,alpha=0.75){
lmax<-l2filter(x,n,r,pmax)
dec<-holesDWT(x,filter,n,r,pmax,len,lmax)# multiresolution decomposition
day<-dec$days-len
b<-r.matrix(dec,r,pmax,day,len)# matrix for modeling
frm<-paste("b[,",(r+1)*pmax+1,"]~0",sep="")
form<-frm
for (i in 1:((r+1)*pmax))
  form<-paste(form,"+b[,",i,"]",sep="")
form<-as.formula(form)# formula for modeling
mdl<-lm(form,data=b)# linear regression model
c<-cor(b)
form2<-frm
beta<-0
for (i in 1:((r+1)*pmax))
  if (abs(c[i,(r+1)*pmax+1])>=alpha)
  {
    form2<-paste(form2,"+b[,",i,"]",sep="")
    beta<-beta+1
  }
if (beta==0)
  stop("The given correlation level 'alpha' is very high.")
form2<-as.formula(form2)# formula for new modeling
mdl2<-lm(form2,data=b)# new linear regression model
index<-c(1:len)# set of labels of the time axis
ymin<-min(min(dec$appr[1,(day+1):dec$days]),min(fitted(mdl)),min(fitted(mdl2)))
ymax<-max(max(dec$appr[1,(day+1):dec$days]),max(fitted(mdl)),max(fitted(mdl2)))
plot(index,dec$appr[1,(day+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",pmax=",pmax,",a=",alpha,",b=",beta,sep=""))# the time series
lines(index,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphic of the fitted values of the model
lines(index,fitted(mdl2),xlab="",ylab="",type='l',col="green")# graphic of the fitted values of the new model
return(mdl2)
}