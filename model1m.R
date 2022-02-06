model1m<-function(x,filter="haar",n=5,r=3,pmax=5,len=20,alpha=0.1){
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
c<-summary(mdl)$coefficients
form1<-frm
for (i in 1:(r*pmax))
  if (c[i,4]<=alpha)
    form1<-paste(form1,"+b[,",i,"]",sep="")
beta<-0
for (i in 1:pmax)
  if (c[r*pmax+i,4]<=alpha)
  {
    form1<-paste(form1,"+b[,",r*pmax+i,"]",sep="")
    beta<-beta+1
  }
if (beta==0)
  form1<-paste(form1,"+b[,",r*pmax+order(c[(r*pmax+1):((r+1)*pmax),4])[1],"]",sep="")
form1<-as.formula(form1)# formula for new modeling
mdl1<-lm(form1,data=b)# new linear regression model
index<-c(1:len)# set of labels of the time axis
ymin<-min(min(dec$appr[1,(day+1):dec$days]),min(fitted(mdl)),min(fitted(mdl1)))
ymax<-max(max(dec$appr[1,(day+1):dec$days]),max(fitted(mdl)),max(fitted(mdl1)))
plot(index,dec$appr[1,(day+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",pmax=",pmax,",a=",alpha,sep=""))# graphic of the time series
lines(index,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphic of the fitted values of the model
lines(index,fitted(mdl1),xlab="",ylab="",type='l',col="green")# graphic of the fitted values of the new model
return(mdl1)
}