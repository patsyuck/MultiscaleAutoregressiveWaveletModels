model<-function(x,filter="haar",n=5,r=3,p=c(3,3,3,3),len=20){
pmax<-max(p)
lmax<-l2filter(x,n,r,pmax)
dec<-holesDWT(x,filter,n,r,pmax,len,lmax)# multiresolution decomposition
day<-dec$days-len
b<-r.matrix(dec,r,p,day,len)# matrix for modeling
form<-paste("b[,",sum(p)+1,"]~0",sep="")
for (i in 1:sum(p))
  form<-paste(form,"+b[,",i,"]",sep="")
form<-as.formula(form)# formula for modeling
mdl<-lm(form,data=b)# linear regression model
index<-c(1:len)# set of labels of the time axis
ymin<-min(min(dec$appr[1,(day+1):dec$days]),min(fitted(mdl)))
ymax<-max(max(dec$appr[1,(day+1):dec$days]),max(fitted(mdl)))
plot(index,dec$appr[1,(day+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red")# graphic of the time series
lines(index,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphic of the fitted values of the model
return(mdl)
}