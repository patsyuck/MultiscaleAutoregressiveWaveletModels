model1b<-function(x,filter="haar",n=5,r=3,pmax=5,len=20,alpha=0.1){
lmax<-l2filter(x,n,r,pmax)
dec<-holesDWT(x,filter,n,r,pmax,len,lmax)# multiresolution decomposition
day<-dec$days-len
b<-r.matrix(dec,r,pmax,day,len)# matrix for modeling
frm<-paste("b[,",(r+1)*pmax+1,"]~0",sep="")
form<-frm
for (i in 1:((r+1)*pmax))
  form<-paste(form,"+b[,",i,"]",sep="")
pold<-vector(mode="numeric",length=r+1)
p<-rep(pmax,r+1)
iter<-0
lc<-left.cut(b,frm,form,r,pmax,p,pold,alpha,iter)# the optimal model
w<-""
for (i in 1:(r+1))
  w<-paste(w,lc$vect[i],sep="")
index<-c(1:len)# set of labels of the time axis
ymin<-min(min(dec$appr[1,(day+1):dec$days]),min(fitted(lc$model)))
ymax<-max(max(dec$appr[1,(day+1):dec$days]),max(fitted(lc$model)))
plot(index,dec$appr[1,(day+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",pmax=",pmax,",a=",alpha,",p=",w,",iter=",lc$iter,sep=""))
lines(index,fitted(lc$model),xlab="",ylab="",type='l',col="blue")# graphics of the time series and fitted values of the model
return(lc$model)
}