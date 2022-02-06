model4<-function(x,filter="haar",n=5,r=3,pmax=3,len=20){
lmax<-l2filter(x,n,r,pmax)
dec<-holesDWT(x,filter,n,r,pmax,len,lmax)# multiresolution decomposition
day<-dec$days-len
b<-r.matrix(dec,r,pmax,day,len)# matrix for modeling
frm<-paste("b[,",(r+1)*pmax+1,"]~0",sep="")
form<-frm
for (j in 1:((r+1)*pmax))
  form<-paste(form,"+b[,",j,"]",sep="")
form<-as.formula(form)# formula for modeling
mdl<-lm(form,data=b)# linear regression model
lx<-length(x)
RSS<-sum((x[(lx-len+1):lx]-fitted(mdl))^2)
t<-vector(mode="numeric",length=(r+1)*pmax)
for (i in 1:((r+1)*pmax))
{
  q<-rep(1,(r+1)*pmax)
  q[i]<-0
  form<-frm
  for (j in 1:((r+1)*pmax))
    if (q[j]==1)
      form<-paste(form,"+b[,",j,"]",sep="")
  form<-as.formula(form)# formula for modeling
  mdl<-lm(form,data=b)# linear regression model
  rss<-sum((x[(lx-len+1):lx]-fitted(mdl))^2)
  t[i]<-sqrt((len-(r+1)*pmax)*abs(rss/RSS-1))
}
ordt<-order(t,decreasing=TRUE)
Cp<-vector(mode="numeric",length=(r+1)*pmax)
q<-vector(mode="numeric",length=(r+1)*pmax)
for (i in 1:((r+1)*pmax))
{
  q[ordt[i]]<-1
  form<-frm
  for (j in 1:((r+1)*pmax))
    if (q[j]==1)
      form<-paste(form,"+b[,",j,"]",sep="")
  form<-as.formula(form)# formula for modeling
  mdl<-lm(form,data=b)# linear regression model
  rss<-sum((x[(lx-len+1):lx]-fitted(mdl))^2)
  Cp[i]<-(len-(r+1)*pmax)*rss/RSS+2*i-len
}
minCp<-order(Cp)[1]
q<-vector(mode="numeric",length=(r+1)*pmax)
for (i in 1:minCp)
  q[ordt[i]]<-1
form<-frm
  for (j in 1:((r+1)*pmax))
    if (q[j]==1)
      form<-paste(form,"+b[,",j,"]",sep="")
form<-as.formula(form)# formula for modeling
mdl<-lm(form,data=b)# linear regression model
w<-""
for (i in 1:((r+1)*pmax))
  w<-paste(w,q[i],sep="")
index<-c(1:len)# set of labels of the time axis
ymin<-min(min(dec$appr[1,(day+1):dec$days]),min(fitted(mdl)))
ymax<-max(max(dec$appr[1,(day+1):dec$days]),max(fitted(mdl)))
plot(index,dec$appr[1,(day+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",p=",pmax,",q=",w,",l=",len,sep=""))# graphic of the time series
lines(index,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphic of the fitted values of the model
return(list(q,mdl))
}