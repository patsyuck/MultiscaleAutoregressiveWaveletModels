progn2m<-function(xy,z=NULL,filter="haar",n=5,r=5,pmax=5,lenmin=NULL,lenmax=NULL,alpha=0.1,fut=10){
lxy<-length(xy)
x<-xy[1:(lxy-fut)]
y<-xy[(lxy-fut+1):lxy]
lmax<-l2filter(x,n,r,pmax)
dec<-holesDWT(x,filter,n,r,pmax,len=NULL,lmax)# multiresolution decomposition
if (is.null(lenmin))
  lenmin<-dec$lenmin
else
  if (lenmin<dec$lenmin)
    stop("Parameter 'lenmin' is very small!")
if (is.null(lenmax))
  lenmax<-dec$lenmax
else
  if (lenmax>dec$lenmax)
    stop("Parameter 'lenmax' is very large!")
day<-dec$days-lenmax
bm<-r.matrix(dec,r,pmax,day,lenmax)# matrix for modeling
mt<-matrix(nrow=lenmax-lenmin+1,ncol=r+3)# box of matrix for the models test
frm<-paste("b[,",(r+1)*pmax+1,"]~0",sep="")
for (len in lenmin:lenmax)
{
  b<-as.data.frame(bm[(lenmax-len+1):lenmax,])
  form<-frm
  for (i in 1:((r+1)*pmax))
    form<-paste(form,"+b[,",i,"]",sep="")
  pold<-vector(mode="numeric",length=r+1)
  p<-rep(pmax,r+1)
  iter<-0
  lc<-left2cut(b,frm,form,r,pmax,p,pold,alpha,iter)# the optimal model
  mt[len-lenmin+1,1:(r+1)]<-lc$vect# the optimal vector for len
  mt[len-lenmin+1,r+2]<-len
  cf<-coef(lc$model)# coefficients of the model
  lz<-dec$days+fut
  pr<-future(dec,lc$vect,n,r,cf,fut,lz)# prognosis
  mt[len-lenmin+1,r+3]<-sum((pr-y[1:fut])^2)
}
nm<-order(mt[,r+3])[1]
p<-mt[nm,1:(r+1)]
nm<-mt[nm,r+2]
b<-as.data.frame(bm[(lenmax-nm+1):lenmax,])
form<-frm
for (i in 1:(r+1))
  if (p[i]>0)
    for (j in 1:p[i])
      form<-paste(form,"+b[,",(i-1)*pmax+j,"]",sep="")
form<-as.formula(form)# formula for optimal modeling
mdl<-lm(form,data=b)# optimal linear regression model
cf<-coef(mdl)# coefficients of the optimal model
lz<-dec$days+2*fut
pr<-future(dec,p,n,r,cf,2*fut,lz)# prognosis
train<-pr[1:fut]
test<-pr[(fut+1):(2*fut)]
dev1<-sd(train)
dev2<-sd(test)
w<-""
for (i in 1:(r+1))
  w<-paste(w,p[i],sep="")
par(mfrow=c(1,2))
ind<-c(1:(nm+fut))# set of labels of the time axis for left graphics
ind1<-c(1:nm)# set of labels of the time axis for regression data set
ind2<-c((nm+1):(nm+fut))# set of labels of the time axis for train set
ymin<-min(min(dec$appr[1,(dec$days-nm+1):dec$days]),min(fitted(mdl)),min(train)-2*dev1)
ymax<-max(max(dec$appr[1,(dec$days-nm+1):dec$days]),max(fitted(mdl)),max(train)+2*dev1)
plot(ind,xy[(lxy-fut-nm+1):lxy],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",p=",pmax,",a=",alpha,",p=",w,",l=",nm,sep=""))
lines(ind1,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphics of the time series and fitted values of the model
lines(ind2,train,xlab="",ylab="",type='l',col="blue")# graphic of the train prognosis
lines(ind2,train-2*dev1,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the bottom boundary of the train prognosis
lines(ind2,train+2*dev1,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the top boundary of the train prognosis
ind<-c(1:fut)# set of labels of the time axis for right graphics
if (is.null(z))
{
  ymin<-min(test)-2*dev2
  ymax<-max(test)+2*dev2
}
else
{
  ymin<-min(min(test)-2*dev2,min(z[1:fut]))
  ymax<-max(max(test)+2*dev2,max(z[1:fut]))
}
plot(ind,test,xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="blue",main="prognosis")# graphic of the prognosis
lines(ind,test-2*dev2,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the bottom boundary
lines(ind,test+2*dev2,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the top boundary
if (!(is.null(z)))
  lines(ind,z[1:fut],xlab="",ylab="",type='l',col="red")# graphic of the time series
return(test)
}