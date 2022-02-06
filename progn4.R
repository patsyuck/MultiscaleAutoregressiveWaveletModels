progn4<-function(x,y=NULL,filter="haar",n=5,r=5,pmax=5,lenmin=NULL,lenmax=NULL,fut=10){
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
lt<-matrix(nrow=lenmax-lenmin+1,ncol=(r+1)*pmax+1)# box of matrix for len test
frm<-paste("b[,",(r+1)*pmax+1,"]~0",sep="")
lx<-length(x)
for (len in lenmin:lenmax)
{
  b<-as.data.frame(bm[(lenmax-len+1):lenmax,])
  form<-frm
  for (j in 1:((r+1)*pmax))
    form<-paste(form,"+b[,",j,"]",sep="")
  form<-as.formula(form)# formula for modeling
  mdl<-lm(form,data=b)# linear regression model
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
  lt[len-lenmin+1,1:((r+1)*pmax)]<-q# the optimal vector for len
  lt[len-lenmin+1,r+2]<-len
}
opt<-optimal(lt)# mode of the optimal models
b<-as.data.frame(bm[(lenmax-opt$m+1):lenmax,])
form<-frm
for (i in 1:((r+1)*pmax))
  if (opt$vect[i]==1)
    form<-paste(form,"+b[,",i,"]",sep="")
form<-as.formula(form)# formula for optimal modeling
mdl<-lm(form,data=b)# optimal linear regression model
cf<-coef(mdl)# coefficients of the optimal model
lz<-dec$days+fut
pr<-future(dec,opt$vect,n,r,cf,fut,lz)# prognosis
dev<-sd(pr)
par(mfrow=c(1,2))
ind<-c(1:opt$m)# set of labels of the time axis for left graphics
ymin<-min(min(dec$appr[1,(dec$days-opt$m+1):dec$days]),min(fitted(mdl)))
ymax<-max(max(dec$appr[1,(dec$days-opt$m+1):dec$days]),max(fitted(mdl)))
plot(ind,dec$appr[1,(dec$days-opt$m+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",p=",pmax,",p=",opt$t,",l=",opt$m,sep=""))
grid(nx=NULL,ny=NULL)
lines(ind,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphics of the time series and fitted values of the model
ind<-c(1:fut)# set of labels of the time axis for right graphics
if (is.null(y))
{
  ymin<-min(pr)-2*dev
  ymax<-max(pr)+2*dev
}
else
{
  ymin<-min(min(pr)-2*dev,min(y[1:fut]))
  ymax<-max(max(pr)+2*dev,max(y[1:fut]))
}
plot(ind,pr,xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="blue",main="prognosis")# graphic of the prognosis
grid(nx=NULL,ny=NULL)
lines(ind,pr-2*dev,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the bottom boundary
lines(ind,pr+2*dev,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the top boundary
if (!(is.null(y)))
  lines(ind,y[1:fut],xlab="",ylab="",type='l',col="red")# graphic of the time series
return(pr)
}