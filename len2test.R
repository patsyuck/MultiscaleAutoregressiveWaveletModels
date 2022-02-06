len2test<-function(x,filter="haar",n=5,r=5,pmax=5,lenmin=NULL,lenmax=NULL,alpha=0.1,matr=FALSE){
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
lt<-matrix(nrow=lenmax-lenmin+1,ncol=r+2)# box of matrix for len test
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
  lt[len-lenmin+1,1:(r+1)]<-lc$vect# the optimal vector for len
  lt[len-lenmin+1,r+2]<-len
}
opt<-optimal(lt)# mode of the optimal models
b<-as.data.frame(bm[(lenmax-opt$m+1):lenmax,])
form<-frm
for (i in 1:(r+1))
  if (opt$vect[i]>0)
    for (j in 1:opt$vect[i])
      form<-paste(form,"+b[,",(i-1)*pmax+j,"]",sep="")
form<-as.formula(form)# formula for optimal modeling
mdl<-lm(form,data=b)# optimal linear regression model
ind<-c(1:opt$m)# set of labels of the time axis
ymin<-min(min(dec$appr[1,(dec$days-opt$m+1):dec$days]),min(fitted(mdl)))
ymax<-max(max(dec$appr[1,(dec$days-opt$m+1):dec$days]),max(fitted(mdl)))
plot(ind,dec$appr[1,(dec$days-opt$m+1):dec$days],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",p=",pmax,",a=",alpha,",p=",opt$t,",l=",opt$m,sep=""))
lines(ind,fitted(mdl),xlab="",ylab="",type='l',col="blue")# graphics of the time series and fitted values of the model
if (matr==FALSE)
  return(list(mode.quant=c(opt$mode,lenmax-lenmin+1),mode.share=opt$mode/(lenmax-lenmin+1),mean.day=opt$m,vect=opt$vect,model=mdl))
else
  return(list(matr=lt,mode.quant=c(opt$mode,lenmax-lenmin+1),mode.share=opt$mode/(lenmax-lenmin+1),mean.day=opt$m,vect=opt$vect,model=mdl))
}