gprint<-function(d,leng,l,TS,TSmol,TSappr,fut,TSfut,TSprogn,filter,r,pmax,alpha,p){
par(mfrow=c(1,2))
ind<-c(1:l)# set of labels of the time axis for left graphics
ymin<-min(min(TS[(leng-l+1):leng]),min(TSappr))
ymax<-max(max(TS[(leng-l+1):leng]),max(TSappr))
plot(ind,TS[(leng-l+1):leng],xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="red",main=paste(filter,",r=",r,",pmax=",pmax,",a=",alpha,",p=",p,",l=",l,sep=""))
grid(nx=NULL,ny=NULL)
if (!(is.null(d)))
  lines(ind,TSmol[(d-l+1):d],xlab="",ylab="",type='l',col="green")
lines(ind,TSappr,xlab="",ylab="",type='l',col="blue")# graphics of the time series and fitted values of the model
ind<-c(1:fut)# set of labels of the time axis for right graphics
dev<-sd(TSprogn)
if (is.null(TSfut))
{
  ymin<-min(TSprogn)-2*dev
  ymax<-max(TSprogn)+2*dev
}
else
{
  ymin<-min(min(TSprogn)-2*dev,min(TSfut[1:fut]))
  ymax<-max(max(TSprogn)+2*dev,max(TSfut[1:fut]))
}
plot(ind,TSprogn,xlab="",ylab="",ylim=c(ymin,ymax),type='l',col="blue",main="prognosis")# graphic of the prognosis
grid(nx=NULL,ny=NULL)
lines(ind,TSprogn-2*dev,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the bottom boundary
lines(ind,TSprogn+2*dev,xlab="",ylab="",type='l',lty=2,col="blue")# graphic of the top boundary
if (!(is.null(TSfut)))
  lines(ind,TSfut[1:fut],xlab="",ylab="",type='l',col="red")# graphic of the time series
return()
}