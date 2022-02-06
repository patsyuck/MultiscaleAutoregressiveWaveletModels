holesDWT<-function(x,filter,n,r,pmax,len=NULL,lmax=NULL){
switch(filter,
  "haar"={h<-c(0.5,0.5)},
  "db2"={h<-c(0.3415,0.5915,0.1585,-0.0915)},
  "db3"={h<-c(0.2352,0.5706,0.3252,-0.0955,-0.0604,0.0249)},
  "db4"={h<-c(0.1629,0.5055,0.4461,-0.0198,-0.1323,0.0218,0.0233,-0.0075)},
  "sym4"={h<-c(0.0228,-0.0089,-0.0702,0.2106,0.5683,0.3519,-0.0210,-0.0536)},
  "coif1"={h<-c(-0.0514,0.2389,0.6029,0.2721,-0.0514,-0.0111)},
  )# what kind of filter we shall use
l<-length(h)# length of the using filter
if (!(is.null(lmax)))
  if (l>lmax)
    stop("Unsufficient data! Use a filter of the less length.")
lx<-length(x)# length of the input time series
lenmin<-len.min(r,pmax)
lenmax<-len.max(lx,l,n,r,pmax)
if (!(is.null(len)))
{
  if (len<lenmin)
    stop("Parameter 'len' is very small!")
  if (len>lenmax)
    stop("Parameter 'len' is very large!")
}
a<-matrix(nrow=n+1,ncol=lx)# box of the approximation matrix 
d<-matrix(nrow=n,ncol=lx)# box of the detalization matrix
for (t in 1:lx)
  a[1,t]<-x[t]
lenA<-lx
for (j in 1:n)
{
  lenA<-lenA-(l-1)*2^(j-1)
  for (t in (lx-lenA+1):lx)
  {
    a[j+1,t]<-0
    for (k in 1:l)
      a[j+1,t]<-a[j+1,t]+h[k]*a[j,t-(k-1)*2^(j-1)]
    d[j,t]<-a[j,t]-a[j+1,t]
  }
}# we have the approximation and detalization matrices
list(f=h,lf=l,lenmin=lenmin,lenmax=lenmax,days=lenA,appr=a[,(lx-lenA+1):lx],detal=d[,(lx-lenA+1):lx])
}