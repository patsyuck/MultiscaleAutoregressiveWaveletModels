holesDWT<-function(x,filter,n,r,pmax,len=NULL,lmax=NULL){
library(wmtsa)# Wavelet Methods for Time Series Analysis
switch(filter,
  "haar"={h<-wavDaubechies("d2")$scaling},
  "db2"={h<-wavDaubechies("d4")$scaling},
  "db3"={h<-wavDaubechies("d6")$scaling},
  "db4"={h<-wavDaubechies("d8")$scaling},
  "db5"={h<-wavDaubechies("d10")$scaling},
  "db6"={h<-wavDaubechies("d12")$scaling},
  "db7"={h<-wavDaubechies("d14")$scaling},
  "db8"={h<-wavDaubechies("d16")$scaling},
  "db9"={h<-wavDaubechies("d18")$scaling},
  "db10"={h<-wavDaubechies("d20")$scaling},
  "sym4"={h<-wavDaubechies("s8")$scaling},
  "sym5"={h<-wavDaubechies("s10")$scaling},
  "sym6"={h<-wavDaubechies("s12")$scaling},
  "sym7"={h<-wavDaubechies("s14")$scaling},
  "sym8"={h<-wavDaubechies("s16")$scaling},
  "sym9"={h<-wavDaubechies("s18")$scaling},
  "sym10"={h<-wavDaubechies("s20")$scaling},
  "coif1"={h<-wavDaubechies("c6")$scaling},
  "coif2"={h<-wavDaubechies("c12")$scaling},
  "coif3"={h<-wavDaubechies("c18")$scaling},
  "coif4"={h<-wavDaubechies("c24")$scaling},
  "coif5"={h<-wavDaubechies("c30")$scaling},
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