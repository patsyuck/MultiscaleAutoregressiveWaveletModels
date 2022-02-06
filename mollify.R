mollify<-function(x,filter){
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
lx<-length(x)# length of the input time series
a<-vector(mode="numeric",length=lx)
for (t in l:lx)
{
a[t]<-0
for (k in 1:l)
  a[t]<-a[t]+h[k]*x[t-k+1]
}# we have the mollifying time series
return(a[l:lx])
}