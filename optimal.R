optimal<-function(lt){
d<-dim(lt)
ls<-vector(mode="character",length=d[1])
for (i in 1:d[1])
{
  ls[i]<-""
  for (j in 1:(d[2]-1))
    ls[i]<-paste(ls[i],lt[i,j],sep="")
}
lf<-as.factor(ls)
tf<-table(lf)
m<-max(tf)# mode quantity
nm<-order(tf)[length(tf)]
w<-levels(lf)[nm]# mode of models (i.e., p-vector in text format)
v<-as.numeric(w)
p<-vector(mode="numeric",length=d[2]-1)
for (i in 1:(d[2]-1))
{
  p[i]<-v%/%10^(d[2]-1-i)
  v<-v%%10^(d[2]-1-i)
}# we have mode of models (i.e., p-vector in numeric format)
nf<-as.integer(lf)
sm<-0
for (i in 1:d[1])
  if (nf[i]==nm)
    sm<-sm+lt[i,d[2]]
mn<-round(sm/m)# mean of len for the optimal model
return(list(mode=m,vect=p,t=w,m=mn))
}