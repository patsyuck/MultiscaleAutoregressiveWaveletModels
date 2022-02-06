progn5<-function(xz,y=NULL,fmol="haar",filter="haar",n=5,r=5,pmax=5,lenmin=NULL,lenmax=NULL,alpha=0.1,fut=10){
leng<-length(xz)
x<-mollify(xz,fmol)
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
cf<-coef(mdl)# coefficients of the optimal model
lz<-dec$days+fut
pr<-future(dec,opt$vect,n,r,cf,fut,lz)# prognosis
gprint(dec$days,leng,opt$m,xz,dec$appr[1,],fitted(mdl),fut,y,pr,filter,r,pmax,alpha,opt$t)
return(pr)
}