m.form<-function(p,r=NULL){
if (length(p)==1)
  s<-(r+1)*p
else
  s<-sum(p)
form<-paste("b[,",s+1,"]~0",sep="")
for (i in 1:s)
  form<-paste(form,"+b[,",i,"]",sep="")
return(as.formula(form))# formula for modeling
}