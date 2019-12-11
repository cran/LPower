LPower=function(rx_effect=NULL,sample_size=NULL,allocationRatio=c(1,1),power=NULL,
                contrast=c(rep(0,dim(xMatrix[[1]])[2]-1),1),xMatrix,
                 vMatrix,attritionRates=0,alpha=0.025,simulate=FALSE,
                 nsims=1000,betas=c(rep(0,dim(xMatrix[[1]])[2]-1),1))
#xMatrix is a list of xMatrixRx0, xMatrixRx1
#att_rates vector of probabilities of attending ith visit 1st=1 last=remaining patients
  {
  drops=function(x) 1-exp(-x)
  m=dim(xMatrix[[1]])
  vsv=matrix(0,m[2],m[2])
  if(class(vMatrix)[1]=='matrix')
                     vMatrix=replicate(length(xMatrix),vMatrix,simplify=FALSE)
  if(length(attritionRates)==1) attritionRates=rep(attritionRates,m[1]-1)
  for (i in 1:length(xMatrix)){
    remains=exp(-attritionRates[1])
    if(class(vMatrix)[1]=='list') vs=vMatrix[[i]] else vs=vMatrix
    for (j in 2:m[1]){
      if(j<m[1]){
      pp=remains*drops(attritionRates[j])
      remains=remains-pp} else pp=remains
      if(pp>0){
      vis=1:j
      x=xMatrix[[i]][vis,]
      vs=vMatrix[[i]][vis,vis]
      vsv=vsv+(t(x)%*%MASS::ginv(vs)%*%x)*pp*allocationRatio[i]
      }
    }
  }
  vsvt=MASS::ginv(vsv)
  #sdb=sqrt(vsvt[[3,3]])
  a=matrix(contrast,1,m[2])
  sdb=sqrt(a%*%vsvt%*%t(a))
  mult=sample_size/sum(allocationRatio);
  if(is.null(power))
  {
    power=pnorm(sqrt(mult)*rx_effect/sdb-qnorm(1-alpha));
  }else if(is.null(rx_effect))
  {
    rx_effect=(qnorm(power)+qnorm(1-alpha))*sdb/sqrt(mult);
  }else
  {
    sample_size=ceiling(sum(allocationRatio)*((qnorm(power)+qnorm(1-alpha))*sdb/rx_effect)^2);
  }

  if(simulate){
  sRatio=sum(allocationRatio)
  mcontrast=matrix(contrast)
  betas=rep(rx_effect/(matrix(betas,1,m[2])%*%mcontrast),m[2])*betas
  criteria=qnorm(1-alpha)
  rej=0
  for(i in 1:nsims){

  df=data.frame()
  for(i in 1:length(xMatrix)){
    nSubjects=floor(allocationRatio[i]*sample_size/sum(allocationRatio))
    y=MASS::mvrnorm(nSubjects,mu=xMatrix[[i]]%*%matrix(betas),vMatrix[[i]])
    dropTimes=approx(cumsum(c(0,attritionRates)),0:(m[1]-1),
                     xout=-log(runif(nSubjects,max=m[1]-1)),rule=2)$y
  start=ifelse(is.null(df$subject),1,max(df$subject)+1)
  df1=data.frame(subject=rep(start:(start+nSubjects-1),each=m[1]),
                 dropTimes=rep(dropTimes,each=m[1]),visits=rep(1:m[1],nSubjects),
              y=c(t(y)),t(matrix(rep(t(xMatrix[[i]]),nSubjects),m[2],nSubjects))
            )
  df=rbind(df,df1[df1$dropTimes>=(df1$visits-1),])
  }

  corForm=as.formula(paste("~","visits","|","subject"))
  varForm=as.formula(paste("~1|","visits"))
  xvars=paste(colnames(df)[5:dim(df)[2]],collapse='+')
  meanForm=as.formula(paste("y~-1+",xvars))
  mod2=nlme::gls(meanForm,
                 correlation=corSymm(form=corForm),weights=varIdent(form=varForm),
                 na.action=na.omit,data=df)
  rej=rej+(sum(contrast*mod2$coefficients)/sqrt( t(mcontrast)%*%mod2$varBeta%*%mcontrast)>criteria)
  }
  spower=rej/nsims
} else spower=NULL

  return(c(rx_effect=rx_effect,sample_size=sample_size,power=power,simulatedPower=spower));
}

randomEffectsMatrix=function(zMatrix,vs,sigma2){
  if(class(zMatrix)[1]=='list'){
  vMatrix=list()
  for (i in 1:length(zMatrix))
  vMatrix[[i]]=sigma2*diag(rep(1,dim(zMatrix[[i]])[1]))+zMatrix[[i]]%*%vs%*%t(zMatrix[[i]])
  }
 else vMatrix=sigma2*diag(rep(1,dim(zMatrix)[1]))+zMatrix%*%vs%*%t(zMatrix)
 return(vMatrix)
 }

randomSlopesMatrix= function(visit,vs,sigma2,dropPerMonth,baselineTreatment=FALSE)
  {
  if(class(visit)[1]!='list') visit=list(visit,visit)
  m=length(visit[[1]])
  xMatrix=z=vMatrix=list()
  for(i in 1:2){
  rxMatrix=matrix(c(1,i-1),1,2)
  visitMatrix=matrix(visit[[i]])
  if(baselineTreatment)  xMatrix[[i]]=cbind(1,matrix(1,m[1])%*%(i-1),visitMatrix%*%rxMatrix)
  else xMatrix[[i]]=cbind(1,visitMatrix%*%rxMatrix)
  z[[i]]=cbind(rep(1,m),matrix(visit[[i]]))}
  vMatrix=randomEffectsMatrix(z,vs,sigma2)
  #vMatrix[[i]]=sigma2*diag(rep(1,m))+z[[i]]%*%vs%*%t(z[[i]])
  attritionRates=dropPerMonth*(visit[[i]][2:m]-visit[[i]][1:(m-1)])
  return(list(xMatrix=xMatrix,vMatrix=vMatrix,attritionRates=attritionRates))
}



