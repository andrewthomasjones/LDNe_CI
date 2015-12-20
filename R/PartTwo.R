library(adegenet)
library(genetics)
library(tcltk)
library(reshape)
library(inline)
library(rbenchmark)
library(compiler)
library(ggplot2)


#will need to set working dir to this file
dire<-paste(dirname(getwd()),"/Results/", sep="") # or similar

#functions for  various calcs
NehatW<-function(r,S){
  r2prime<-(r)-((1/S)+(3.19/S^2))
  return(((1/3)+sqrt((1/9)-2.76*r2prime))/(2*r2prime))
}

NehatSved<-function(r,SvCorr){
  return(1/(3*(r-SvCorr)))
}

parametric <- function(r,S,J,alpha){
  lower<- NehatWC(r*J/qchisq(0.5*alpha, J),S)
  upper<-NehatWC(r*J/qchisq(1-(0.5*alpha), J),S)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

Chijackknife <- function(rT,S,J1,alpha){
  upper<-NehatW(rT*J1 / qchisq(0.5*alpha, J1),S)
  lower<-NehatW(rT*J1 / qchisq((1-alpha/2), J1),S)
  ConfInt<-cbind(lower,upper)
  return(ConfInt)
}

ChijackknifeSv <- function(rT,S,J1,alpha,Sv){
  
  upper<-NehatSved(rT*J1 / qchisq(0.5*alpha, J1),Sv)
  lower<-NehatSved(rT*J1 / qchisq((1-alpha/2), J1),Sv)
  ConfInt<-cbind(lower,upper)
  return(ConfInt)
}

Normjackknife <- function(rT,S,Se,alpha,jkmean){
  
  mult<-qnorm(alpha/2)  
  upper<-NehatW(jkmean+mult*Se,S)
  lower<-NehatW(jkmean-mult*Se,S)
  ConfInt<-cbind(lower,upper)
  return(ConfInt)
}

NormjackknifeSv <- function(rT,S,Se,alpha,Sv,jkmean){
  
  mult<-qnorm(alpha/2)  
  upper<-NehatSved(jkmean+mult*Se,Sv)
  lower<-NehatSved(jkmean-mult*Se,Sv)
  ConfInt<-cbind(lower,upper)
  return(ConfInt)
}   


files<-dir(dire)

csvFiles<-paste0(dire,files[grep('*.csv', files)])

csvList<-list()


for(i in 1:length(csvFiles)){
  temp1<-read.csv(csvFiles[i])
  temp1$run<-rep(unlist(strsplit(csvFiles[i],".csv")),nrow(temp1))
  csvList[[length(csvList)+1]]<-temp1
}

complete.dat <- do.call(rbind, csvList)
complete.dat$refNe<- 1/(( (1/complete.dat$trueNe)+ 0.5*(1/complete.dat$trueNe2) + 0.25*(1/complete.dat$trueNe3) + 0.125*(1/complete.dat$trueNe4))/1.875)

complete.dat$alpha2<-as.numeric(as.character((complete.dat$alpha)))
complete.dat$alpha<-factor(complete.dat$alpha)
complete.dat$pop_size<-factor(complete.dat$pop_size)
complete.dat$sample_size2<-as.numeric(as.character(((complete.dat$sample_size))))
complete.dat$sample_size<-factor(complete.dat$sample_size)
complete.dat$n_loci<-factor(complete.dat$n_loci)
complete.dat$al_per_loc<-factor(complete.dat$al_per_loc)
complete.dat$run<-factor(basename(as.character(complete.dat$run)))
complete.dat$pcrit<-factor(complete.dat$pcrit)
#adjust variance here
fud<-(0.84)^2
complete.dat$Se <- fud*sqrt(complete.dat$Va)
complete.dat$Js <- 2/ ((complete.dat$Se)^2 / (complete.dat$r2^2))


#redo CI based on new J*
 complete.dat$Cjup <- Chijackknife(complete.dat$r2, complete.dat$sample_size2, complete.dat$Js,  complete.dat$alpha2)[,1]
 complete.dat$Cjdown <-Chijackknife(complete.dat$r2,complete.dat$sample_size2, complete.dat$Js,  complete.dat$alpha2)[,2]
 
 complete.dat$Njup <-Normjackknife(complete.dat$r2,complete.dat$sample_size2, complete.dat$Se , complete.dat$alpha2, complete.dat$r2)[,2]
 complete.dat$Njdown <- Normjackknife(complete.dat$r2,complete.dat$sample_size2, complete.dat$Se , complete.dat$alpha2, complete.dat$r2)[,1] 
 
 complete.dat$CjupSv <- ChijackknifeSv(complete.dat$r2,complete.dat$sample_size2, complete.dat$Js,  complete.dat$alpha2,complete.dat$Sved)[,1]
 complete.dat$CjdownSv <-ChijackknifeSv(complete.dat$r2,complete.dat$sample_size2, complete.dat$Js,  complete.dat$alpha2,complete.dat$Sved)[,2]
 
 complete.dat$NjupSv <-NormjackknifeSv(complete.dat$r2,complete.dat$sample_size2, complete.dat$Se , complete.dat$alpha2, complete.dat$Sved, complete.dat$r2)[,2]
 complete.dat$NjdownSv <- NormjackknifeSv(complete.dat$r2,complete.dat$sample_size2, complete.dat$Se , complete.dat$alpha2, complete.dat$Sved, complete.dat$r2)[,1] 



#set all neg estimates to be VERY LARGE
complete.dat[(complete.dat<0)]<-10^9

splt.by <- c('run','pcrit','alpha')
splitData<-split(complete.dat, complete.dat[,splt.by] )






CICalc<-function(data){
  
  N<-dim(data)[1]
  
  Cj<-sum(data$trueNe>data$Cjdown & data$trueNe<data$Cjup)/N
  Nj<-sum(data$trueNe>data$Njdown & data$trueNe<data$Njup)/N
  
  lj<-sum(data$trueNe>data$Ljdown & data$trueNe<data$lJup)/N
  p<-sum(data$trueNe>data$pdown & data$trueNe<data$pup)/N
  
  Nj2<-sum(data$trueNe>data$Nj2down & data$trueNe<data$Nj2up)/N
  CjSv<-sum(data$trueNe>data$CjdownSv & data$trueNe<data$CjupSv)/N
  NjSv<-sum(data$trueNe>data$NjdownSv & data$trueNe<data$NjupSv)/N
  ljSv<-sum(data$trueNe>data$LjdownSv & data$trueNe<data$lJupSv)/N
  pSv<-sum(data$trueNe>data$pdownSv & data$trueNe<data$pupSv)/N
  
  alpha<-as.numeric(as.character(data$alpha[1]))
  run<-as.character(data$run[1])
  pcrit<-as.numeric(as.character(data$pcrit[1]))
  pop_size<-as.numeric(as.character(data$pop_size[1]))
  sample_size<-as.numeric(as.character(data$sample_size[1]))
  n_loci<-as.numeric(as.character(data$n_loci[1]))
  al_per_loc<-as.numeric(as.character(data$al_per_loc[1]))
  
  data$NeSv2<-data$NeSv
  data$Ne2<-data$Ne
  
  data$Ne2[data$Ne2<0]<-NA
  data$NeSv2[data$NeSv2<0]<-NA
  
  biasSv<-mean(data$trueNe-data$NeSv2, na.rm=T)
  bias<-mean(data$trueNe-data$Ne2, na.rm=T)
  varSv<-var(data$trueNe-data$NeSv2, na.rm=T)
  var<-var(data$trueNe-data$Ne2, na.rm=T)
  
  J<-mean(as.numeric(as.character(data$IndComp)))
  Jprime<-mean(as.numeric(as.character(data$Jp)),na.rm=T)
  Jstar<-mean(as.numeric(as.character(data$Js)),na.rm=T)
  
  
  CjA<-sum(data$trueNe>data$Cjdown )/N
  NjA<-sum(data$trueNe>data$Njdown )/N
  
  ljA<-sum(data$trueNe>data$Ljdown )/N
  pA<-sum(data$trueNe>data$pdown )/N
  
  
  CjB<-sum( data$trueNe<data$Cjup)/N
  NjB<-sum( data$trueNe<data$Njup)/N
  
  ljB<-sum(data$trueNe<data$lJup)/N
  pB<-sum(data$trueNe<data$pup)/N
  
  return(list(Cj=Cj,Nj=Nj, lj=lj,p=p, CjSv=CjSv,NjSv=NjSv,ljSv=ljSv,pSv=pSv, alpha=alpha, pcrit=pcrit, J=J, Jprime=Jprime, Jstar = Jstar, run=run,pop_size=pop_size,sample_size=sample_size,n_loci=n_loci,al_per_loc=al_per_loc,  bias= bias, biasSv=biasSv, var=var, varSv=varSv,CjB=CjB,NjB=NjB, ljB=ljB,pB=pB,CjA=CjA,NjA=NjA, ljA=ljA,pA=pA ))
}


#actually combine it all
outList<-list()
for(i in 1:length(splitData)){
  outList[[length(outList)+1]] <- CICalc(splitData[[i]]) 
}

#produce final data frame
final <- data.frame(matrix(unlist(outList), nrow=length(splitData), byrow=T))
names(final)<-names(CICalc(splitData[[1]]))
levels(final$run)<-levels(splitData[[1]]$run)
write.csv(final,paste(dirname(dire), "/CombinedResults", ".csv", sep=""))
