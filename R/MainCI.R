library(adegenet)
library(genetics)
library(tcltk)
library(Rcpp)
library(RcppEigen)
library(inline)
library(rbenchmark)
library(compiler)
library(gsubfn)
library(tools)
source("http://gsubfn.googlecode.com/svn/trunk/R/list.R")
#for C compiler
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

#will need to set working dir to this file
Rcpp::sourceCpp(paste(getwd(), '/CIcodeDec.cpp', sep=""))


##Set in and out directories########
#########################################
direIn<-paste(dirname(getwd()),"/Samples/", sep="") # or similar
direOut<-paste(dirname(getwd()),"/Results/", sep="")  # or similar

dir.create(direOut, showWarnings = FALSE) # creates output dir if doesnt already exist
################################################
#number of perms for Sved method 
M=10
#Confidence level for CI
alpha<-c(0.01,0.05,0.10)
# pcrits to use
pcrits<-c(0,0.01,0.02,0.05,0.1)
#############################################


#functions for  various calcs
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case,include.dirs=F)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}


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
parametricSv <- function(r,S,J,alpha,Sv){
  lower<- NehatSved(r*J/qchisq(0.5*alpha, J),Sv)
  upper<-NehatSved(r*J/qchisq(1-(0.5*alpha), J),Sv)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

#coeff. variation
CV<-function(x) 2/(var(x)/(mean(x)^2))
 

#'jackknife'
nprime<-function(rvals){
  len<-length(rvals)
  cvari<-array(0,len)
  for(i in 1:len){                                                                                                                                                                                                             
    cvari[i]<-CV(rvals[-i])
  }
  return(sum(cvari))
}

#'jackknife'
nprime2<-function(rvals){
  len<-length(rvals)
  cvari<-array(0,len)
  for(i in 1:len){
    cvari[i]<-CV(rvals[-i])
  }
  return(sum(cvari))
}


LDjackknife <- function(r,S,nprime,alpha){
  lower<- NehatWC(r*nprime/qchisq(0.5*alpha, nprime),S)
  upper<-NehatWC(r*nprime/qchisq(1-(0.5*alpha), nprime),S)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

LDjackknifeSv <- function(r,S,nprime,alpha,Sv){
  lower<- NehatSved(r*nprime/qchisq(0.5*alpha, nprime),Sv)
  upper<-NehatSved(r*nprime/qchisq(1-(0.5*alpha), nprime),Sv)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}


Chijackknife <- function(rT,S,Va,alpha){
  scalechisqa<-Va/(2*rT)
  scalechisqb<-rT/scalechisqa
  
  upper<-NehatWC(scalechisqa*qchisq(alpha/2, scalechisqb),S)
  lower<-NehatWC(scalechisqa*qchisq((1-alpha/2), scalechisqb),S)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

ChijackknifeSv <- function(rT,S,Va,alpha,Sv){
  scalechisqa<-Va/(2*rT)
  scalechisqb<-rT/scalechisqa
  
  upper<-NehatSved(scalechisqa*qchisq(alpha/2, scalechisqb),Sv)
  lower<-NehatSved(scalechisqa*qchisq((1-alpha/2), scalechisqb),Sv)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

Normjackknife <- function(rT,S,Se,alpha){
  
  mult<-qnorm(alpha/2)  
  upper<-NehatWC(rT+mult*Se,S)
  lower<-NehatWC(rT-mult*Se,S)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

NormjackknifeSv <- function(rT,S,Se,alpha,Sv){
  
  mult<-qnorm(alpha/2)  
  upper<-NehatSved(rT+mult*Se,Sv)
  lower<-NehatSved(rT-mult*Se,Sv)
  ConfInt<-c(lower,upper)
  return(ConfInt)
}

#standard error and varaince jacknife
jackfun<-function(data3,pcrit,M){
  indexing=data3@loc.nall
  N=dim(data3@tab)[1]
  Nj=N-1
  L= length(data3@loc.names)
  true=rBurrrowsS(data3@tab, N, pcrit, L, indexing,M)
  r=true[1]
  J=true[2]
  sved=true[3]
  rvals=true[-(1:3)]
  
  JKvalues= array(0,N)
  
  for(i in 1:N){
    data4<-data3[-i]
    data5<- data4@tab[,colSums(data4@tab)!=0]
    reIndex<- data4@loc.nall-unlist(lapply(split(colSums(data4@tab)==0, data4@loc.fac),sum))
    temp=rBurrrowsS( data5, Nj, pcrit, L,reIndex,0)
    JKvalues[i]=temp[1]
  }
  
  JKVar<- ((N-1)/N)*sum((JKvalues-mean(JKvalues,na.rm = T))^2,na.rm = T) 
  JKSe<- (((N-1)/N)*sum((JKvalues-mean(JKvalues,na.rm = T))^2,na.rm = T))^0.5
  
  
  return(c(r,J,sved, JKSe, JKVar, rvals))
}


#CI util fun

#returns true if lower, upper, true, fit
CIcalc<-function(l,u,N){
  return(((l<N) & (u>N)))
  
}
#compiled versions of some of the above

NehatWC<-cmpfun(NehatW)
jackfunC<-cmpfun(jackfun)
ChijackknifeC<-cmpfun(Chijackknife)
NormjackknifeC<-cmpfun(Normjackknife)
LDjackknifeC<-cmpfun(LDjackknife)
parametricC<-cmpfun(parametric)
ChijackknifeSvC<-cmpfun(ChijackknifeSv)
NormjackknifeSvC<-cmpfun(NormjackknifeSv) 
LDjackknifeSvC<-cmpfun(LDjackknifeSv)
parametricSvC<-cmpfun(parametricSv)


######################################################
#actual runs
#
#####################################################

dirlist<-list.dirs(direIn, full.names=T)
outputNames<-list.dirs(direIn, full.names=F)

################## code here so it doesnt redo when i add more folders of sims

currentFiles<-file_path_sans_ext(list.files(direOut))
dirlist<-dirlist[!outputNames%in%currentFiles]
outputNames<-outputNames[!outputNames%in%currentFiles]

####################
len_dir<-(length(dirlist))
for(y in 1:length(dirlist)){
  df<-data.frame()
  print(dirlist[y])
  filelist<-list.files(dirlist[y])
  write.csv(df,paste(direOut, outputNames[y],".csv", sep=""))
  dirlist<-list.dirs(direIn, full.names=T)
  
  currentFiles<-file_path_sans_ext(list.files())
  dirlist<-dirlist[!outputNames%in%currentFiles]
  
  

  
  
  
  parameters<-read.csv(paste(dirlist[y],"/", filelist[1], sep=""))
  
  runInfo<-read.csv(paste(dirlist[y],"/",filelist[2], sep=""))
  #otherwise wrong order
  runInfo<-runInfo[order(runInfo$file),]
  
  #dont want param files
  filelist2<-filelist[c(-1,-2)]




    someData <- rep(NaN, length(filelist2)* length(alpha)* length(pcrits))  
    someData2 <- rep(NaN, length(filelist2)*length(alpha)* length(pcrits)*2) 
    
    Ne<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    NeSv<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Ne2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    NeSv2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    trueNe1<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    trueNe2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    trueNe3<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    trueNe4<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    r2store<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Sved<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Vastore<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    
    Chijack<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Chijack2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Normjack<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Normjack2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    LDjack<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    LDjack2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    para<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    para2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
      
    ChijackSv<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    NormjackSv<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    LDjackSv<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    paraSv<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    ChijackSv2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    NormjackSv2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    LDjackSv2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    paraSv2<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    JPrime<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    JOrig<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    alphastore<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    pcritstore<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))
    Sestore<-array(someData, c(length(filelist2), length(pcrits), length(alpha)))

  
    #main estimation loop
    for(i in 1:length(filelist2)){
      #setTkProgressBar(pb, i/10, title = NULL, label = NULL)
      data<-read.fstat(paste(dirlist[y],"/", filelist2[i], sep=""),missing=NA,quiet=TRUE)
      
      #data<-read.fstat(testFile,missing=NA,quiet=TRUE)
      
      for(j in 1:length(pcrits)){
        
      x=jackfunC(data, pcrits[j], M)
      #return(c(r,J,sved, JKSe, JKVar, rvals))
      r<-x[1]
      J<-x[2]
      SvCorr<-x[3]
      Se<-x[4]
      Va<-x[5]
      #r2<-x[3]
      #J2<-x[4]
      
      S<-parameters$sample_size
      #S=dim(data@tab)[1]
      #Ne<-NehatWC(r,S)
      
     
      
      rvals<-x[-(1:5)]
      npr<-nprime(rvals)
           
      
      
      
      
      for(k in 1:length(alpha)){
        
        alp<-alpha[k]
        JPrime[i,j,k]<-npr
        JOrig[i,j,k]<-J
        Ne[i,j,k]<-NehatWC(r,S)
        NeSv[i,j,k]<-NehatSved(r,SvCorr)
        #Ne2[i,j,k]<-NehatWC(r2,S)
        #NeSv2[i,j,k]<-NehatSved(r2,SvCorr)
        
        list[Chijack[i,j,k],Chijack2[i,j,k]]<-ChijackknifeC(r,S,Va,alp)
        list[Normjack[i,j,k],Normjack2[i,j,k]]<-NormjackknifeC(r,S,Se,alp)
        
        list[LDjack[i,j,k],LDjack2[i,j,k]]<-LDjackknifeC(r,S,npr,alp)
        list[para[i,j,k],para2[i,j,k]]<-parametricC(r,S,J,alp)
        
        list[ChijackSv[i,j,k],ChijackSv2[i,j,k]]<-ChijackknifeSvC(r,S,Va,alp,SvCorr)
        list[NormjackSv[i,j,k],NormjackSv2[i,j,k]]<-NormjackknifeSvC(r,S,Se,alp,SvCorr)
        list[LDjackSv[i,j,k],LDjackSv2[i,j,k]]<-LDjackknifeSvC(r,S,npr,alp,SvCorr)
        list[paraSv[i,j,k],paraSv2[i,j,k]]<-parametricSvC(r,S,J,alp,SvCorr)
        
        alphastore[i,j,k]<-alp
        pcritstore[i,j,k]<-pcrits[j]
        trueNe1[i,j,k]<-runInfo$Ne1[i]
        trueNe2[i,j,k]<-runInfo$Ne2[i]
        trueNe3[i,j,k]<-runInfo$Ne3[i]
        trueNe4[i,j,k]<-runInfo$Ne4[i]
        r2store[i,j,k]<-r
        Sved[i,j,k]<-SvCorr
        Vastore[i,j,k]<-Va
        Sestore[i,j,k]<-Se
        
      }
      
      
      #
      
    }
  }


bigL<-list()
M<-length(someData)
for(i in 1:M){
  #w=c(1,2)
  w <- list( Cjup =Chijack2[i] , Cjdown=Chijack[i] , Njup = Normjack2[i], Njdown = Normjack[i], lJup =LDjack2[i] , Ljdown =LDjack[i] , pup =para2[i], pdown=para[i],  CjupSv =ChijackSv2[i] , CjdownSv=ChijackSv[i] , NjupSv = NormjackSv2[i], NjdownSv = NormjackSv[i], lJupSv =LDjackSv2[i] , LjdownSv =LDjackSv[i] , pupSv =paraSv2[i], pdownSv=paraSv[i] )
  v <- list(trueNe = trueNe1[i],trueNe2 = trueNe2[i],trueNe3 = trueNe3[i],trueNe4 = trueNe4[i], Ne=Ne[i], NeSv=NeSv[i], Ne2=Ne2[i], NeSv2=NeSv2[i],alpha=alphastore[i], pcrit=pcritstore[i], r2=r2store[i],Va=Vastore[i],Se=Sestore[i], Sved=Sved[i], Jp=JPrime[i], IndComp=JOrig[i])
  bigL[[length(bigL)+1]] <- list(c(v,w,parameters))
}

df <- data.frame(matrix(unlist(bigL), nrow=length(bigL), byrow=T))
names(df)<-c(names(v),names(w), names(parameters))

write.csv(df,paste(direOut,outputNames[y],".csv", sep=""))


}


