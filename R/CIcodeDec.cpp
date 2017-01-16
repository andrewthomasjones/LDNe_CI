#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <array>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixBase;
using Eigen::ArrayBase;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXi;
using Eigen::VectorXi;
using Eigen::Map; 
using Eigen::Ref;

using namespace std;


// [[Rcpp::export]]
ArrayXXi countTable(ArrayXXd data, int i, int j, int N, VectorXi loc, double pcrit) {
  //data=(data);
  int starti,startj, kj,ki;

  starti=loc.head(i).sum();
  startj=loc.head(j).sum();
  
  ki=loc(i);
  kj=loc(j);
  
  ArrayXXi datai2(N,ki);
  ArrayXXi dataj2(N,kj);

  ArrayXXi table(3,3);
  ArrayXXi bigTable(3,3*ki*kj);
  
    datai2=data.block(0,starti,N,ki).cast<int>();
    dataj2=data.block(0,startj,N,kj).cast<int>();

  ArrayXXd ps = datai2.colwise().sum().cast<double>()/N;
  ArrayXXd qs = dataj2.colwise().sum().cast<double>()/N;
//   
// //  cout<< N*(ps*(ps < pcrit).cast<double>()).sum() << endl;
// //  cout<< (qs < pcrit).sum()<< endl << endl;
//   
  for(int i=0; i< ki;i++){
    if(ps(i)<pcrit){
      datai2.col(i).fill(0);
    }
  }
  for(int i=0; i< kj;i++){
    if(qs(i)<pcrit){
      dataj2.col(i).fill(0);
    }
  }
//   

//   
//   
//   
//   
    int l = 0;
//   
  for(int a=0;a<ki;a++){
    for(int b=0;b<kj;b++){
       table.fill(0);
       for(int n=0;n<N;n++){
         //Rcpp::Rcout<<  datai2.col(a)(n) <<' ' << dataj2.col(b)(n) <<' ' << std::endl;
         
          table(datai2.col(a)(n),dataj2.col(b)(n))++;
       }
       Rcpp::Rcout << std::endl;
      // 
      bigTable.block(0,0+(l*3),3,3)=table;
      l++;

    }
  }
//   

  
  return bigTable;
  
}

double variance(ArrayXXd vec){
  double mean = (vec.sum())/(vec.size());
  ArrayXXd means(1,vec.size());
  means.fill(mean);
  return(((vec - means)*(vec - means)/vec.size()).cast<double>().sum());
}

// [[Rcpp::export]] 
ArrayXXi countTableJ(ArrayXXd data, int i, int j, int N, VectorXi loc, ArrayXXi jack) {
 
  
  
  data=(data)*2;
  int starti,startj, kj,ki, loc1, loc2;
  
  starti=loc.head(i).sum();
  startj=loc.head(j).sum();
  
  ki=loc(i);
  kj=loc(j);
  
  ArrayXXi datai2(N,ki);
  ArrayXXi dataj2(N,kj);
    
  datai2=data.block(0,starti,N,ki).cast<int>();
  dataj2=data.block(0,startj,N,kj).cast<int>();
  
  
  
  
  for(int a=0;a<ki;a++){
    for(int b=0;b<kj;b++){
      
      for(int n=0;n<N;n++){
          jack(n,starti+a)=datai2.col(a)(n);
          jack(n,startj+b)=dataj2.col(b)(n);
      }
     
    }
  }
  
  return jack;
  
}


//N is aample size, L is num loci, loc is allele counts,

// [[Rcpp::export]]
ArrayXXi pairsCombo(int L){
    
    ArrayXXi combo(2,(L*(L-1)/2));
    int p = 0;
    
    for(int i =0; i<L; i++){
        for(int j=0; j<i;j++){
          combo(0,p)=i;
          combo(1,p)=j;
          p++;
        }
          
    }
    return combo;
}

// [[Rcpp::export]]
int indCombCount(int L, VectorXi loc){
    int ind =0;
    
    for(int i=0; i<L;i++){
      for(int j=0; j<i;j++){
        ind+=(loc(i)-1)*(loc(j)-1);
      }
    }  
    return ind;
} 

// [[Rcpp::export]]
vector<int>   rIndex(int  n){
    
  std::vector<int> indexes;
  indexes.reserve(n);
  
  for (int i = 0; i < n; ++i){
      indexes.push_back(i);
  }
  std::random_shuffle(indexes.begin(), indexes.end());
  
  //std::cout << std::endl;
      
  return(indexes);    
}

// [[Rcpp::export]]
ArrayXXd rShuff(ArrayXXd data, VectorXi loc){
    //double tempArray[data.rows()];
    ArrayXXd data2(data.rows(),data.cols());
    data2.fill(0);
    //Eigen::Map< ArrayXXd > (tempArray, data.rows(),1) = data.col(0);
    //vector< double > tempForShuffle;(tempArray, tempArray + sizeof tempArray / sizeof tempArray[0]);
    vector<int> loc2;
    vector<int> indexes;
    
    for(int i=0; i< loc.size();i++){
      loc2.push_back(loc.head(i).sum());
    }
    
    
    for(int i=0; i<data.cols();i++){
      
      if(std::find(loc2.begin(), loc2.end(), i) != loc2.end()) {
             indexes = rIndex(data.rows());
              
      } 
      
    
      for (std::vector<int>::iterator it1 = indexes.begin(); it1 != indexes.end(); ++it1 ){
    
           data2(it1 - indexes.begin(),i) = data(indexes[*it1],i);
      }
    }
    return(data2);    
}

// [[Rcpp::export]]
vector<double> mainLoop(ArrayXXd data, int N, double pcrit, int L, VectorXi loc){
   int C, indComb,lengthAB,ki,kj, whichA, whichB;
   double p,q,hi,hj, Delta, DeltaAdj, rDelta, grandMean,grandMean2, ptot, qtot;
   
    C = (L*(L-1)/2); 
    ArrayXXi combinations(2,C);
    ArrayXXi AlleleComb(1,C);
    ArrayXXd rFinal(1,N);
    ArrayXXi nABList;
    ArrayXXd nAB(3,3);
    
    ArrayXXd rMeans(1,C);
    ArrayXXd DMeans(1,C);
    ArrayXXd actualPairs(1,C);
    
    combinations = pairsCombo(L);
    indComb=indCombCount(L,loc);
    
      // 
      // 
      // 
      // 
      // 
      // 
      // 
      // 
       for(int i =0; i<C; i++){
      //   
        whichA = combinations(0,i);
        whichB = combinations(1,i);
        
        //Rcpp::Rcout<<  N <<' ' <<  loc <<' ' << pcrit << std::endl;
        //Rcpp::Rcout<<  whichA  <<' ' <<  whichB << std::endl;
        
      nABList=countTable(data, whichA , whichB, N, loc, pcrit);
// //
      kj=loc(whichB);
      ki=loc(whichA);
// 
// 
      lengthAB = ki*kj;
      AlleleComb(i)=(ki-1)*(kj-1);
      ArrayXXd rhat(1,lengthAB);
      ArrayXXd rhat2(1,lengthAB);
      ArrayXXd effSamp(1,lengthAB);
      rhat.fill(0);
      rhat2.fill(0);

        double dropped =0.0;
        int pc=0;
        int qc=0;

        for(int a =0; a<lengthAB; a++){

          nAB=nABList.block(0,a*3,3,3).cast<double>();


          p = (2*nAB.row(0).sum()+nAB.row(1).sum())/(2*N);
          q = (2*nAB.col(0).sum()+nAB.col(1).sum())/(2*N);
          hj= nAB.row(0).sum()/(N);
          hi= nAB.col(0).sum()/(N);

          Delta = (2*nAB(0,0)+nAB(1,0)+nAB(0,1)+(nAB(1,1)/2))/(N)-(2*p*q);
          DeltaAdj = N*Delta/(N-1);
          rDelta = (DeltaAdj*DeltaAdj)/((p*(1-p)+(hj-p*p))*(q*(1-q)+(hi-q*q)));

          if((1-p)<pcrit||(1-q)<pcrit ){
            rDelta=0.0;
            dropped++;
          }

          if((1-p)<pcrit){
            pc++;
          }

          if((1-q)<pcrit ){
            qc++;
          }


          rhat(0,a)=rDelta;


        }

          rMeans(i) = rhat.sum()/(lengthAB-dropped);

          double part1 =0.0;
          double part2 =0.0;

          if(pc==0){
            part1 = (ki-1);
          }else{
            part1 = (ki-pc/kj);
          }
          if(qc==0){
            part2 = (kj-1);
          }else{
            part2 = (kj-qc/ki);
          }
          actualPairs(i) = part1*part2;
//           Rcpp::Rcout<<  " 4 " << std::endl;

        }
       
      vector<double> out;

      grandMean=((rMeans*actualPairs.cast<double>()).sum())/(actualPairs.cast<double>().sum());
      out.push_back(grandMean);
      out.push_back(actualPairs.cast<double>().sum());

      for(int i =0; i<C; i++){
        out.push_back(rMeans(i));
      }
      return(out);
}

//N is sample size, L is num loci, loc is allele counts,
// [[Rcpp::export]]
vector <double> testS(ArrayXXd data, int N, double pcrit, int L, VectorXi loc, int nPerm) {
  vector<double> output;
  Rcpp::Rcout<< N << " " << L << " " << pcrit << " " << nPerm << " " << std::endl;
  //cout<< loc << endl;
  
  
  

  Rcpp::Rcout<< "main again" << std::endl;
  
  //mainloop
  vector<double> mainCalc;//= mainLoop(data, N, pcrit,  L,  loc);
  
  Rcpp::Rcout<< "vec alloc" << std::endl;    
  //r
  output.push_back(mainCalc.at(0));
  
  return(output);
  
}


//N is sample size, L is num loci, loc is allele counts,
// [[Rcpp::export]]
vector <double> rBurrrowsS(ArrayXXd data, int N, double pcrit, int L, VectorXi loc, int nPerm) {
    vector<double> output;
    vector<double> perms;
    
   // Rcpp::Rcout<< N << " " << L << " " << pcrit << " " << nPerm << " " << std::endl;
    //cout<< loc << endl;
    
      
    //get sved correction
    for(int j =0; j<nPerm; j++){  
      ArrayXXd data2 = rShuff(data, loc);
      Rcpp::Rcout<< j << " " << std::endl;
      perms.push_back(mainLoop(data2, N, pcrit,  L,  loc).at(0));
     }
       
    double sum = std::accumulate(perms.begin(), perms.end(), 0.0);
    double mean = sum / perms.size();
        
    //Rcpp::Rcout<< "main again" << std::endl;    
    //mainloop
    vector<double> mainCalc= mainLoop(data, N, pcrit,  L,  loc);
    
    //Rcpp::Rcout<< "vec alloc" << std::endl;    
    //r
    output.push_back(mainCalc.at(0));
    
    //J
    output.push_back(mainCalc.at(1));
    
    //sved
    output.push_back(mean);
        
    //individual r vals
    for(int i=2;i<mainCalc.size();i++){ 
      output.push_back(mainCalc.at(i));
    }
    return(output);
}

