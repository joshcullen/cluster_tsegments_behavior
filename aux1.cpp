#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int cat1(double value, NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector z(prob.nrow());

  for(int i=0; i<prob.nrow();i++){
    z[i]=cat1(randu[i],prob(i,_));
  }
  return z;
}

//' This function converts v into phi
// [[Rcpp::export]]
NumericVector GetPhi(NumericVector vec, int nclustmax) {
  NumericVector res(nclustmax);
  NumericVector prod1(nclustmax);
  double prod=1.0;
  
  for(int j=0; j<nclustmax;j++){
    res[j]=vec[j]*prod;
    prod=prod*(1-vec[j]);
  }
  
  return res;
}
