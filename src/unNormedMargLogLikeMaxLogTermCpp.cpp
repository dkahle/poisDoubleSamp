#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double unNormedMargLogLikeMaxLogTermCpp(double phi, 
  int z11, int z12, int z21, int z22, 
  int m011, int m012, int m021, int m022, 
  int y01, int y02, 
  int N1, int N2, int N01, int N02
){

  
  double biggestLogTerm = Rf_lchoose(z12, 0) + Rf_lchoose(z22, 0) + 
    (m011 + z11 + 0)*log(phi) + 
    lgamma(m012 + z12 - 0 + 1) -
    (m012 + z12 - 0 + 1)*log(N1 + N01 + 0.) -
    (m022 + z22 - 0 + 1)*log(N2 + N02 + 0.) +
    lgamma(m022 + z22 - 0 + 1) +
    lgamma(m011 + m021 + z11 + z21 + 0 + 0 + 1) -
    (m011 + m021 + z11 + z21 + 0 + 0 + 1) * 
      log((N1 + N01)*phi + (N2 + N02)) +
    Rf_lbeta(y01 + 0 + 1, m011 - y01 + z11 + 1) +
    Rf_lbeta(y02 + 0 + 1, m021 - y02 + z21 + 1);
  
  
  // find the biggest variable to scale the others
  for(int k = 0; k < z12+1; ++k){
    for(int l = 0; l < z22+1; ++l){      
      biggestLogTerm = std::max(
        biggestLogTerm,
        Rf_lchoose(z12, k) + Rf_lchoose(z22, l) + 
        (m011 + z11 + k)*log(phi) + 
        lgamma(m012 + z12 - k + 1) -
        (m012 + z12 - k + 1)*log(N1 + N01 + 0.) -
        (m022 + z22 - l + 1)*log(N2 + N02 + 0.) +
        lgamma(m022 + z22 - l + 1) +
        lgamma(m011 + m021 + z11 + z21 + k + l + 1) -
        (m011 + m021 + z11 + z21 + k + l + 1) * 
          log((N1 + N01)*phi + (N2 + N02)) +
        Rf_lbeta(y01 + k + 1, m011 - y01 + z11 + 1) +
        Rf_lbeta(y02 + l + 1, m021 - y02 + z21 + 1)
      );      
    }
  }  
 
  return biggestLogTerm;
}
