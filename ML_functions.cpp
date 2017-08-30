#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double GHeq_LL(double Z, NumericVector Lbar, NumericVector ss,
                 double Linf, double K, double Lc) {
  int count = Lbar.size();
  
  double Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));

  double sum_square_Lpred = 0.;
  double nyear = 0.;
  double sigma;
  double nLL;
  
  for(int m=0;m<count;m++) {
    if(ss[m]>0) {
      sum_square_Lpred += ss[m] * pow(Lbar[m]-Lpred, 2.);
      nyear += 1.;
    }
  }
  sigma = sqrt(sum_square_Lpred/nyear);
  nLL = -nyear * log(sigma) - 0.5 * sum_square_Lpred/(sigma*sigma);
  nLL *= -1;
  return nLL;
}



// [[Rcpp::export]]
double GHeffort(NumericVector stpar, NumericVector Lbar, NumericVector ss, 
                NumericVector eff, double Linf, double K, double a0, 
                double Lc, double eff_init, int n_age) {

  double ac = a0 - log(1 - Lc/Linf)/K;
  
  int n_yr = Lbar.size();
  int y;
  int a;
  double ndata = 0.;
  
  double q = stpar[0];
  double M = stpar[1];
  
  double Z_init = q * eff_init + M;
  NumericVector Z(n_yr);
  NumericMatrix N(n_yr,n_age);
  NumericVector La(n_age);
  NumericVector age(n_age);
  
  NumericVector Lpred(n_yr);
  double sum_square = 0.;
  double sigma;
  double nLL;
  
  for(a=0;a<n_age;a++) {
    double ageD = a;
    age[a] = ac + ageD;
    La[a] = Linf * (1 - exp(-K*(age[a] - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;
  
  N(0,0) = 1.;
  for(a=1;a<n_age;a++) N(0,a) = N(0,a-1) * exp(-Z_init);
  for(y=1;y<n_yr;y++) {
    N(y,0) = 1.;
    for(a=1;a<n_age;a++) N(y,a) = N(y-1,a-1) * exp(-Z(y-1));
    Lpred(y) = sum(N(y,_)*La)/sum(N(y,_));
    if(ss[y]>0) {
      ndata += 1.;
      sum_square += ss[y] * pow(Lbar[y] - Lpred[y], 2);
    }
  }
  sigma = pow(sum_square/ndata, 0.5);
  nLL = -ndata * log(sigma) - 0.5 * sum_square/(sigma*sigma);
  nLL *= -1;
  return nLL;
}
