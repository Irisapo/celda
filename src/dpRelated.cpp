// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector countwFixedL(IntegerVector labelVec, int Length) {
  IntegerVector Vcount(Length, 0);
  int L = labelVec.length();
  for (int i=0; i<L; ++i){
     Vcount[ labelVec[i]-1 ] += 1;
  }

  return(Vcount);
} 


// [[Rcpp::export]]
IntegerMatrix countNFixedCol(IntegerMatrix M, IntegerVector labelVec,int nG, int nCol, int RnCol) {
  IntegerMatrix countN(nG, RnCol);
  int label=0;
  for(int c=0; c<nCol; ++c) { 
    label = labelVec[c] - 1;
    for(int g=0; g<nG; ++g) {
      countN(g, label) += M(g,c); //operator () to index
    }
  }
  return(countN);
}



// [[Rcpp::export]]
IntegerVector reverseCumSum(IntegerVector v) {
  IntegerVector vrcount(v.size(), 0);
  int c = 0;
  for(int i = v.size(); i>0; --i) {
    vrcount[i-1] += c;
    c += v[i-1];
  }
  return(vrcount);
}


// [[Rcpp::export]]
NumericVector SBP2Prop(NumericVector a, NumericVector b, double alpha) {
  double sum_digamma_b = 0;
  double sum_digamma_ab = 0;
  int K = a.size();
	int i = 0;
  NumericVector logprob(K, 0.0);

	for (i = 0; i < K; i++) {
    // Add prior pseudo account to alpha
		a[i] += 1;
		b[i] += alpha;
	}
	NumericVector digamma_a = digamma(a);
	NumericVector digamma_b = digamma(b);
	NumericVector digamma_ab = digamma(a + b);

  for (i = 0; i < K; i++) {
    sum_digamma_ab += digamma_ab[i];
    logprob[i] = digamma_a[i] + sum_digamma_b - sum_digamma_ab;
    sum_digamma_b += digamma_b[i];
  }
  return(logprob);
}
