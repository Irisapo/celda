// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <Rcpp.h>
using namespace Rcpp;


IntegerVector sampleLl(const NumericVector& llProbs) {
    Function f(".sampleLl");
    return f(llProbs); // vector of length 1 and the index is 1 base 
}
    
//This is to replace the .cCCalcGibbsProbZ() function within the while loop in .celda_C 
// [[Rcpp::export]]
void celdaC_GibbsUpdate(const IntegerMatrix & counts, IntegerMatrix & mCPByS, IntegerMatrix & nGByCP, const IntegerVector & nByC, 
                    IntegerVector & nCP, IntegerVector & z, const IntegerVector & s, const int& K, const int& nG, const int& nM, 
                    const double& alpha, const double& beta, bool doSample = true) { 

    // iterate through each of the cells 
    int s_c = 0;
    int k_c = 0;
    for (int i=0; i<nM; ++i) { 
        NumericVector lprobK(K);  // to store log-probability of the K clusters for current cell 
        k_c = z[i] - 1 ;    // reindex k(c-cluster) into 0 base 
        s_c = s[i] - 1 ;    // reindex s(sample-index) into 0 base
        for (int k=0; k<K; ++k) { 
            if ( k != k_c ) { 
                lprobK[k] = log(mCPByS(k, s_c) + alpha);  // Theta simplified
                lprobK[k] += sum(lgamma(nGByCP( _, k) + counts( _, i) + beta)); //
                lprobK[k] -= lgamma(nCP[k] + nByC[i] + nG * beta); //
                lprobK[k] -= sum(lgamma(nGByCP( _, k) + beta)); //
                lprobK[k] += lgamma(nCP[k] + nG * beta); //
            } else {  // when k == k_c 
                lprobK[k] = log(mCPByS(k, s_c) - 1 + alpha);  // Theta simplified
                lprobK[k] += sum(lgamma(nGByCP( _, k) + beta)); //
                lprobK[k] -= lgamma(nCP[k] + nG * beta); //
                lprobK[k] -= sum(lgamma(nGByCP( _, k) - counts( _, i) + beta)); //
                lprobK[k] += lgamma(nCP[k] - nByC[i] + nG * beta); //
            }
        }
        if (doSample) { 
            int newk_c =0; 
            newk_c = sampleLl(lprobK)[0] - 1;  // IntegerVector to int, then reindex k(c-cluster) to 0 base 

            if (newk_c != k_c) { 
                mCPByS(k_c, s_c) -= 1; 
                mCPByS(newk_c, s_c) += 1; 
                
                nGByCP( _, k_c) = nGByCP( _, k_c) - counts( _, i); 
                nGByCP( _, newk_c) = nGByCP( _, newk_c) + counts( _, i);
                
                nCP[k_c] -= nByC[i];
                nCP[newk_c] += nByC[i];

                z[i] = newk_c + 1; // reindex k(c-cluster) into 1 base
            }
        }
    }

    // calcualte and return log-likelihood
}

// This is to calculate the probability of each cell in each of the k(c-cluster)
// Notice the log-probablity will never be the same as from calculated in the original .cCCalcGibbsProbZ()$probs with the same z, 
// this is because the math formula is not the same. This version is simplified (same terms have been calceled out). 
// [[Rcpp::export]]
NumericMatrix cC_calProb(const IntegerMatrix & counts, const IntegerMatrix & mCPByS, const IntegerMatrix & nGByCP, const IntegerVector & nByC,
                         const IntegerVector & nCP, const IntegerVector & z, const IntegerVector & s, const int& K,
                         const int& nG, const int& nM, const double& alpha, const double& beta) { 

    NumericMatrix lprobs(nM, K); 
    for (int i=0; i<nM; ++i) { 
        int k_c = z[i] -1;   // reindex to 0 base
        int s_c = s[i] -1;  // reindex to 0 base 
        for(int k=0; k<K; ++k) { 
            if (k != k_c) { 
                lprobs(i, k) = log(mCPByS(k, s_c) + alpha);  // Theta simplified
                lprobs(i, k) += sum(lgamma(nGByCP( _, k) + counts( _, i) + beta)); //
                lprobs(i, k) -= lgamma(nCP[k] + nByC[i] + nG * beta); //
                lprobs(i, k) -= sum(lgamma(nGByCP( _, k) + beta)); //
                lprobs(i, k) += lgamma(nCP[k] + nG * beta); //
            } else {
                lprobs(i, k) = log(mCPByS(k, s_c) - 1 + alpha);  // Theta simplified
                lprobs(i, k) += sum(lgamma(nGByCP( _, k) + beta)); //
                lprobs(i, k) -= lgamma(nCP[k] + nG * beta); //
                lprobs(i, k) -= sum(lgamma(nGByCP( _, k) - counts( _, i) + beta)); //
                lprobs(i, k) += lgamma(nCP[k] - nByC[i] + nG * beta); //           
            }
        }
    }
    return lprobs;
}

// [[Rcpp::export]]
NumericMatrix cC_calProbT(const IntegerMatrix & counts, const IntegerMatrix & mCPByS, const IntegerMatrix & nGByCP, const IntegerVector & nByC,
                         const IntegerVector & nCP, const IntegerVector & z, const IntegerVector & s, const int& K,
                         const int& nG, const int& nM, const double& alpha, const double& beta) { 

    NumericMatrix lprobs(K, nM); 
    for (int i=0; i<nM; ++i) { 
        int k_c = z[i] -1;   // reindex to 0 base
        int s_c = s[i] -1;  // reindex to 0 base 
        for(int k=0; k<K; ++k) { 
            if (k != k_c) { 
                lprobs(k, i) = log(mCPByS(k, s_c) + alpha);  // Theta simplified
                lprobs(k, i) += sum(lgamma(nGByCP( _, k) + counts( _, i) + beta)); //
                lprobs(k, i) -= lgamma(nCP[k] + nByC[i] + nG * beta); //
                lprobs(k, i) -= sum(lgamma(nGByCP( _, k) + beta)); //
                lprobs(k, i) += lgamma(nCP[k] + nG * beta); //
            } else {
                lprobs(k, i) = log(mCPByS(k, s_c) - 1 + alpha);  // Theta simplified
                lprobs(k, i) += sum(lgamma(nGByCP( _, k) + beta)); //
                lprobs(k, i) -= lgamma(nCP[k] + nG * beta); //
                lprobs(k, i) -= sum(lgamma(nGByCP( _, k) - counts( _, i) + beta)); //
                lprobs(k, i) += lgamma(nCP[k] - nByC[i] + nG * beta); //           
            }
        }
    }
    return lprobs;
}
     
// This is to replace .cCCalcLL() function
// [[Rcpp::export]]
double celdaC_llh(const IntegerMatrix & mCPByS, const IntegerMatrix & nGByCP, const IntegerVector & s, const IntegerVector & z,
                  const int& K, const int& nS, const int& nG, const double& alpha, const double& beta){ 
    
    // calculate for "Theta" component
    double a = nS * lgamma(K * alpha); 
    double b = 0.0; // in R:  b <- sum(lgamma(mCPByS + alpha))   
    double c = nS * K * lgamma(alpha); 
    double d = 0.0; 
    
    for (int s=0; s<nS; ++s) { 
        double temd = K * alpha;
        for (int k=0; k<K; ++k) { 
            b += lgamma(mCPByS(k, s) + alpha);   // or use Armadillo for vectorization --> need to check datatype conversion etc
            temd += mCPByS(k, s); 
        }
        d += lgamma(temd); 
    }

    double thetaLl = a + b - c - d; 

    // calculate for "Phi" component 
    a = K * lgamma(nG * beta); 
    b = 0.0; 
    c = K * nG * lgamma(beta); 
    d = 0.0; 

    for (int k=0; k<K; ++k) { 
        double temd = nG * beta; 
        for (int g=0; g<nG; ++g) {
            b += lgamma(nGByCP(g, k) + beta); 
            temd += nGByCP(g, k); 
        }
        d += lgamma(temd); 
    }

    double phiLl = a + b - c - d; 

    return thetaLl + phiLl; 
}


NumericMatrix fastMultMat(const NumericMatrix & phi, const IntegerMatrix & counts) { 

    const Eigen::Map<Eigen::MatrixXd> Phi(as<Eigen::Map<Eigen::MatrixXd>>(phi)); 
    const Eigen::Map<Eigen::MatrixXi> Counts(as<Eigen::Map<Eigen::MatrixXi>>(counts));

    Eigen::MatrixXd Probs = Phi.transpose() * Counts.cast<double>();
    
    return wrap(Probs); 
}

IntegerVector updateZ_EM(const IntegerMatrix & counts, const IntegerMatrix & mCPByS, const IntegerMatrix & nGByCP, 
                         const IntegerVector & s, const int& nM, const double& alpha, const double& beta) { 
    // calculate theta
    int rowSize = mCPByS.nrow(), colSize = mCPByS.ncol();  
    NumericMatrix theta(rowSize, colSize); 

    for (int i=0; i<colSize; ++i) { 
        double i_sum = sum(mCPByS( _, i)) + rowSize * alpha;
        if( i_sum <= 0 ) { 
            stop("Division by 0. Make sure colSums of counts does not contain 0 after rounding counts to integers."); 
        }
        for(int j=0; j<rowSize; ++j) {
            theta(j,i) = log((mCPByS(j, i) + alpha) / i_sum); 
        }
    }
    
    // calculate phi
    rowSize = nGByCP.nrow();
    colSize = nGByCP.ncol(); 
    NumericMatrix phi(rowSize, colSize); 

    for(int i=0; i<colSize; ++i){ 
        double i_sum = sum(nGByCP( _, i)) + rowSize * beta; 
        if ( i_sum <=0 ) { 
            stop("Division by 0. Make sure colSums of counts does not contain 0 after rounding counts to integers."); 
        }
        for(int j=0; j<rowSize; ++j) {
            phi(j,i) = log((nGByCP(j,i) + beta) / i_sum);
        }
    } 

    NumericMatrix probs = fastMultMat(phi, counts); 

    IntegerVector newZ(nM);
    for (int i=0; i<nM; ++i) {
        newZ[i] = 1 + which_max( probs( _, i) + theta( _, s[i] - 1));  // 1 base k(c-cluster)
    }

    return newZ; 
}

// [[Rcpp::export]]
void celdaC_EMUpdate(const IntegerMatrix & counts, IntegerMatrix & mCPByS, IntegerMatrix & nGByCP, const IntegerVector & nByC, 
                    IntegerVector & nCP, IntegerVector & z, const IntegerVector & s, const int& K, const int& nG, const int& nM, 
                    const double& alpha, const double& beta, bool doSample = true) { 
    
    // update z and others
    if (doSample) {
        IntegerVector previousZ = clone(z); 
        IntegerVector newZ = updateZ_EM(counts, mCPByS, nGByCP, s, nM, alpha, beta);
        // update z(in 1 base),  nGByCP, nCP and mCPByS
        for( int i=0; i<nM; ++i) {
            if (previousZ[i] != newZ[i]) {
                z[i] = newZ[i];
                nGByCP( _, previousZ[i] - 1) = nGByCP( _, previousZ[i] - 1) - counts( _, i);   // 1-base index to 0-base for slicing
                nGByCP( _, z[i] - 1 ) = nGByCP( _, z[i] -1 ) + counts( _, i);  
                nCP[previousZ[i] - 1] -= nByC[i];
                nCP[z[i] - 1] += nByC[i];
                mCPByS(previousZ[i] - 1, s[i] - 1) -= 1; 
                mCPByS(z[i] - 1, s[i] -1) += 1; 
            }
        }
    }
}


