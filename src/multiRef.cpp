#include <iostream>
#include <fstream>
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

// Function to update the P matrix for each reference;
// arma::mat getcomplexPart(mat Z,List PList,List SList, int r, int nCT){
//     arma::mat S1 = SList(0);
//     arma::mat totalSum = zeros<mat>(Z.n_rows,S1.n_rows);
//     int numRef = SList.size(); // number of references
//     for(int i = 0; i < numRef; i++){
//         if (i != r) {
//             arma::mat Qi = diagmat(Z.col(i));
//             arma::mat Pi_noNCT = PList(i);
//             arma::mat Si_noNCT = SList(i);
//             Si_noNCT.shed_col(nCT);
//             Pi_noNCT.shed_col(nCT);
//             //totalSum = totalSum + Qi * Pi * SitSi.col(nCT) - Qi * Pi.col(nCT) * diag_SitSi(nCT);
//             totalSum = totalSum + Qi * Pi_noNCT * Si_noNCT.t();

//         }
//     }
//     return(totalSum);
// }

arma::mat getWr(mat Z,List PList,List SList, int r, int mGene, int nSample){
  arma::mat totalSum = zeros<mat>(mGene, nSample);
  int numRef = SList.size(); // number of references
  for(int i = 0; i < numRef; i++){
    if (i != r) {
      arma::mat Qi = diagmat(Z.col(i));
      arma::mat Pi = PList(i);
      arma::mat Si = SList(i);
      totalSum = totalSum + Si * Pi.t() * Qi;
    }
  }
  return(totalSum);
}

arma::mat getMainDiff(mat Yinput, mat Z, List PList, List SList){
  arma::mat totalSum = Yinput;
  int numRef = SList.size(); // number of references
  for(int i = 0; i < numRef; i++){
    arma::mat Qi = diagmat(Z.col(i));
    arma::mat Pi = PList(i);
    arma::mat Si = SList(i);
    totalSum = totalSum - Si * Pi.t() * Qi;
  }

  return(totalSum);

}


// [[Rcpp::export]]
SEXP SWIMRMultiRef(SEXP YinputIn, SEXP SListIn, SEXP KIn, SEXP phi1In, SEXP phi2In, SEXP max_iterIn, SEXP epsilonIn, SEXP initPList, SEXP initsMatList, SEXP initSigma_e2, SEXP initLambdaList, SEXP initZ)
 {
   try {
     // read in the data
     arma::mat Yinput = as<mat>(YinputIn);
     List PList = Rcpp::clone(initPList);
     List SList = Rcpp::clone(SListIn);
     List sMatList = Rcpp::clone(initsMatList);
     int numRef = SList.size(); // number of references
     int nSample = (int)Yinput.n_cols; // number of spatial sample points
     int mGene = (int)Yinput.n_rows; // number of genes in spatial deconvolution
     arma::vec vecOne = ones<vec>( nSample);
     //arma::mat S = as<mat>(UIn);
     arma::mat K = as<mat>(KIn);
     arma::mat L = zeros<mat>(nSample,nSample);
     arma::mat D = zeros<mat>(nSample,nSample);
     arma::vec colsum_K = zeros<vec>(nSample);
     colsum_K = sum(K,1);
     D =  diagmat(colsum_K);// diagnol matrix whose entries are column
     double phi1 = as<double>(phi1In);
     double phi2 = as<double>(phi2In);
     arma::mat H1 = D -  phi1*K; // graph laplacian
     arma::mat H2 = D -  phi2*K; // graph laplacian
     double accu_H1 = accu(H1);
     double accu_H2 = accu(H2);
     int max_iter = Rcpp::as<int>(max_iterIn);
     double epsilon = as<double>(epsilonIn);
     //arma::mat P = as<mat>(initP);
     //arma::vec s = as<vec>(inits); // each reference should have a vector of s_rk
     double sigma_e2 = as<double>(initSigma_e2);
     List LambdaList = Rcpp::clone(initLambdaList); // in case the number of cell types are different
     //arma::mat lambda = as<mat>(initLambda);
     //arma::mat bMat = as<mat>(initbMat); // initialize the b matrix, which is R by K matrix
     //int k1 = (int)U.n_cols; // number of cell types in reference 1
     //int k2 = (int)U.n_cols; // number of cell types in reference 2
     // gamma prior for lambda_rk
     double alpha1 = 1.0;
     double beta1 = nSample / 2.0;
     // gamma prior for lambda_r
     double alpha2 = 1.0;
     double beta2 = nSample / 2.0;
     //double trac_xxt = accu(Yinput % Yinput);
     //initialize objective function;
     double obj = 0;
     double obj_old = 0;
     //PList_old = PList;
     arma::mat Z = as<mat>(initZ); // the weights for each reference per location
     bool logicalLogL = FALSE;
     // iteration starts
     for(int i = 1; i <= max_iter; ++i) {
       // cout << i <<endl;
       // update P_rk, the k-th cell type proportion corresponding to reference r
       double logLambda = 0.0;
       double logLZr = 0.0;
       double logLlambdar = 0.0;
       double logP = 0.0;
       for(int r = 0; r < numRef; r++){
         arma::vec Zr = Z.col(r);
         arma::mat Qr = diagmat(Zr);
         arma::mat Sr = SList(r);
         arma::mat Pr = PList(r);
         //arma::vec s_r = sMat.col(r);
         arma::mat Yt = Yinput.t();
         arma::mat QrYt = Yt.each_col() % Zr;
         arma::mat QrQrt = diagmat(Zr % Zr);
         // int k_r = Sr.n_cols;
         arma::mat SrtSr = Sr.t() * Sr;
         arma::vec diag_SrtSr = SrtSr.diag();
         arma::vec s_rk = sMatList(r);
         arma::vec lambda = LambdaList(r);
         // update Zr
         arma::mat PrSrt = Pr * Sr.t();
         arma::mat Urt = PrSrt * Yinput; // the diagonal the same as Ur as it is square matrix
         arma::mat Mr = PrSrt * PrSrt.t();
         arma::mat Wr = PrSrt * getWr(Z, PList, SList, r, mGene, nSample);
         // extract the diagonal lines
         arma::vec ur = arma::diagvec(Urt);
         arma::vec wr = arma::diagvec(Wr);
         arma::vec mr = arma::diagvec(Mr);
         // create Br
         arma::mat Br = diagmat(mr);
         // update sr
         arma::vec sr = sum(Zr.t() * H2, 1) / accu_H2;
         // update lambdar
         double temp_r = as_scalar((Zr.t() - sr * vecOne.t()) * H2 * (Zr - vecOne * sr.t()));
         double lambda_r = (temp_r / 2.0 + beta2 ) / (double(numRef) / 2.0 + alpha2 + 1.0);
         // update Zr
         arma::vec updateZr_num = lambda_r * ur + sigma_e2 * (phi2 * K * Zr + colsum_K * sr.t());
         arma::vec updateZr_den = lambda_r * (wr + Br*Zr) + sigma_e2 * (D * Zr + phi2 * colsum_K * sr.t());
         Zr = Zr % (updateZr_num / updateZr_den);
         Z.col(r) = Zr;
         // the  loglike.ihood for the first lambda
         logLambda = logLambda - (alpha1 + 1.0) * sum(log(lambda)) - sum(beta1 / lambda);
         // loglikeliohood for Zr
         temp_r = as_scalar((Zr.t() - sr * vecOne.t()) * H2 * (Zr - vecOne * sr.t()));
         logLZr = logLZr - (double)(nSample) * 0.5 * log(lambda_r) - 0.5 / lambda_r * temp_r;
         arma::mat temp_Pr = (Pr.t() - s_rk * vecOne.t()) * H1 * (Pr - vecOne * s_rk.t());
         logP = logP - (double)(nSample) * 0.5 * sum(log(lambda)) - 0.5 * (sum(temp_Pr.diag() / lambda));
         logLlambdar = logLlambdar - (alpha2 + 1.0) * log(lambda_r) - beta2 / lambda_r; // Should add a log sign to lambda_r, lambda_r is a double scalar, so logLlambdar = logLlambdar - (alpha2 + 1.0) * log(lambda_r) - beta2 / lambda_r;
       }
       // update sigma_e2
       arma::mat mainDiff = getMainDiff(Yinput, Z, PList, SList);
       double normNMF = accu(mainDiff % mainDiff);
       sigma_e2 = normNMF / (double)(mGene * nSample);

       // calculate likelihood function
       double logX = -(double)(nSample * mGene) * 0.5 * log(sigma_e2) - 0.5 * (double)(mGene * nSample); // Ying's correct: normNMF / sigma_e2 = mGene * nSample), so no need to change;
       obj = logX + logP + logLambda + logLZr + logLlambdar;
       logicalLogL = (obj > obj_old) && (abs(obj - obj_old) * 2.0 / abs(obj + obj_old) < epsilon);
       if(isnan(obj) || logicalLogL){
         if(i > 5){ // run at least 5 iterations
           break;
         }
       }else{
         obj_old = obj;
         //V_old = V;
       }
     }
     return List::create(Named("PList") = PList,
                         Named("Z") = Z,
                         Named("sigma_e2") = sigma_e2,
                         Named("lambda") = LambdaList,
                         Named("sMatList") = sMatList,
                         Named("Obj") = obj);
   }//end try
   catch (std::exception &ex)
   {
     forward_exception_to_r(ex);
   }
   catch (...)
   {
     ::Rf_error("C++ exception (unknown reason)...");
   }
   return R_NilValue;
 } // end funcs


