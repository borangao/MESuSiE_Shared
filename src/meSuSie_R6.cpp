#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::export]]
double loglik_cpp_R6(arma::vec V, const arma::mat& betahat, const arma::mat& shat2,const arma::vec& prior_weight, const int nancestry, arma::uvec diag_index){		
  mat Xty = betahat/shat2;
  
  
  mat cor_mat(nancestry,nancestry); 
  uvec upper_indices = trimatu_ind(size(cor_mat));
  cor_mat(upper_indices) = V;
  cor_mat = symmatu(cor_mat);
  cor_mat.diag().ones();
  
  mat se_mat = eye(nancestry,nancestry);
  V.elem(diag_index-1) =  sqrt(exp(V.elem(diag_index-1))) ;
  se_mat.diag() = V.elem(diag_index-1);
  
  arma::mat V_mat = se_mat*cor_mat*se_mat;
 // cout<<V_mat<<endl;
 // cout<<V_mat<<endl;
  vec lbf_multi(shat2.n_rows);
  double lbf_part_1,lbf_part_2;
  mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
  rowvec Xty_var;
  //cout<<shat2.n_rows<<endl;
  for(int i=0;i<shat2.n_rows;i++){
    SIGMA_INV = diagmat(shat2.row(i));
    SIGMA = inv(SIGMA_INV);
    V_SIGMA = V_mat*SIGMA;
    // post_var = diagmat(shat2.row(i))-inv_sympd(SIGMA+SIGMA*V_SIGMA);
    post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
    Xty_var = Xty.row(i)*post_var;	
    lbf_part_1 = 0.5*dot(Xty_var.t(),Xty.row(i));
    lbf_part_2 = 0.5*log(det(V_SIGMA+eye(nancestry,nancestry)));
    lbf_multi(i) = lbf_part_1 - lbf_part_2;
  }
  
  double maxlbf,weighted_sum_w;
  vec lbf_1,lbf_2,w_multi,w_1,w_2,max_vec;
  vec mu = zeros<vec>(shat2.n_rows);

    maxlbf = lbf_multi.max();
    w_multi = exp(lbf_multi - maxlbf)%prior_weight;
    weighted_sum_w = accu(w_multi);
  double lbf_model = maxlbf + log(weighted_sum_w);
  return(-1.0*lbf_model);
}
// [[Rcpp::export]]
SEXP mvlmm_reg(arma::mat betahat,arma::mat shat2, arma::mat V_mat){
  mat Xty = betahat/shat2;
  double nancestry = shat2.n_cols;
  vec lbf_multi(shat2.n_rows);
  mat post_mean_wmulti(shat2.n_rows,nancestry),post_mean2_wmulti(shat2.n_rows,nancestry);
  
  double lbf_part_1,lbf_part_2;
  mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
  rowvec Xty_var;
  
  for(int i=0;i<shat2.n_rows;i++){
    SIGMA_INV = diagmat(shat2.row(i));
    SIGMA = inv(SIGMA_INV);
    V_SIGMA = V_mat*SIGMA;
    // post_var = diagmat(shat2.row(i))-inv_sympd(SIGMA+SIGMA*V_SIGMA);
    post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
    Xty_var = Xty.row(i)*post_var;	
    lbf_part_1 = 0.5*dot(Xty_var.t(),Xty.row(i));
    lbf_part_2 = 0.5*log(det(V_SIGMA+eye(nancestry,nancestry)));
    lbf_multi(i) = lbf_part_1 - lbf_part_2;
    post_mean_wmulti.row(i) = Xty_var;
    post_mean2_wmulti.row(i)=square(Xty_var)+post_var.diag().t();
  }
  List res = List::create(Named("lbf") = lbf_multi , Named("post_mean")=post_mean_wmulti,Named("post_mean2")=post_mean2_wmulti);
  
  return(res);
  
}
/*Note that mu1,mu2 are cube of dimension (p, L, Nancestry) */

SEXP test_ELBO(arma::mat alpha,arma::cube mu1,arma::cube mu2,arma::cube XtX,arma::mat XtX_diag, arma::mat Xty,arma::vec yty,arma::vec N_vec ,arma::vec sigma2,double KL){
  
  arma::cube B_1 = mu1.each_slice()%alpha; //(p,L,nancestry)
  arma::cube B_2 = mu2.each_slice()%alpha; //(p,L,nancestry)
  arma::mat B_1_bar = sum(B_1,1); //Row sums of the elements returns p*nancestry matrix
  arma::mat B_2_bar = sum(B_2,1); //Row sums of the elements returns p*nancestry matrix
  
  vec BXXB(B_1.n_slices),BbarXXBbar(B_1.n_slices);
  for ( int i=0; i<B_1.n_slices; i++ ){
    BXXB(i) = accu(trans(B_1.slice(i))*XtX.slice(i)%trans(B_1.slice(i)));
    BbarXXBbar(i) =accu(trans(B_1_bar.col(i))*XtX.slice(i)*B_1_bar.col(i));
  }
  
  rowvec XXB2 = sum(XtX_diag%B_2_bar,0);
  rowvec BbarXty = 2*sum(B_1_bar%Xty,0);
  vec intermediate = yty - BbarXty.t()+BbarXXBbar+XXB2.t()-BXXB;
  vec sigma_update = intermediate/N_vec;
  double elbo =accu(-0.5*yty%log(2*datum::pi*sigma2) -0.5/sigma2%intermediate)-KL;
  List res = List::create(Named("ELBO") = elbo , Named("sigma2_update")=sigma_update);
  return(res);
}



