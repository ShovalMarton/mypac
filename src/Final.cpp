// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

//' @export
// [[Rcpp::export(kernel)]]
Rcpp::NumericVector kernel(Rcpp::NumericVector x, int krnl_type){
Rcpp::NumericVector W;

if (krnl_type==1) {;
        W = Rcpp::dnorm( x, 0.0,1.0);//Gause

} else if (krnl_type==2) {;
    Rcpp::NumericVector O;
    O = Rcpp::ifelse(Rcpp::abs(x) < 1,1,0);
    W = 0.75*(1-Rcpp::pow(x, 2)) * O ; //Epan

} else if (krnl_type==3) {;
    W =  ( Rcpp::abs ( x ) <= 1); //box

}else if (krnl_type==4) {;
  W =  ( 1/(Rcpp::exp(x)+2+Rcpp::exp(-x))); //Logistic
  
}else if (krnl_type==5) {;
  W =  (1/(Rcpp::exp(x)+2+Rcpp::exp(-x)))*(2/3.141592653589793238463); //Sigmoid function
} ;
/*int W;

//if(krnl_type == 1) {W = 1 ;

} else if (krnl_type == 0){ ;

W = 0;

};
//gause



*/
return W;

}

//' @export
// [[Rcpp::export(getBeta0)]]
double getBeta0(const arma::vec X, const arma::vec Y, const arma::vec w){
    // this function takes two numeric vectors
    // checks to make sure the dimensions match
    // returns the coefficients from a bivariate regression
    // TODO: add an intercept switch ... noInt

    int n = X.size();
    double x_w_sum = arma::as_scalar(w.t()*X);
    double y_w_sum = arma::as_scalar(w.t()*Y);

    arma::vec w_x_squ_vec(n);
    arma::vec w_y_x_vec(n) ;

    for(int i; i<n; i++){
      w_x_squ_vec(i) = X(i)*X(i)*w(i);
       w_y_x_vec(i) = X(i)*Y(i)*w(i);
    }

    double w_x_squ_sum = arma::accu(w_x_squ_vec);
    double w_y_x_sum = arma::accu(w_y_x_vec);
    double w_sum = arma::accu(w);

    double num  = w_x_squ_sum*y_w_sum - x_w_sum*w_y_x_sum;
    double dom  = w_sum*w_x_squ_sum - x_w_sum*x_w_sum;
    // find coefficients
    //double ans;

    //arma::mat ans = arma::zeros(1,1);
    double ans = num/dom;
    return ans;
}

//' @export
// [[Rcpp::export(ifelseYt)]]
arma::vec ifelseYt(arma::vec x,double num){

int n;
n = x.size();
bool state;
arma::vec x_ret (n);

for(int i = 0; i < n ; i++){
 state =  x(i)>=num ;
 if(state){x_ret(i) = 1;

 }else{
x_ret(i) = 0;
 }
}

return x_ret;
}

//' @export
// [[Rcpp::export(ifelseNt)]]
arma::vec ifelseNt(arma::vec x,arma::vec delta1,double num){

int n;
n = x.size();
bool state1, state2;
arma::vec x_ret (n);

for(int i = 0; i < n ; i++){
 state1 =  x(i)==num ;
 state2 =  delta1(i)==1 ;
 if(state1 && state2 ){x_ret(i) = 1;

 }else{
x_ret(i) = 0;
 }
}

return x_ret;
}

//' @export
// [[Rcpp::export(bw_RT_bool)]]
double bw_RT_bool(arma::vec x,arma::vec y,int krnl, double deg){

    double C_vp,p,bw,denom,numer,denom_numer;
    arma::mat x1, beta ,deriv,res,B,A,C;
    p = 1/(2 * deg + 3);

    for(int ii = 0 ; ii <= deg+3 ; ii++ ){
            x1 = arma::join_horiz(x1,pow(x,ii));

     }


   // W = arma::diagmat(weig);
    A = x1.t()* x1 ;
    B =  arma::inv(B,A);


//Finding the appropriate C_vp argument:

if (krnl==1) {;
        C_vp = 0.7763884;//Gause

} else if (krnl==2) {;
    C_vp = 1.719; //Epan

} else if (krnl==3) {;
    C_vp =  1.351 ; //Uniform
};
if(arma::as_scalar(B)==1){

    beta = (x1.t() * x1).i() * x1.t() * y ;

     deriv = 2*beta[2] + 6*beta[3]*x+12*beta[4]*pow(x,2);
      res = beta.t() * x1.t();
      res = y - res.t();

      denom = as_scalar(deriv.t() * deriv);
      numer = as_scalar(res.t()*res/x.size());
      denom_numer  = numer/denom;
      bw = C_vp*pow(denom_numer,p);

      }else{

          bw = R_NaN;
}


return bw;

}

//' @export
// [[Rcpp::export(bw_RT)]]
double bw_RT(arma::vec x,arma::vec y,int krnl, double deg){
  
  double bw;
  bw = bw_RT_bool(x = x, y = y, krnl = krnl, deg= deg);
  

  return bw;
  
}

//' @export
// [[Rcpp::export(est_h)]]
arma::mat est_h(arma::mat data, int krnl,arma::vec S_s, arma::vec T1_s ) {

  int LU_S,LU_T1,n_x;
  double bw;

  LU_S = S_s.size();
  LU_T1 = T1_s.size();

  arma::mat est_int1(LU_S,LU_T1);
  Rcpp::NumericVector x_var_we_s1;
  arma::vec y_var,x_var,Yvec_temp,Nvec,Yvec,Nvec_temp,weig,weights,V1,delta1,death_time;

  //est_int1 = arma::zeros(LU_S,LU_T1);

//v.fill(datum::nan);
  delta1 = data.col(6);
  V1 =  data.col(4);
  death_time = data.col(10);

  //Rcpp::NumericVector Yvec_temp,delta1,V1,Nvec,Yvec,death_time,x_varR;

  for(int ii = 0 ; ii < LU_S ; ii++) { // s
    for(int jj = 0; jj < LU_T1 ; jj++){ //  t1, t1 < s

        if(S_s(ii)<T1_s(jj)){
        continue;
        }

        Nvec_temp = ifelseNt(V1,delta1,T1_s(jj)) ;
        Yvec_temp = ifelseYt(V1 ,T1_s(jj));

        arma::uvec pos = find(Yvec_temp > 0);

        Yvec = Yvec_temp.elem(pos);
        Nvec = Nvec_temp.elem(pos);
        x_var = death_time.elem(pos)- S_s(ii);
        y_var = Nvec/Yvec ;

        n_x = x_var.size();

        //weig = arma::ones(n_x,1);

        bw = bw_RT_bool(x_var,y_var,krnl,1);


  // This is strange that the program not throw me out in a singular matrix case. It returns NaN instead. But R's results are identical.
        if(arma::is_finite(bw)){

                x_var_we_s1 = x_var/bw;
                weights = kernel(x_var_we_s1,krnl);
                est_int1(ii,jj) = getBeta0(x_var,y_var,weights);

                }else{;

                continue;

                }
            }
        }

  return est_int1;
}

//' @export
// [[Rcpp::export(Lambda_given_s)]]
arma::mat Lambda_given_s(arma::vec S_s,arma::vec T1_s,arma::mat est_int1) {
  
  int LU_S,LU_T1;
  arma::mat Surv_given_s;
  
  LU_S = S_s.size();
  LU_T1 = T1_s.size();
  
  Surv_given_s = arma::zeros(LU_S,LU_T1);
  
  //Rcpp::NumericVector Yvec_temp,delta1,V1,Nvec,Yvec,death_time,x_varR;
  
for(int ii = 0 ; ii < LU_S ; ii++) { // s
    
    for(int jj = 0; jj < LU_T1 ; jj++){ //  t1, t1 < s
      
      Surv_given_s(ii,jj) = arma::accu(est_int1(ii, arma::span(0, jj)));
      
      
    }
 }
  
  return  Surv_given_s;
}

//' @export
// [[Rcpp::export(final_CIF)]]
arma::mat final_CIF(arma::mat Surv_given_s,arma::vec T1_s,arma::vec dval) {

  arma::vec Final_est;
  int LU_T1;
  double col_jj;
  LU_T1 = T1_s.size();
  Final_est = arma::zeros(LU_T1);

            for(int jj = 0; jj < LU_T1 ; jj++){ //  t1, t1 < s
                    col_jj = arma::as_scalar((1-Surv_given_s.col(jj).t())*dval);

                    Final_est(jj) = col_jj;
                     }

  return  Final_est;
}

//' @export
// [[Rcpp::export(new_CIF1)]]
arma::mat new_CIF1(arma::mat data, int krnl,arma::vec S_s, arma::vec T1_s ) {
  
  
  int LU_S,LU_T1,n_x,mints,mintsFinal;
  double bw;
  
  LU_S = S_s.size();
  LU_T1 = T1_s.size();
  
  arma::mat est_int1(LU_S,LU_T1), Surv_given_s(LU_S,LU_T1) ;
  Rcpp::NumericVector x_var_we_s1;
  arma::vec y_var,x_var,Yvec_temp,Nvec,Yvec,Nvec_temp,weig,weights,V1,delta1,death_time;
  
  //est_int1 = arma::zeros(LU_S,LU_T1);
  
  //v.fill(datum::nan);
  delta1 = data.col(6);
  V1 =  data.col(4);
  death_time = data.col(10);
  
  //Rcpp::NumericVector Yvec_temp,delta1,V1,Nvec,Yvec,death_time,x_varR;
  
  for(int ii = 0 ; ii < LU_S ; ii++) { // s
    for(int jj = 0; jj < LU_T1 ; jj++){ //  t1, t1 < s
      
      if(S_s(ii)<T1_s(jj)){
        continue;
      }
      
      Nvec_temp = ifelseNt(V1,delta1,T1_s(jj)) ;
      Yvec_temp = ifelseYt(V1 ,T1_s(jj));
      
      arma::uvec pos = find(Yvec_temp > 0);
      
      Yvec = Yvec_temp.elem(pos);
      Nvec = Nvec_temp.elem(pos);
      x_var = death_time.elem(pos)- S_s(ii);
      y_var = Nvec/Yvec ;
      
      n_x = x_var.size();
      
      //weig = arma::ones(n_x,1);
      
      bw = bw_RT_bool(x_var,y_var,krnl,1);
      
      
      // This is strange that the program not throw me out in a singular matrix case. It returns NaN instead. But R's results are identical.
      if(arma::is_finite(bw)){
        
        x_var_we_s1 = x_var/bw;
        weights = kernel(x_var_we_s1,krnl);
        est_int1(ii,jj) = getBeta0(x_var,y_var,weights);
        
      }else{;
        
        continue;
        
      }
    }
  }
  
  arma::vec min_vec;
  
  min_vec = arma::zeros(2);
  
  //Rcpp::NumericVector Yvec_temp,delta1,V1,Nvec,Yvec,death_time,x_varR;
  
  for(int ii = 0 ; ii < LU_S ; ii++) { // s
    
    for(int jj = 0; jj < LU_T1 ; jj++){ //  t1, t1 < s
      
      min_vec(0) = S_s(ii);
      min_vec(1) = T1_s(jj);
      
      mints = min_vec.index_min();
      
      if(mints==0){mintsFinal = ii;
        
      }else{
        
        mintsFinal = jj;
        
      }
      
      Surv_given_s(ii,jj) = arma::sum(est_int1(ii, arma::span(0, mintsFinal)));
      
    }
  }
  Surv_given_s = arma::expmat(Surv_given_s);
  return  Surv_given_s;
  
}