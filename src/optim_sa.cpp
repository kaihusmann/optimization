#include <Rcpp.h>
#include <math.h>
#include <iostream> // only needed for system output
using namespace std;  // only needed for system output
using namespace Rcpp;

// [[Rcpp::export]]
double cfun (double x, double y) {
  return ( ( x * x + y - 11.0 ) * ( x * x + y - 11.0 )
             + ( x + y * y - 7.0 ) * ( x + y * y - 7.0 ) );
}

NumericVector func (NumericVector para, Function fun) {
  NumericVector loss_i_temp = fun(para);
  return loss_i_temp;
}

// [[Rcpp::export]]
NumericVector var_funcc (NumericVector para_0, int fun_length, NumericVector rf) {
  NumericVector ret_var_func(fun_length);
  for(int k = 0; k < (fun_length); k++) {
    ret_var_func[k] = para_0[k] + R::runif(0.00000000001, rf[k]) * ((R::rbinom(1, 0.5) * -2) + 1);
  }

  // ret_var_func <- para_0 + runif(fun_length, 0.000001, rf) *  ((rbinom(fun_length, 1, 0.5) * -2) + 1)
  return ret_var_func;
}

// [[Rcpp::export]]
List main_loop (double temp, double t_min, double r, int fun_length, int nlimit, NumericVector para_0, NumericVector para_i, Function var_func, bool vf_user, NumericVector rf, NumericVector lower, NumericVector upper, Function fun, double loss_0, double k, double loss_opt, NumericVector para_opt, bool dyn_rf, double maxgood, double ac_acc, int stopac) {
  // Initializating variables outside the while loop
  IntegerVector n_oob(fun_length);
  int n_outer = 0;
  int savei = 0;
  double savet = 0;
  int ac = 0;


  // The outer while loop: Number of repeatitions depends on cooling function and the temp. limit.
  while (temp > t_min){


    // Initializing and resetting variables
    int goodcounter = 0;
    int n_inner = 0;
    n_outer++;

    std::fill(n_oob.begin(), n_oob.end(), 0);

    for (int i = 0; i < nlimit; i++) { // Inner loop, no. of repeatitions depends on the break criteria or on nlimit if no break criterion stops the loop.
      // Changing the parameters
      n_inner++;
      if(!vf_user){ // Variation of the parameters...
        para_i = var_funcc(para_0, fun_length, rf); // ...by the default function
      } else {
        para_i = var_func(para_0, fun_length, rf); // ...by a user declared function. This is an SEXP. The algorithm is therefore much slower with it.
      }

      // Counting the parameters which are out of bounds and change them.
      for(int j = 0; j < fun_length; j++){
        if(para_i[j] < lower[j] || para_i[j] > upper[j]){
          n_oob[j]++;
          // Generate new values for the variable until it is within the boundaries.
          int emergency_stop = 0;
          while (para_i[j] < lower[j] || para_i[j] > upper[j]) {
            emergency_stop++;
            NumericVector temp_para_i(1);

            if(!vf_user){ // Variation of the parameters.
              NumericVector para_0_j(1);
              para_0_j = para_0[j];
              NumericVector rf_j(1);
              rf_j = rf[j];
              temp_para_i = var_funcc(para_0_j, 1, rf_j); // By the default function
              //temp_para_i = var_func(para_0[j], 1, rf[j]);
            } else {
              temp_para_i = var_func(para_0[j], 1, rf[j]); // By a user declared function. This is an SEXP. The algorithm is therefore much slower with it.
            }
            // NumericVector temp_para_i = var_func(para_0[i], 1, rf[i]); // MUST BE UPDATED: C FUN NEEDED

            para_i[j] = temp_para_i[1];
            if (emergency_stop > 10000){stop("The restrictions cannot be hold. Try different combination of starting values, boundaries or random factor.");}
          }

        }

      }

      // Calculate the result of the loss function at recent parameter combination.

      // NumericVector loss_i_temp = fun(para_i); // Must be a vector for technical reasons (in rcpp)

      //double x = para_i[0]; double y = para_i[1];
      //NumericVector loss_i_temp = cfun(x, y);


      //NumericVector loss_i_temp = fun(para_i);
      NumericVector loss_i_temp = func(para_i, fun); // TBD: If simple function: C Function not SEXP

      double loss_i = loss_i_temp[0];
      double delta = loss_i - loss_0;
      // Check, if the loss has improved
      if (delta < 0){
        loss_0 = loss_i;
        para_0 = para_i;
      } else{ // This is the difference between Sim. Ann. and other Algorithms. It ist the prob. of accepting the worse loss.

        if (R::runif(0, 1) < exp (- fabs (delta) / (k * temp) )){
          loss_0 = loss_i;
          para_0 = para_i;
        }
      }
      if (loss_0 < loss_opt) {
        goodcounter++;
        loss_opt = loss_0;
        para_opt = para_0;
        savei = n_outer;
        savet = temp;
      }
      // Check for break criterions.
      if (goodcounter > maxgood) {break;}
      if (fabs(loss_0 - loss_opt) < ac_acc){
        ac++;

      }else{
        ac = 0;
      }
    if (ac > stopac){break;}
    } // End of the inner loop.


    temp = temp * r; // Temperature reduction.

    // Trace could be inserted here.

    // Calculation of rf for the next iteration step according to the ratio of random values out of bounds (Corana et al. 1987).
    if (dyn_rf == true){
      NumericVector ratio_noob (n_oob.size());
      for(int j = 0; j < n_oob.size(); j++){
        ratio_noob[j] = ( (double) n_inner - (double) n_oob[j]) / (double) n_inner;
        if (ratio_noob[j] < 0.4 || ratio_noob[j] > 0.6) {
          if (ratio_noob[j] < 0.4) {
            rf[j] = rf[j] * (1.0 / (1.0 + 2.0 * ((0.4 - (double) ratio_noob[j]) / 0.4)));
          }else{
            rf[j] = rf[j] * (1.0 + (2.0 * (( (double) ratio_noob[j] - 0.6) / 0.4)));
          }
        }


      }

      // Downscaling of the rf makes the algorithm more efficient (Pronzato et al. 1984)
      // Could be relative instead of 5
      NumericVector ds (n_oob.size());
      for(int j = 0; j < n_oob.size(); j++){
        if (n_outer <= 5) {
          ds[j] = 1.0;
       }else {
          ds[j] = 1.0 / ( (double) n_outer  / 5.0);
        }

        if (rf[j] * ds[j] <= 0.1) {
          rf[j] = 0.1;
        }else{
          rf[j] = rf[j] * ds[j];
        }

     }


    }



  }
  NumericMatrix trace_array (4, 4);
  List ret;
  ret["savei"] = savei;
  ret["savet"] = savet;
  ret["loss_opt"] = loss_opt;
  ret["para_opt"] = para_opt;

  return ret;
}
