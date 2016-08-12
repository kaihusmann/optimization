#include <Rcpp.h>
#include <math.h>
#include <iostream> // only needed for system output
using namespace std;  // only needed for system output
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector func (NumericVector para, Function fun) {
  NumericVector ret = fun(para);
  return ret[0];
}

// [[Rcpp::export]]
List main_loop (double temp, double t_min, double r, int fun_length, int nlimit, NumericVector para_0, NumericVector para_i, Function var_func, NumericVector rf, NumericVector lower, NumericVector upper, Function fun, double loss_0, double k, double loss_opt, NumericVector para_opt, bool dyn_rf, double maxgood, double ac_acc, int stopac) {
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

      n_inner++;
      para_i = var_func(para_0, fun_length, rf); // Variation of the parameters. Dieser Aufruf einer R Funktion verdoppelt die Laufzeit des Algorithmus, sollte noch irgendwie veraendert werden



      // Changing the parameters

      // Counting the parameters which are out of bounds
      for(int i = 0; i < n_oob.size(); i++){
        if(para_i[i] < lower[i] || para_i[i] > upper[i]){
          n_oob[i]++;
          // Generate new values for the variable until it is within the boundaries.
          int emergency_stop = 0;
          while (para_i[i] < lower[i] || para_i[i] > upper[i]) {
            emergency_stop++;

            NumericVector temp_para_i = var_func(para_0[i], 1, rf[i]);
            para_i[i] = temp_para_i[0];
            if (emergency_stop > 10000){stop("The restrictions cannot be hold. Try different combination of starting values, boundaries or random factor.");}
          }

        }

      }

      // Calculate the result of the loss function at recent parameter combination.
      //NumericVector loss_i_temp = fun(para_i);
      NumericVector loss_i_temp = func(para_i, fun);
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
      for(int i = 0; i < n_oob.size(); i++){
        ratio_noob[i] = ( (double) n_inner - (double) n_oob[i]) / (double) n_inner;
        if (ratio_noob[i] < 0.4 || ratio_noob[i] > 0.6) {
          if (ratio_noob[i] < 0.4) {
            rf[i] = rf[i] * (1.0 / (1.0 + 2.0 * ((0.4 - (double) ratio_noob[i]) / 0.4)));
          }else{
            rf[i] = rf[i] * (1.0 + (2.0 * (( (double) ratio_noob[i] - 0.6) / 0.4)));
          }
        }


      }

      // Downscaling of the rf makes the algorithm more precise and efficient (Pronzato et al. 1984)
      // Could be relativr instead of 5
      NumericVector ds (n_oob.size());
      for(int i = 0; i < n_oob.size(); i++){
        if (n_outer <= 5) {
          ds[i] = 1.0;
       }else {
          ds[i] = 1.0 / ( (double) n_outer  / 5.0);
        }

        if (rf[i] * ds[i] <= 0.1) {
          rf[i] = 0.1;
        }else{
          rf[i] = rf[i] * ds[i];
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
