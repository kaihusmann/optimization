#--------------------------------------------------#
#### Simulated Annealing Optimization Algorithm ####
#--------------------------------------------------#

# Updated 15.01.2016 #

optim_sa <- function (fun, start, maximization = FALSE, trace = FALSE ,lower, upper, control = list(), ...) {
#   two <- 2
#   four <- .Call('optimization_testd', two ,PACKAGE = 'optimization')
#   print(four)

  #------------------#
  ## Initialisation ##
  #------------------#

  ## User declarations ##

  # start: Vector with starting values
  # fun: Function to be optimized
  # trace: Trace Matrix
  # lower, upper: Boundaries of the x variables
  # control: List with optional assignments


  ## Declaration of variables, vectors and matrices ##


  # Default control assignments:
  con <- list(
    vf      = NULL,
    rf      = 1,
    dyn_rf  = TRUE,
    t0      = 10000,
    nlimit  = 1000,
    r       = 0.9,
    k       = 1,
    t_min   = 0.01,
    maxgood = 100,
    stopac  = 300,
    ac_acc  = NA
  )
  corr_names <- names(con) # Saving the correct names for later consistency checks.

  con[names(control)] <- control # Overwriting the default assignments with user declarations.

  # Saving the list objects as stand alone variables.
  vf      <- con$vf
  rf      <- con$rf
  dyn_rf  <- con$dyn_rf
  t0      <- con$t0
  nlimit  <- con$nlimit
  r       <- con$r
  k       <- con$k
  t_min   <- con$t_min
  maxgood <- con$maxgood
  stopac  <- con$stopac
  ac_acc  <- con$ac_acc


  # Declaration of Intern variables:
  fun_length   <- length(start) # Number of parameters of loss function.
  temp         <- t0 # Starting temperature.
  para_0       <- start # Storing vector for the current parameter combination. Initialy with starting values.
  para_opt     <- start # Storing vector for current best parameter combination. Initialy with starting values.
  para_i       <- rep(NA, fun_length) # Vector with length = number of parameters
  loss_i       <- NA # Result of the loss function
  loss_0       <- fun(start) # Result of the loss function at initially the parameter comb.
  delta        <- NA
  goodcounter  <- 0
  ac           <- 0
  n_outer      <- 0 # Counter vector for the outer while loop
  n_inner      <- 0 # Counter vector for the inner loop (necessary since the inner loop has stop criteria and is not always of length 'nlimit').
  n_oob        <- rep(0, fun_length) # Counter vector for parameters out of bounds after changing function.
  ratio_noob   <- rep(0, fun_length) # Ratio between iterations of parameters NOT out of bounds and total number of iterations in the inner loop.

  ifelse(maximization, loss_opt <- -Inf, loss_opt <- +Inf)
  if (trace) {
    trace_array <- matrix(nrow = 0, ncol = 5 + fun_length, dimnames = list(character (0),
                                                                           c("iteration_outer",
                                                                             "function_value",
                                                                             paste("x", c(1 : fun_length), sep="_"),
                                                                             "n_iterations_inner",
                                                                             "temperature",
                                                                             "good_counter")
                                                                           )
                          ) # Saving matrix for the trace.
  } else {
    trace_array <- NULL # Empty variable. Necessary because the output is expected to have the same dimension even if trace is FALSE.
  }


  ## Check for consistency of user declarations ##

  # Check if there are start values and boundaries.
  if (!exists("start") &&  (length(lower) == 0 || length(upper) == 0)) {
    stop ("Starting values or boundaries are not defined.")
  }

  if (length(lower) != length(upper) || length(start) != length(upper) || length(start) != length(lower)) {
    stop ("Starting values or boundaries vector do not have the same length.")
  }

  if (any(start < lower | start > upper)) {
    start[start < lower | start > upper] <- apply(cbind(lower[start < lower | start > upper], upper[start < lower | start > upper]), 1, mean)
    warning ("At least one starting value was out of bounds. These were set to mean of their boundaries.", call. = FALSE)
    }

  if (mode(vf) == "function") {
    var_func <- vf
  } else {
    var_func <- function (para_0, fun_length, rf) {
      ret_var_func <- para_0 + runif(fun_length, 0.000001, rf) *  ((rbinom(fun_length, 1, 0.5) * -2) + 1)
      return (ret_var_func)
      }
  }

  if (r >=1) {
    r <- 0.9
    warning("r must be < 1. It is set to 0.9", call. = FALSE)
  }

  # The user can choose wheather rf is a scalar or a vector. A scalar will be extended (repeated) to the required length.
  # A vector of wrong length conditions a warning message and only the first entry is considered.
  if (length(rf) == 1) {
    rf <- rep(rf,fun_length)
  } else {
      if (!length (rf) == fun_length) {
        rf <- rep(rf[1], fun_length)
        warning ("rf was of wrong length. Only first value was conisidered.", call. = FALSE)
        }
  }

  # Are there wrong variable names in the control list?
  testnames <- names (control)
  if (length(testnames[!testnames %in% corr_names])!=0) {
    warning ("Control contains wrong names.", call. = FALSE)
    }

  # If ac_acc is not initialized by user, it will be stated relatively to the dimension of the initial y.
  if (is.na (ac_acc)) {ac_acc <- fun(start) / 10000}

  #----------------#
  ## Optimization ##
  #----------------#

  # Calling the Cpp source
  result <- .Call('optimization_main_loop',
              temp = temp,
              t_min = t_min,
              r = r,
              fun_length = fun_length,
              nlimit = nlimit,
              para_0 = para_0,
              para_i = para_i,
              var_func = var_func,
              rf = rf,
              lower = lower,
              upper = upper,
              fun = fun,
              loss_0 = loss_0,
              k = k,
              loss_opt = loss_opt,
              para_opt = para_opt,
              dyn_rf = dyn_rf,
              maxgood = maxgood,
              ac_acc,
              stopac = stopac,
              package = 'optimization')

#   if (! maximization) {
#     ## Minimization ##
#     while (temp > t_min) { # Outer while loop, number of repeatitions depends on cooling function.
#       goodcounter <- 0
#       n_outer     <- n_outer + 1
#       n_inner     <- 0
#       n_oob       <- rep(0, fun_length)
#
#       for (i in c(1 : nlimit)) { # Inner loop, no of repeatitions depends on the break criteria or on nlimit if no break criterion breaks the loop.
#
#         n_inner <- n_inner + 1 # Count the inner loop.
#         para_i <- var_func(para_0, fun_length, rf) # Variation of the parameters.
#         n_oob[para_i < lower | para_i > upper] <- n_oob[para_i < lower | para_i > upper] + 1 # Count the parameters which are out of bounds.
#         emergency_stop <- 0
#
#         while (any(para_i < lower) | any(para_i > upper)) { # Generate new values for these parameters.
#           emergency_stop <- emergency_stop + 1
#           para_i[para_i < lower | para_i > upper] <- var_func(para_0[para_i < lower | para_i > upper], length((para_i[para_i < lower | para_i > upper])), rf)
#           if (emergency_stop > 10000) {
#             stop ("The restrictions cannot be hold. Try different combination of starting values, boundaries or random factor.\n")
#             }
#         }
#
#         loss_i <- fun(para_i) # Result of the loss function at recent parameter comb.
#         delta  <- loss_i - loss_0
#
#         if (delta < 0) { # Comparison
#           loss_0 <- loss_i; para_0 <- para_i
#         } else { # This is the difference between Sim. Ann. and other Algorithms. It ist the prob. of accepting the worse loss.
#             if (runif(1) < exp(- (abs(delta) / (k * temp)))) {
#             para_0 <- para_i
#             loss_0 <- loss_i
#           }
#         }
#
#         if (loss_0 < loss_opt) {
#           goodcounter <- goodcounter + 1
#           loss_opt    <- loss_0
#           para_opt    <- para_0
#           savei       <- i; savet <- temp
#         }
#
#         if (goodcounter > maxgood) {print("bm"); break}
#         ifelse(abs(loss_0 - loss_opt) < ac_acc, ac <- ac + 1, ac <- 0) # Break statement for y-values oscillating around one eqilibrium.
#         #print(ac)
#         if (ac >= stopac) {print("ba");break}
#       } # End of the inner loop.
#
#       temp <- temp * r # Temperature reduction
#
#       if (trace) {
#         trace_array <- rbind(trace_array, c(n_outer,  loss_i, para_i, n_inner, temp, goodcounter))
#       }
#
#       if (dyn_rf) {
#         ratio_noob <- (n_inner-n_oob) / rep(n_inner, fun_length)
#
#
#         # Calculation of rf for the next iteration step according to the ratio of random values out of bounds (Corana et al. 1987)
#         # (assignment at end of the loop for the next step because necessary informations are lacking in the 1st iteration)
#         rf <- rf * ifelse(
#           ratio_noob >= 0.4 & ratio_noob <= 0.6,
#           1,
#           ifelse(
#             ratio_noob < 0.4,
#             1 / (1 + 2 * ((0.4 - ratio_noob) / 0.4)),
#             1 + (2 * ((ratio_noob - 0.6) / 0.4))
#             )
#           )
#
#         # Downscaling of the rf makes the algorithm more precise and efficient (Pronzato et al. 1984)
#         ifelse(n_outer <= 5, ds <- 1, ds <- 1 / (n_outer / 5))
#         ifelse(any(rf * ds <= 0.1), rf[rf * ds <= 0.1] <- 0.1, rf <- rf * ds)
#       }
#
#     }
#   } else {
#     ## Maximization ##
#
#     while (temp > t_min) { # Outer while loop, number of repeatitions  depends on cooling function.
#
#       goodcounter <- 0
#       n_outer     <- n_outer + 1
#       n_inner     <- 0
#       n_oob       <- rep(0, fun_length)
#
#       for (i in c(1 : nlimit)) { # Inner loop, no of repeatitions depends on the break criteria or on nlimit if no break criterion breaks the loop.
#
#         n_inner <- n_inner + 1 # Count the inner loop.
#         para_i <- var_func(para_0, fun_length, rf) # Variation of the parameters.
#         n_oob[para_i < lower | para_i > upper] <- n_oob[para_i < lower | para_i > upper] + 1 # Count the parameters which are out of bounds.
#         emergency_stop <- 0
#
#         while (any(para_i < lower) | any(para_i > upper)) { # Generate new values for these parameters.
#           emergency_stop <- emergency_stop + 1
#           para_i[para_i < lower | para_i > upper] <- var_func(para_0[para_i < lower | para_i > upper], length((para_i[para_i < lower | para_i > upper])), rf)
#           if (emergency_stop > 10000) {
#             stop ("The restrictions cannot be hold. Try different combination of starting values, boundaries or random factor.\n")
#             }
#         }
#
#         loss_i <- fun(para_i) # Result of the loss function at recent parameter comb.
#         delta  <- loss_i - loss_0
#
#         if (delta > 0) { # Comparison
#           loss_0 <- loss_i; para_0 <- para_i
#         } else { # This is the difference between Sim. Ann. and other Algorithms. It ist the prob. of accepting the worse loss.
#           if (runif(1) < exp(- (abs(delta) / (k * temp)))) {
#             para_0 <- para_i
#             loss_0 <- loss_i
#           }
#         }
#
#         if (loss_0 > loss_opt) {
#           goodcounter <- goodcounter + 1
#           loss_opt    <- loss_0
#           para_opt    <- para_0
#           savei       <- i; savet <- temp
#         }
#
#         if (goodcounter > maxgood) {maxgood_save <- rbind(maxgood_save, n_outer); break}
#         ifelse(abs(loss_0 - loss_opt) < ac_acc, ac <- ac + 1, ac <- 0) # Break statement for y-values oscillating around one eqilibrium.
#         if (ac >= stopac) {break}
#       } # End of the inner loop.
#
#       temp <- temp * r # Temperature reduction
#
#
#       if (trace) {
#         trace_array <- rbind(trace_array, c(n_outer,  loss_i, para_i, n_inner, temp, goodcounter))
#       }
#
#
#       if (dyn_rf) {
#         ratio_noob <- (n_inner-n_oob) / rep(n_inner, fun_length)
#
#         # Calculation of rf for the next iteration step according to the ratio of random values out of bounds (Corana et al. 1987)
#         # (assignment at end of the loop for the next step because necessary informations are lacking in the 1st iteration)
#         rf <- rf * ifelse(ratio_noob >= 0.4 & ratio_noob <= 0.6,
#                           1,
#                           ifelse(ratio_noob < 0.4,
#                                  1 / (1 + 2 * ((0.4 - ratio_noob) / 0.4)),
#                                  1 + (2 * ((ratio_noob - 0.6) / 0.4))
#                                  )
#                           )
#         # Downscaling of the rf makes the algorithm more precise and efficient (Pronzato et al. 1984)
#         ifelse(n_outer <= 5, ds <- 1, ds <- 1 / (n_outer / 5))
#         ifelse(any(rf * ds <= 0.1), rf[rf * ds <= 0.1] <- 0.1, rf <- rf * ds)
#       }
#
#     }
#   }


  #

  #----------------------------------------#
  ## Postprocessing and output generation ##
  #----------------------------------------#

  status_message <- tryCatch (data.frame(iteration = result[["savei"]], function_value = round(result[["loss_opt"]], 5), t(round(result[["para_opt"]], 5)), temperature = round(result[["savet"]], 5)),
    error = function (e) {
      warning ("Algorithm did not converge. Try different random factor or staring values.", call. = FALSE)
      return (data.frame(iteration = NA, temp = NA, y = NA, t(rep(NA, fun_length))))
    }
    )

  if (!is.na(status_message$iteration)) {
    cat("Algorithm converged.\n")
    names(status_message)[3 : (2 + fun_length)] <- paste("x", c(1 : fun_length), sep = "_")
    print(data.frame(status_message, row.names = ""),digits = 3)
  }

  output <- list(par            = result[["para_opt"]],
                 function_value = result[["loss_opt"]],
                 trace          = result[["trace_array"]],
                 fun            = fun,
                 start          = start,
                 lower          = lower,
                 upper          = upper,
                 control        = list(
                   rf      = rf,
                   t0      = t0,
                   nlimit  = nlimit,
                   r       = r,
                   k       = k,
                   t_min   = t_min,
                   maxgood = maxgood,
                   stopac  = stopac,
                   ac_acc  = ac_acc
                   )
                 )

  class(output) <- append(class(output), "optim_nmsa")
  return (output)
}
