\name{optim_sa}
\alias{optim_sa}
\title{
Simulated Annealing Optimization
}
\description{
Numerical optimization algorithm that searches for the global optimum of non-linear functions. The function is allowed to be non-differentiable and multimodal.
}
\usage{
optim_sa(fun, start, maximization = FALSE, trace = FALSE,
         lower, upper, control = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun}{
  Function to be optimized. It should return a scalar value. The parameters should be assigned as a vector.
  }
  \item{start}{
  Vector of initial values for the function parameters. Must be of same length as the parameter vector of the function.
  }
  \item{maximization}{
  Logical. Has to be TRUE if the function is to be maximized. Default is FALSE.
  }
  \item{trace}{
  Logical. When TRUE, interim results are stored in a matrix. Necessary for the plot function. Default is FALSE.
  }
  \item{lower}{
  Vector of lower boundaries for the function parameters. Must be of same length as the parameter vector of the function and the initial values.
  }
  \item{upper}{
  Vector of upper boundaries for the function parameters. Must be of same length as the parameter vector of the function and the initial values.
  }
  \item{...}{
  Additional arguments defined by the user. May be useful for complex functions.
  }
  \item{control}{
  List with optional further arguments that can be used to control the behavior of the algorithm:
    \describe{
      \item{\code{vf}}{
        Function that determines the variation of the function parameters. Default is a uniform distributed random number.
      }
      \item{\code{rf}}{
        Numeric. Random factor that determines the range of the random number of \code{vf} in relation to the dimension of the function parameters. Default is 1. \code{rf} depends on the number of function calls which are out of the user specific bounds and the temperature. It is thus dynamic over time.
      }
      \item{\code{dyn_rf}}{
        Logical. It is sometimes useful to disable the dynamic \code{rf} changing. Especially when unimodal functions are solved.  As \code{rf} dynamics usually improve the performance as well as the accuracy it should be enabled when possible. Default is TRUE.
      }
      \item{\code{t0}}{
        Numeric. Initial temperature. Default is 10000.
      }
      \item{\code{nlimit}}{
        Integer. Maximum number of iterations of the inner loop. Default is 1000.
      }
      \item{\code{r}}{
        Numeric. Temperature reduction in the outer loop. Default is 0.6.
      }
      \item{\code{k}}{
        Numeric. Boltzmann constant in the Metropolis function. Default is 1.
      }
      \item{\code{t_min}}{
        Numeric. Limit temperature where the algorithm stops. Default is 0.01.
      }
      \item{\code{maxgood}}{
        Integer. Break criterion to improve the algorithm performance. Maximum number of function value improvements in the inner loop. Breaks the inner loop and starts the next. Default is 100.
      }
      \item{\code{stopac}}{
        Integer. Break criterion to improve the algorithm performance. Maximum number of repetitions where the function value improvement is less than a certain accuracy. Breaks the inner loop and starts the next. Default is 300.
      }
      \item{\code{ac_acc}}{
        Numeric. Accuracy of the \code{stopac} break criterion. Default is 1/10000 of the function value at initial parameter combination.
      }
    }
  }

}
\details{
Simulated Annealing is an optimization algorithm for solving complex functions that may have several minima (or maxima respectively) or are non-differentiable. The algorithms performance highly depends on the settings. It is therefore not recommended to use the Simulated Annealing Optimization Algorithm for non-complex functions. The control list is parameterized for solving a function of medium complexity. It should be changed when solving a relatively simple function (e. g. a three dimensional function) to improve the performance. For relatively complex functions the settings should be changed to improve the accuracy. The performance highly depends on \code{nlimit} and \code{r}. \code{t0} should increase with models complexity since the possibility of leaving a local optimum increases with the temperature.
}
\value{
The output is a nmsa_optim list object with following entries:
  \describe{
    \item{\code{par}}{
      Function parameters after optimization.
    }
    \item{\code{function_value}}{
      Function value after optimization.
    }
    \item{\code{trace}}{
      Matrix with interim results. NULL if \code{trace} was not activated.
    }
    \item{\code{fun}}{
      The optimized function.
    }
    \item{\code{start}}{
      The initial function parameters.
    }
    \item{\code{lower}}{
      The lower boundaries of the function parameters.
    }
    \item{\code{upper}}{
      The upper boundaries of the function parameters.
    }
    \item{\code{control}}{
      control arguments after optimization.
    }


  }
}

\references{
Corana, A., Marchesi, M., Martini, C. and Ridella, S. (1987), Minimizing Multimodal Functions of Continuous Variables with the 'Simulated Annealing' Algorithm. ACM Transactions on Mathematical Software, 13(3): 262-280.

Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by Simulated Annealing. Science, 220(4598): 671-680.
}
\author{
Kai Husmann
}

\seealso{
\code{\link{optim_nm}}, \code{\link{optim}}, \code{\link{plot}}
}

\examples{
# Himmelblau's function
hi <- function(x){(x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2}
# minimum at f(3, 2) = 0
 # f(-2.804, -3.131) = 0
 # f(-3.779, -3.283) = 0
 # f( 3.584, -1.848) = 0

# Random start values. Control arguments for relatively simple functions.
out <- optim_sa(fun = hi, start = c(runif(2, min = -1, max = 1)),
                trace = TRUE, lower = c(-4, -4) ,upper=c(4, 4),
                control = list(t0 = 1000, nlimit = 1500, r = 0.8))

# Visual inspection of algorithm behavior.
plot(out)
plot(out, type = 'contour')

ho <- function(x){
  x1 <- x[1]
  x2 <- x[2]

  fact1 <- sin(x1) * cos(x2)
  fact2 <- exp(abs(1 - sqrt(x1^2 + x2^2) / pi))
  y <- -abs(fact1 * fact2)
}

# Control arguments for functions with medium complexity.
x <- optim_sa(fun = ho, start = c(1, 1), lower = c(-10, -10),
              upper = c(10, 10), trace = TRUE)

co <- function(x){ # Colville function
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]

  term1 <- 100 * (x1^2 - x2)^2
  term2 <- (x1 - 1)^2
  term3 <- (x3-1)^2
  term4 <- 90 * (x3^2 - x4)^2
  term5 <- 10.1 * ((x2 - 1)^2 + (x4 - 1)^2)
  term6 <- 19.8 * (x2 - 1)*(x4-1)

  y <- term1 + term2 + term3 + term4 + term5 + term6
}

# Control arguments for relatively complex functions
x <- optim_sa(fun = co, start = c(4, 4, 4, 4),lower = c(-5, -5, -5, -5),
              trace = FALSE, upper = c(5, 5, 5, 5),
              control = list(r = 0.999, nlimit = 50))
}
