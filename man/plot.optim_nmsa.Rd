\name{plot.optim_nmsa}
\alias{plot}

\title{
Plot method for objects of class "optim_nmsa"
}
\description{
Creates convergence or contour plot for visual inspection of algorithm behavior.
}
\usage{
# S3 method for class 'optim_nmsa'
plot(x, type = 'convergence', lower = NA, upper = NA, ...)
}
\arguments{
  \item{x}{
  Object of type 'nmsa' to be plotted. Trace is expected to be TRUE.
  }
  \item{type}{
  Character string which determines the plot type. Either 'convergence' or 'contour' is allowed. 'convergence' is default.
  }
  \item{lower}{
  Vector containing the lower limits for the plotting window. Only usable with type = 'contour'.
  }
  \item{upper}{
  Vector containing the upper limits for the plotting window. Only usable with type = 'contour'.
  }
  \item{...}{
  Additional arguments defined by the user.
  }

}

\author{
Kai Husmann, Alexander Lange
}
\seealso{
\code{\link{optim_nm}}, \code{\link{optim_sa}}
}

\examples{
# Himmelblau's function
hi <- function(x){(x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2}

out_nm <- optim_nm(hi, k = 2, trace = TRUE)
out_sa <- optim_sa(fun = hi, start = c(runif(2, min = -1, max = 1)),
                   trace = TRUE, lower = c(-4, -4) ,upper=c(4, 4),
                   control = list(t0 = 1000, nlimit = 1500,r = 0.8))

# Visual inspection of algorithm behavior for Nelder-Mead.
plot(out_nm)
plot(out_nm, type = "contour", lower = c(-4, -4), upper = c(4, 4))

# Visual inspection of algorithm behavior for Simulated Annealing.
plot(out_sa)
plot(out_sa, type = "contour")
}
