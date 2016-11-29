\name{optimization-package}
\alias{optimization-package}
\docType{package}
\title{
\packageTitle{optimization}
}
\description{
\packageDescription{optimization}
}
\details{
Multi-purpose optimization with Nelder-Mead and Simulated Annealing methods.
}
\author{
\packageAuthor{optimization}

Maintainer: \packageMaintainer{optimization}
}
\references{
Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by Simulated Annealing. Science, 220(4598): 671-680.

Nelder, J. and Mead, R. (1965). A simplex method for function minimization.
Computer Journal, 7(4).
}

\seealso{
\code{\link{optim_nm}}, \code{\link{optim_sa}}, \code{\link{optim}},  \code{\link{plot}}
}


\examples{
hi <- function(x){(x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2}
optim_nm(fun = hi, k = 2)
optim_sa(fun = hi, start = c(runif(2, min = -1, max = 1)),
  trace = FALSE,
  lower = c(-4, -4),
  upper = c(4, 4),
  control = list(dyn_rf = FALSE,
    rf = 1.2,
    t0 = 10,
    nlimit = 100,
    r = 0.6,
    t_min = 0.1
  )
)
}
