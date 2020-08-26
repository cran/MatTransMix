\name{print.object}
\alias{print.EM}
\alias{summary.EM}
\title{Functions for Printing or Summarizing Objects}
\description{
  \code{EM} classes for printing and summarizing objects.
}
\usage{
\method{print}{EM}(x, ...)
\method{summary}{EM}(object, ...)
}
\arguments{
  \item{x}{an object with the 'EM' class attributes.}
  \item{object}{an object with the 'EM' class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  Some useful functions for printing and summarizing results.
}
\examples{
set.seed(123)
data(IMDb)
Y <- IMDb$Y/100
p <- dim(Y)[1]
T <- dim(Y)[2]
n <- dim(Y)[3]
K <- 3
init <- MatTrans.init(Y, K = K, n.start = 2)
M <- MatTrans.EM(Y, initial = init, model = "X-VVV-VV", 
long.iter = 1000, silent = FALSE)
print.EM(M)
}
\seealso{
  \code{MatTrans.EM}.
}

