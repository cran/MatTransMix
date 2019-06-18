\name{MatTrans.EM}
\alias{MatTrans.EM}
% Also NEED an '\alias' for EACH other topic documented here.
\title{ EM algorithm for matrix clustering}
\description{
Runs the EM algorithm for matrix clustering}
\usage{
MatTrans.EM(Y, initial = NULL, id = NULL, la = NULL, 
nu = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, 
model = NULL, trans = "Power", la.type = 0, tol = 1e-05, 
max.iter = 1000, size.control = 0, silent = TRUE)
}
% maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y }{ dataset of random matrices (p x T x n), n random matrices of dimensionality (p x T)}
  \item{initial }{ initialization parameters provided by function MatTrans.init()}
  \item{id }{ initial membership vector}
  \item{tau }{ initial mixing proportions (length K)}
  \item{la }{ initial skewness for rows (K x p)}
  \item{nu }{ initial skewness for columns (K x T)}
  \item{Mu }{ initial mean matrices (p x T x K)}
  \item{Sigma }{ initial array of sigma (p x p x K)}
  \item{Psi }{ initial array of psi (T x T x K)}
  \item{model }{ parsimonious model type, if null, then all 98 models are run}
  \item{trans }{ transformation method}
  \item{la.type }{ lambda type}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
  \item{size.control }{ minimum size of clusters allowed for controlling spurious solutions}
  \item{silent }{ whether to produce output of steps or not}
}
\details{
Runs the EM algorithm for modeling and clustering matrices for a provided dataset. Both matrix Gaussian mixture and matrix transformation mixture with given explanatory variables (data matrix X) or without explanatory variables (X is null) can be employed.  
A user has three options to initialize the EM algorithm. The user can use the MatTrans.init() function to get initial parameters and input them as 'initial'. The second  choice is to specify either initial id vector 'id' and transformation parameters 'la'. The third option is to input initial mode parameters 'la', 'tau', 'Mu', and 'Sigma' and 'Psi'. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a Gaussian mixture model. If some transformation parameters have to be excluded from the consideration, in the corresponding positions of matrix 'la', the user has to specify value 0.
A user also has three options to specify the covariance structure of the 'Psi' matrices, including unrestricted case, spherical matrices and autoregressive structured matrices. 
Notation: n - sample size, p x T - dimensionality of the random matrices, K - number of mixture components.
}
\value{
  \item{result }{ parsimonious models}
  \item{model }{ model types}
  \item{loglik }{ log likelihood values}
  \item{bic }{ bic values}
  \item{best.result }{ best parsimonious model}
  \item{best.model }{ best model type}
  \item{best.loglik }{ best logliklihood}
  \item{best.bic }{ best bic}
}

\examples{
set.seed(123)

data(crime)


Y <- crime$Y[c(2,7),,] / 1000

p <- dim(Y)[1]
T <- dim(Y)[2]
n <- dim(Y)[3]
K <-  2

init <- MatTrans.init(Y, K = K, n.start = 3)
Gauss <- MatTrans.EM(Y, initial = init, la = matrix(1, K, p), 
nu = matrix(1, K, T), max.iter = 1000, model = "G-EII-UI")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }% __ONLY ONE__ keyword per line
