\name{MatTrans.EM}
\alias{MatTrans.EM}
% Also NEED an '\alias' for EACH other topic documented here.
\title{ EM algorithm for matrix clustering}
\description{
Runs the EM algorithm for matrix clustering}
\usage{
MatTrans.EM(Y, initial = NULL, la = NULL, nu = NULL, 
model = NULL, trans = "Gaussian", la.type = 0, tol = 1e-05, 
max.iter = 1000, size.control = 0, silent = TRUE)
}
% maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y }{ dataset of random matrices (p x T x n), n random matrices of dimensionality (p x T)}
  \item{initial }{ initialization parameters provided by function MatTrans.init()}
  \item{la }{ initial skewness for rows (K x p)}
  \item{nu }{ initial skewness for columns (K x T)}
   \item{model }{ parsimonious model type, if null, then all 196 models are run}
  \item{trans }{ transformation method: Gaussian, Power, Manly}
  \item{la.type }{ lambda type 0 or 1, 0: unrestricted, 1: same lambda across all variables}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
  \item{size.control }{ minimum size of clusters allowed for controlling spurious solutions}
  \item{silent }{ whether to produce output of steps or not}
}
\details{
Runs the EM algorithm for modeling and clustering matrices for a provided dataset. Both matrix Gaussian mixture, matrix Power mixture and matrix Manly transformation mixture can be employed. The user should use the MatTrans.init() function to get initial parameters and input them as 'initial'. In the case when transformation parameters are not provided but 'trans' is specified to be 'Power' or 'Manly', 'la' and 'nu' take value of 0.5. 'model' can be specified as 'X-XXX-XX'. The first digit 'X' stands for the mean structure. It is either 'G': general mean  or 'A': additive mean. The second 'XXX' specifies the variance-covariance Sigma. There are 14 options including EII, VII, EEI, VEI, EVI, VVI, EEE, EVE, VEE, VVE, EEV, VEV, EVV and VVV with detailed explanation as follows:
"EII" spherical, equal volume
"VII" spherical, unequal volume
"EEI" diagonal, equal volume and shape
"VEI" diagonal, varying volume, equal shape
"EVI" diagonal, equal volume, varying shape
"VVI" diagonal, varying volume and shape
"EEE" ellipsoidal, equal volume, shape, and orientation
"EVE" ellipsoidal, equal volume and orientation (*)
"VEE" ellipsoidal, equal shape and orientation (*)
"VVE" ellipsoidal, equal orientation (*)
"EEV" ellipsoidal, equal volume and equal shape
"VEV" ellipsoidal, equal shape
"EVV" ellipsoidal, equal volume (*)
"VVV" ellipsoidal, varying volume, shape, and orientation
 The last 2-digit 'XX' specifies the variance-covariance Psi. There are 7 options including II, EI, VI, EE, VE, EV, VV. The user can specify the 'model' to be for example 'X-VVV-EV', then both 'G' and 'A' mean structures will be fitted while Sigma and Psi are fixed at 'VVV' and 'EV', respectively. Similarly, 'model' can be specified as 'G-XXX-EV' or 'G-VVV-XX' for selection of Sigma and Psi structures.
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
Gauss <- MatTrans.EM(Y, initial = init, max.iter = 1000, model = "G-EII-UI")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }% __ONLY ONE__ keyword per line
