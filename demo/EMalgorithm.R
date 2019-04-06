set.seed(123)


#Application to dataset IMDb
data("IMDb")
X <- IMDb$Y /100


n <- dim(X)[3]
p <- dim(X)[1]
T <- dim(X)[2]
K <- 2

#Run the EM algorithm with parsimonious models

initial <- MatTrans.init(X, K, n.start = 2, scale = 1)

Power <- MatTrans.EM(X, initial = initial, la.type = 0, la = matrix(0.9, K, p), nu = matrix(0.9, K, T), silent = FALSE)

Gauss <- MatTrans.EM(X, initial = initial, la.type = 0, la = matrix(1, K, p), nu = matrix(1, K, T), silent = FALSE)

print(Power$best.model)

print(Gauss$best.model)
