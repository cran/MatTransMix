
data("crime")
X <- crime$Y[c(2,7),,] / 1000

p <- dim(X)[1]
T <- dim(X)[2]
n <- dim(X)[3]
K <-  2

initial <- MatTrans.init(X, K = K, n.start = 3)

Manly <- MatTrans.EM(X, initial = initial, la.type = 1, trans = "Manly", la = matrix(0.1, K, p), nu = matrix(0.1, K, T), max.iter = 1000, model = "A-VVV-VV")

Power <- MatTrans.EM(X, initial = initial, la.type = 0, la = matrix(1.1, K, p), nu = matrix(1.1, K, T), max.iter = 1000, model = "A-EVI-UI")

Gauss <- MatTrans.EM(X, initial = initial, la.type = 0, la = matrix(1, K, p), nu = matrix(1, K, T), max.iter = 1000, model = "G-VII-EI")


print(Manly$bic)

print(Power$bic)

print(Gauss$bic)
