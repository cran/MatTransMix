
data("crime")
X <- crime$Y[c(2,7),,] / 1000

p <- dim(X)[1]
T <- dim(X)[2]
n <- dim(X)[3]
K <-  2

initial <- MatTrans.init(X, K = K, n.start = 1)

Manly <- MatTrans.EM(X, initial = initial, la.type = 1, trans = "Manly", la = matrix(0.1, K, p), nu = matrix(0.1, K, T), max.iter = 1000, model = "A-VVV-VV")

Power <- MatTrans.EM(X, initial = initial, la.type = 0, trans = "Power", la = matrix(1.1, K, p), nu = matrix(1.1, K, T), max.iter = 1000, model = "A-EVI-UI")

Gauss <- MatTrans.EM(X, initial = initial, trans = "Gaussian", max.iter = 1000, model = "G-VII-EI")


print(Manly$bic)

print(Power$bic)

print(Gauss$bic)
