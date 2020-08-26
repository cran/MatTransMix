########################
#### EMdemoIris #######
########################

set.seed(1)
data(iris)

data <- as.matrix(iris[,-5])
n <- 150
p <- 4
T <- 1

X <- array(NA, dim = c(p, T, n))
for(i in 1:n){
  X[,,i] <- data[i,]/10
}
K <- 3
init <- MatTrans.init(X, K = K, n.start = 10)

Trans <- MatTrans.EM(X, initial = init, model = "G-VVV-VV",
                    row.skew = TRUE, col.skew = TRUE,
                    trans = "Manly", silent = FALSE)





